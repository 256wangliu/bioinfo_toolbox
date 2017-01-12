#!/bin/bash

# Mapping (STAR) and transcript quantification (RSEM) pipeline.
# This is customized version of ENCODE long RNA-seq pipeline.
# The original shell script can be found from here:
# https://github.com/ENCODE-DCC/long-rna-seq-pipeline/

# Output:
# Aligned.sortedByCoord.out.bam                    alignments, standard sorted BAM, agreed upon formatting
# Signal.{Unique,UniqueMultiple}.unstranded.bw     2 bigWig files for unstranded data
# Signal.{Unique,UniqueMultiple}.{plus, minus}.bw  4 bigWig files for stranded data
# Quant.genes.results                              RSEM gene quantifications, tab separated text, RSEM formatting
# Quant.isoforms.results                           RSEM transcript quantifications, tab separated text, RSEM formatting
# Quant.pdf                                        RSEM diagnostic plots
# Log.final.out                                    mapping statistics to be used for QC, text, STAR formatting

#### Software path
export PATH=$PATH:/home/tsujij/soft/samtools-1.0/
export PATH=$PATH:/home/tsujij/soft/rsem-1.2.22/
export PATH=$PATH:/home/tsujij/soft/STAR_2.4.2a/source/
export PATH=$PATH:/home/tsujij/soft/ucsc_tools/  # bedGraphToBigWig & bigBedToBed

#### Usage
function usage {
cat <<EOF
Usage: $0 [options] -D <data-type> -a <star-dir> -b <rsem-dir/prefix> -f <fastq1> -r <fastq2>

Map reads with STAR and calculate transcript abundance with RSEM.

Required arguments:
  -D <data-type>         RNA-seq type:
                            - 'str_SE': stranded single-end RNA-seq
                            - 'str_PE': stranded paired-end RNA-seq
                            - 'unstr_SE': unstranded single-end RNA-seq
                            - 'unstr_PE': unstranded paired-end RNA-seq
  -a <star-dir>          STAR index. Only specify directory name.
  -b <rsem-dir/prefix>   RSEM index. Specify prefix name of index files.
  -f <fastq1>            FASTQ containing single-end reads or paired-end upstream mates.
  -r <fastq2>            Specify -1 if single-end. FASTQ containing paired-end downstream mates.
                         Both <fastq1> and <fastq2> should be gzipped.

Options:
  -h                     Show this message and exit.
  -p <cpu-num>           Number of CPUs (default: 8 cores).
  -o <output-dir>        Output directory (default: prefix of <fastq1>).
EOF
}

#### Parameters
NTHREADS=8
unset DTYPE
unset STAR_INDEX
unset RSEM_INDEX
unset FASTQ_1
unset FASTQ_2
unset OFPREFIX

#### Read arguments and options
[[ $# -eq 0 ]] && usage && exit 1;
while getopts "hD:a:b:f:r:o:p" OPT
do
  case ${OPT} in
    "D") DTYPE=${OPTARG}      ;;
    "a") STAR_INDEX=${OPTARG} ;;
    "b") RSEM_INDEX=${OPTARG} ;;
    "f") FASTQ_1=${OPTARG}    ;;
    "r") FASTQ_2=${OPTARG}    ;;
    "o") OFPREFIX=${OPTARG}   ;;
    "p") NTHREADS=${OPTARG}   ;;
    "h") usage && exit 1      ;;
     * ) usage && exit 1      ;;
  esac
done
[[ $# -lt 10 ]] && echo "$0: need more arguments" && exit 1;

#### Check arguments
if [ "${DTYPE}" != "str_SE" ] && [ "${DTYPE}" != "str_PE" ] && [ "${DTYPE}" != "unstr_SE" ] && [ "${DTYPE}" != "unstr_PE" ]
  then echo "$0: unknown data type: -D ${DTYPE}" && exit 1; fi
if [ ! -d "${STAR_INDEX}" ]
  then echo "$0: input STAR index" && exit 1; fi
if [ $( ls -1 ${STAR_INDEX}/* | wc -l ) -le 4 ]
  then echo "$0: can't interpret STAR index" && exit 1; fi
if [ -z "${RSEM_INDEX}" ]
  then echo "$0: input RSEM index" && exit 1; fi
if [ $( ls -1 ${RSEM_INDEX}.* | wc -l ) -le 4 ]
  then echo "$0: can't interpret RSEM index" && exit 1; fi
if [ "${DTYPE}" == "str_PE" ] || [ "${DTYPE}" == "unstr_PE" ]
then
  if [ ! -f "${FASTQ_1}" ] || [ ! -f "${FASTQ_2}" ]
    then echo "$0: can't read FASTQ" && exit 1; fi
  FASTQ_1=`readlink -f ${FASTQ_1}`
  FASTQ_2=`readlink -f ${FASTQ_2}`
else
  if [ ! -f "${FASTQ_1}" ]
    then echo "$0: can't read FASTQ" && exit 1; fi
  FASTQ_1=`readlink -f ${FASTQ_1}`
  FASTQ_2=" "
fi
if [ $( echo ${NTHREADS} | sed 's/^[+0-9][0-9]*//' | wc -c ) -ne 1 ]
  then echo "$0: can't interpret -p ${NTHREADS}" && exit 1; fi
if [ -z "${OFPREFIX}" ]
  then OFPREFIX=`echo ${FASTQ_1} | awk -F ".fq" '{print $1}' | awk -F ".fastq" '{print $1}'`; fi
STAR_INDEX=`readlink -f ${STAR_INDEX}`
RSEM_INDEX=`readlink -f ${RSEM_INDEX}`

#### Make work directory
mkdir ${OFPREFIX} && cd ${OFPREFIX}

#### Set STAR paramters
STAR_COMMON_PARAM=" --genomeDir ${STAR_INDEX}  --readFilesIn ${FASTQ_1} ${FASTQ_2} \
 --outSAMunmapped Within  --outFilterType BySJout  --outSAMattributes NH HI AS NM MD \
 --outFilterMultimapNmax 20  --outFilterMismatchNmax 999  --outFilterMismatchNoverReadLmax 0.04 \
 --alignIntronMin 20  --alignIntronMax 1000000  --alignMatesGapMax 1000000 \
 --alignSJoverhangMin 8  --alignSJDBoverhangMin 1  --sjdbScore 1  --readFilesCommand zcat"
STAR_RUN_PARAM=" --runThreadN ${NTHREADS}"
STAR_BAM_PARAM=" --quantMode TranscriptomeSAM --outSAMtype BAM SortedByCoordinate --outSAMheaderHD @HD VN:1.4 SO:coordinate"
unset STAR_STRAND_PARAM
unset STAR_WIG_PARAM
case "${DTYPE}" in
  "str_SE"|"str_PE")
   STAR_STRAND_PARAM=""; STAR_WIG_PARAM=" --outWigStrand Stranded";;
  "unstr_SE"|"unstr_PE")
   STAR_STRAND_PARAM=" --outSAMstrandField intronMotif"; STAR_WIG_PARAM=" --outWigStrand Unstranded";;
esac

#### Map reads with STAR
STAR ${STAR_COMMON_PARAM} ${STAR_RUN_PARAM} ${STAR_BAM_PARAM} ${STAR_STRAND_PARAM}

#### Make bedGraph
mkdir Signal
STAR --runMode inputAlignmentsFromBAM --inputBAMfile Aligned.sortedByCoord.out.bam --outWigType bedGraph \
     ${STAR_WIG_PARAM} --outFileNamePrefix ./Signal/ --outWigReferencesPrefix chr
mv Signal/Signal*bg .

#### Convert bedGraph to bigWig
grep -v chrEBV ${STAR_INDEX}/chrNameLength.txt > chrNL.txt # removing chrEBV
LC_COLLATE=C
case "${DTYPE}" in
  "str_SE"|"str_PE")
   STRAND=("minus" "plus")
   for ISTR in 1 2; do
     for MAPSTAT in Unique UniqueMultiple; do
       grep -v chrEBV Signal.${MAPSTAT}.str${ISTR}.out.bg | sort -k1,1 -k2,2n -T ./ > sig.tmp
       bedGraphToBigWig sig.tmp chrNL.txt Signal.${MAPSTAT}.${STRAND[ ($ISTR-1) ]}.bw
     done
   done ;;
  "unstr_SE"|"unstr_PE")
   for MAPSTAT in Unique UniqueMultiple; do
      grep -v chrEBV Signal.${MAPSTAT}.str1.out.bg | sort -k1,1 -k2,2n -T ./ > sig.tmp
      bedGraphToBigWig sig.tmp chrNL.txt Signal.${MAPSTAT}.unstranded.bw
   done ;;
esac
rm chrNL.txt sig.tmp

#### Prepare BAM for RSEM
SORT_RAM=30G
samtools view -H Aligned.toTranscriptome.out.bam > Tr.sam
case "${DTYPE}" in
  "str_SE"|"unstr_SE")
   samtools view Aligned.toTranscriptome.out.bam | sort -S ${SORT_RAM} -T ./ >> Tr.sam ;;
  "str_PE"|"unstr_PE")
   samtools view Aligned.toTranscriptome.out.bam | \
   awk '{printf "%s", $0 " "; getline; print}' | sort -S ${SORT_RAM} -T ./ | tr ' ' '\n' >> Tr.sam ;;
esac
samtools view -bS Tr.sam > Aligned.toTranscriptome.out.bam
rm Tr.sam

#### Set RSEM parameters
RSEM_RAM=30000  # RAM in MB
RSEM_COMMON_PARAM=" --bam --estimate-rspd --calc-ci --no-bam-output --seed 12345"
RSEM_RUN_PARAM="-p ${NTHREADS} --ci-memory ${RSEM_RAM}"
unset RSEM_DTYPE
case "${DTYPE}" in
  "str_SE")   RSEM_DTYPE="--forward-prob 0"              ;;
  "str_PE")   RSEM_DTYPE="--paired-end --forward-prob 0" ;;
  "unstr_SE") RSEM_DTYPE=""                              ;;
  "unstr_PE") RSEM_DTYPE="--paired-end"                  ;;
esac

#### Caluclate expression
rsem-calculate-expression ${RSEM_COMMON_PARAM} ${RSEM_RUN_PARAM} ${RSEM_DTYPE} \
                          Aligned.toTranscriptome.out.bam ${RSEM_INDEX} Quant >& Log.rsem
#### Diagnostic plot creation
rsem-plot-model Quant Quant.pdf

