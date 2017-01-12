#!/bin/bash

# Make index for STAR and RSEM. This is customized verison
# of ENCODE long RNA-seq pipeline. The original shell script
# can be found from here: https://github.com/ENCODE-DCC/long-rna-seq-pipeline/

#### Software path
export PATH=$PATH:/home/tsujij/soft/rsem-1.2.22/
export PATH=$PATH:/home/tsujij/soft/STAR_2.4.2a/source/

#### Usage
function usage {
cat <<EOF
Usage: $0 -a <star-index> -b <rsem-index> -g <genome> -t <annotation.gtf>

Prepare genome index for STAR and RSEM.

Required arguments:
  -a <star-index>       STAR index directory name with *full* path.
  -b <rsem-index>       RSEM index directory name with *full* path.
                        Both directories shouldn't be created beforehand,
                        i.e., the shell script will create both.
  -g <genome>           Reference genome in FASTA.
  -t <annotation.gtf>   Transcript annotation in GTF.

Options:
  -h                    Show this help message and exit.
  -p <cpu-num>          Number of CPUs for STAR (default: 8 cores).
  -o <prefix>           Prefix of index files inside of the STAR '-s' and
                        RSEM '-r' directories (default: prefix of genome FASTA).
EOF
}

#### Parameters
NTHREADS=8
unset STAR_INDEX
unset RSEM_INDEX
unset GENOME
unset OFPREFIX
unset GTF

#### Read arguments and options
[[ $# -eq 0 ]] && usage && exit 1;
while getopts "ha:b:g:t:p:o:" OPT
do
  case ${OPT} in
    "a") STAR_INDEX=${OPTARG} ;;
    "b") RSEM_INDEX=${OPTARG} ;;
    "g") GENOME=${OPTARG}     ;;
    "t") GTF=${OPTARG}        ;;
    "p") NTHREADS=${OPTARG}   ;;
    "o") OFPREFIX=${OPTARG}   ;;
    "h") usage && exit 1      ;;
     * ) usage && exit 1      ;;
  esac
done
[[ $# -lt 8 ]] && echo "$0: need more arguments" && exit 1;

#### Check arguments
if [ -d "${STAR_INDEX}" ]
  then echo "$0: ${STAR_INDEX} directory exists" && exit 1; fi
if [ -d "${RSEM_INDEX}" ]
  then echo "$0: ${RSEM_INDEX} directory exists" && exit 1; fi
if [ ! -f "${GENOME}" ]
  then echo "$0: can't read genome" && exit 1; fi
if [ ! -f "${GTF}" ]
  then echo "$0: can't read annotation file" && exit 1; fi
if [ $( echo ${NTHREADS} | sed 's/^[+0-9][0-9]*//' | wc -c ) -ne 1 ]
  then echo "$0: can't interpret -p ${NTHREADS}" && exit 1; fi
if [ -z "${OFPREFIX}" ]
  then OFPREFIX=`basename ${GENOME} | awk -F ".fasta" '{print $1}' | awk -F ".fa" '{print $1}'`; fi

GTF=`readlink -f ${GTF}`
GENOME=`readlink -f ${GENOME}`

#### RSEM
mkdir ${RSEM_INDEX} && cd ${RSEM_INDEX}
rsem-prepare-reference --gtf ${GTF} ${GENOME} ${RSEM_INDEX}/${OFPREFIX}

#### STAR
mkdir ${STAR_INDEX} && cd ${STAR_INDEX}
STAR --runThreadN ${NTHREADS} --runMode genomeGenerate --genomeDir ${STAR_INDEX} \
     --outFileNamePrefix ${OFPREFIX} --genomeFastaFiles ${GENOME} --sjdbGTFfile ${GTF} \
     --sjdbOverhang 100
