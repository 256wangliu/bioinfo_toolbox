# bioinfo_toolbox

The directory contains scripts or some other information that might be
useful for bioinformatic analysis.  The detailed description of each
program can be obtained by typing `<script-name> -h` or `<script-name>
--help`.

## Contents
* [ChIP-seq](https://github.com/weng-lab/junko/tree/master/bioinfo_toolbox/chip_seq)
* [RNA-seq](https://github.com/weng-lab/junko/tree/master/bioinfo_toolbox/rna_seq)
* [Small RNA analysis](https://github.com/weng-lab/junko/blob/master/bioinfo_toolbox/smallrna)
* [BigWig for UCSC Genome Browser](https://github.com/weng-lab/junko/blob/master/bioinfo_toolbox/README.md#bigwig)
* [Sequences](https://github.com/weng-lab/junko/blob/master/bioinfo_toolbox/README.md#sequences)
* [FASTA](https://github.com/weng-lab/junko/blob/master/bioinfo_toolbox/README.md#fasta)
* [Protein](https://github.com/weng-lab/junko/blob/master/bioinfo_toolbox/README.md#protein)

## BigWig
To visualize read pileup in UCSC Genome Browser, we need BigWig files.

* Raw read counts
```shell
bedtools genomecov -bg -i ${BED} -g ${GENOME} | sort -k1,1 -k2,2n > ${OFPREFIX}.bg
bedGraphToBigWig ${OFPREFIX}.bg ${GENOME} ${OUTPUT}.bw
```

* Normalized read counts by mapped reads (i.e. RPM)
```shell
# obtain mapped reads from bed-converted reads
MAPPED_READS=`wc -l ${BED} | awk '{print $1}'`

# calculate factor to multiply
FACTOR=`echo ${MAPPED_READS} | awk '{printf "%.4f\n", 1000000/$1}'`

bedtools genomecov -bg -i ${BED} -g ${GENOME} -scale ${FACTOR} | sort -k1,1 -k2,2n > ${OFPREFIX}.rpm.bg
bedGraphToBigWig ${OFPREFIX}.rpm.bg ${GENOME} ${OUTPUT}.rpm.bw
```

## Sequences
* `hsa_rrna.fa` contains human nuclear rRNA sequences.
* [hg38/hg19 sponge sequences]
  (http://nar.oxfordjournals.org/content/early/2015/07/09/nar.gkv671.full)

###### Example: Remove contaminated reads

```shell
    # remove reads by mapping to the sequences in blacklist
    $ bowtie-build hsa_rrna.fa rRNA
    $ bowtie -f -S -v 0 -k 1 -p 4 --un ${CLEAN_FASTA} rRNA ${FASTA} &> /dev/null
```

## FASTA
* `rm-n-fa.py` removes sequences composed of only 'N' bases.

## Protein
* `aa-kda.py` calculates protein moledular weights in kDa from FASTA
  file(s).

## LegoPlot
* `lego-plot.py` displays somatic point mutations per three-nucleotide context
