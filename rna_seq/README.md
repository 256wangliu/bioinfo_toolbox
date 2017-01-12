# RNA-seq

## Annotation
STAR and RSEM parse "exon" feature in GTF files. You can include
custom annotations such as chromosomal tRNAs (since GENCODE only
includes mitochondrial tRNAs), and transposons.

## Genome index preparation
Before running `STARSEM_prep.sh`, make sure the following software
packages are installed or included in the path.
* [STAR](https://github.com/alexdobin/STAR)
* [RSEM](http://deweylab.biostat.wisc.edu/rsem/)

Make sure to edit the first few lines in the shell scirpt.
* `STARSEM_prep.sh`

## Mapping and transcript quantification
For running `STARSEM_run.sh`, the following software packages are
required.
* [STAR](https://github.com/alexdobin/STAR)
* [RSEM](http://deweylab.biostat.wisc.edu/rsem/)
* [Samtools](http://www.htslib.org/)
* [bedGraphToBigWig](http://hgdownload.cse.ucsc.edu/admin/exe/)
* [bigBedToBed](http://hgdownload.cse.ucsc.edu/admin/exe/)

Make sure to edit the first lines in the shell script to include
the paths of the above tools.
* `STARSEM_run.sh`

#### Note
For large annotation such as genome-wide transposons, remove
`--calc-ci` from RSEM arguments. The option makes the computation
too slow and possibly harmful to servers.