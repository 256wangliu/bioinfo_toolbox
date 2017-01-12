# ChIP-seq
## Read mapping and post-alignment filtering
For read mapping, the following tools will be used in the pipeline:
* [BWA](http://bio-bwa.sourceforge.net/)
* [Samtools](http://www.htslib.org/)
* [bedtools](http://bedtools.readthedocs.org/en/latest/)
* [Picard tools](http://broadinstitute.github.io/picard/)

Before using the mapping pipeline, make sure edit the software paths in
the shell script.

* `map_se_chip.sh` is for **single-end**.
* `map_pe_chip.sh` is for **paired-end**.

>Shell scirpts are customized version of [ENCODE 3 ChIP-seq pipeline]
(https://docs.google.com/document/d/1lG_Rd7fnYgRpSIqrIfuVlAz2dW1VaSQThzk836Db99c/edit#).

## Strand cross-correlation and fragment length estimation (optional)
For computing strand cross-correlation and estimating fragment length,
the following tools will be used in the pipeline:
* [phantompeakqualtools](https://code.google.com/p/phantompeakqualtools/)
* [Samtools](http://www.htslib.org/)
* [bedtools](http://bedtools.readthedocs.org/en/latest/)

Make sure to edit the sotware paths in the shell script.
* `cross_corr_chip.sh`

If you don't know the length of ChIPped protein sites, use this
fragment length for peak calling.

## Peak calling
For calling peaks, use a relaxed P-value setting for subsequent IDR
analysis (examples with MACS2 as follows). It is not tested, but if
the library is paired-end, specify `BAMPE` to use real insert length.

```shell
    ## Peak calling
    $ macs2 callpeak -t ${CHIP} -c ${CTRL} -f ${FORMAT} -n ${PREFIX} -g ${GENOME} \
                     -p 1e-2 --keep-dup all -B --SPMR

    ## Peak calling with an estimated fragment length by 'phantompeakqualtools'
    ## For single-end data, it may be more accurate?
    $ macs2 callpeak -t ${CHIP} -c ${CTRL} -f ${FORMAT} -n ${PREFIX} -g ${GENOME} \
                     --shift 0 --nomodel --extsize ${FRAGLEN} \
                     -p 1e-2 --keep-dup all -B --SPMR    
```

For broad peaks, simply add `--broad` to the above commands.

## Irreproducible Discovery Rate (IDR)
After calling peaks with a relaxed P-value setting, calculate IDRs to
extract reproducible peaks among samples. For computing IDRs of each
pair and filtering peaks, the following tools will be needed for the pipeline:
* [bedtools](http://bedtools.readthedocs.org/en/latest/)
* [IDR code](https://sites.google.com/site/anshulkundaje/projects/idr/idrCode.tar.gz?attredirects=0)

Make sure to edit the software paths in the shell script.
* `peak_idr_chip.sh`

## Master peak set
When the curated peaks in each sample are ready, simply merge them to make a master peak set:

```shell
    $ cat *.narrowPeak | sort -k1,1 -k2,2n | \
      bedtools merge -i stdin -c 4 -o distinct > master.narrowPeak
```

## Peak annotation
Although there are bunch of tools that can annotate detected peaks,
but most of them have some limitations. For example, a peak is
associated only one gene by a peak annotation software package
 although the peak locates in two genes transcribed bidirectionally.

So I decided to use just `bedtools closest` to annotate peaks.  To
annotate peaks, annotation files, i.e. TSS loci should be
prepared. Since GENCODE doesn't contain chromosomal tRNAs for human
annotation, prediction results of tRNAscan-SE can be merged to GENCODE
annotation. If you don't care tRNAs, I think it's okay to use just
genes or transcripts.

```shell
   # Prepare annotation (GENCODE + tRNA)
   # All BED files should be sorted first
   sort -k1,1 -k2,2n TSS_1base.bed > TSS_1base.sorted.bed
   sort -k1,1 -k2,2n master.narrowPeak | awk 'BEGIN{OFS="\t"}{print $1,$2,$3}' > master.sorted.bed

   # Find the closest gene to a peak
   # If there is no gene close to a peak, the position of the record will be "-1"
   bedtools closest -D b -t all -a master.sorted.bed -b TSS_1base.sorted.bed > annotated_peaks.bed
```

## Generate coverage tracks
For visualizing ChIP-seq tracks, there are several ways in terms of
what values you would like to display. If you want to see log
likelihood, P-values, Q-values, etc. between ChIPped and control
samples, `macs2 bdgcmp` gives you what you want.

In this section, I would like to focus on just read counts, i.e.
genome coverage tracks.

```shell
# simply counting read coverage
bedtools genomecov -bg -i ${BED} -g ${GENOME_LEN} | sort -k1,1 -k2,2n > ${OFPREFIX}.bg
bedGraphToBigWig ${OFPREFIX}.bg ${GENOME} ${OFPREFIX}.bw

# normalized read coverage
MAPPED_READS=`wc -l ${BED} | awk '{print $1}'`
FACTOR=`echo ${MAPPED_READS} | awk '{printf "%.4f\n", 1000000/$1}'`
bedtools genomecov -bg -i ${BED} -g ${GENOME} -scale ${FACTOR} | sort -k1,1 -k2,2n > ${OFPREFIX}.rpm.bg
bedGraphToBigWig ${OFPREFIX}.rpm.bg ${GENOME} ${OFPREFIX}.rpm.bw
```

## Differential peak analysis
Before counting reads mapped to peak regions or TSSs, you may want to
extend the length of reads (say, from 36bp to 147bp). Make sure to
convert the BAM files generated the mapping pipeline to BED file, and
simply run `ext-bed-reads.py`.

It's fine to sum extention length simply to the coordinate.
`ext-bed-reads.py` takes into account cigar hard/soft masks etc.

#### Obtaining fragment length from BEDPE

```shell
zcat ${BEDPE} | \
awk 'BEGIN{ OFS="\t"; FS="\t" }
          { chrom=$1; beg=$2; end=$6; readname=$7;
            if( $2 > $5 ){ beg = $5 }
            if( $3 > $6 ){ end = $3 }
            print chrom,beg,end,readname,".","." }' - > "${OUTPUT}.fragment.bed"
```

#### Count reads
```shell
bedtools intersect -a ${FRAGBED} -b ${PEAKS} -wo | \
awk 'BEGIN{ OFS="\t" }{ print $7,$8,$9,$4 }' | \
sort -k1,1 -k2,2n - | bedtools merge -i stdin -c 4 -o count_distinct > ${COUNT}.count
```
