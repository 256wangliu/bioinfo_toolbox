# Small RNA analysis
The analysis flows for small RNA sequencing.
* [Adapter removal](https://github.com/weng-lab/junko/tree/master/bioinfo_toolbox/smallrna#adapter-removal)
* [(Possible) contamination removal](https://github.com/weng-lab/junko/tree/master/bioinfo_toolbox/smallrna#possible-contamination-removal)
* [miRNA](https://github.com/weng-lab/junko/tree/master/bioinfo_toolbox/smallrna#mirna)
* [tRNA](https://github.com/weng-lab/junko/tree/master/bioinfo_toolbox/smallrna#trna)
* [Small RNA length distribution](https://github.com/weng-lab/junko/tree/master/bioinfo_toolbox/smallrna#small-rna-length-distribution)

## Adapter removal
Small RNA sequencing libraries usually contain 3'adapters (and maybe
also a few 5'adapters). Before mapping, those sequences need to be
removed.

###### If you don't know the 3'adapter sequence
Try `appsr` or `dnap` program to predict 3'adapters, and then move
to the following "If you already know the 3'adapter sequence".

###### If you already know the 3'adapter sequence
Try `cutadapt`, other publicly available adapter clipping tools, or
`adapt-clip`, an in-house script.

## (Possible) contamination removal
~90% of total RNAs are rRNAs. Removing rRNA mapping reads can avoid contaminations
unless you want to investigate rRNA expression. Use `hsa_rrna.fa` (human rRNA sets)
as decoy sequences during mapping.

## miRNA
* `mir-iso.py` outputs miRNA primary transcripts in BED format with
  the miRBase GFF3 file. The script adds flanking regions
  (i.e. offset; +/-5nt in default) to search isomiRs.
* `mir-bed.py` outputs relative positions of mature miRNAs in the
  primary transcripts in BED format. The script uses the miRBase GFF3
  file. **The offsets added by `mir-iso.py` need to be specified in
  its run.**
* `mir-flt.py` assigns reads mapped to multiple miRNA species to miRNA
  species with minimum offsets of 5' and 3' ends. Input BED files
  should be pre-sorted by read and miRNA names, i.e. `sort -k4,4 -k5,5`.
* `len-dist.py` counts reads that fall in specific lengths. It would
   be useful when you want to see length distribution of mapped or a
   specific small RNA reads.

###### Example: miRNA pipeline
Before using this pipeline, make sure reads do not contain any
adapters and pieces of blacklisted sequences (e.g. rRNA).

```shell
    ## Extract miRNA primary transcripts with +/-5nt flanking regions from the genome
    $ python mir-iso.py -u 5 -d 5 hsa.gff3 | \
      bedtools getfasta -fi ${GENOME} -s -name -bed stdin -fo stdout | \
      python rm-n-fa.py - > hsa_precursor_ext5nt.fa

    ## Generate mature miRNA loci in the primary transcripts
    $ python mir-bed.py -u 5 -d 5 hsa.gff3 > hsa_mature_ext5nt.bed
    
    ## Map reads in FASTA to hairpins (any adapters should not be included in reads)
    ## The reads should be collapsed like below:
    #     >AAAAGCTGGGTTGAGAGGGCG_28
    #     AAAAGCTGGGTTGAGAGGGCG
    $ bowtie-build hsa_precursor_ext5nt.fa hsa_precursor
    $ bowtie -f -S -v 0 -a -p 4 --best --strata hsa_precursor ${CLEAN_FASTA} 2> /dev/null | \
      samtools view -F20 -b - | bedtools bamtobed -i stdin > ${BED_OUTPUT}
    
    ## Intersect mapped reads by allowing -/+5bp 5' starts and +5bp 3' ends, and then
    ## assign reads mapped to multiple miRNA species to miRNA species with minimum 5'/3' offsets
    ## Each field in the bed file will be:
    #    [0] miRBase ID for a miRNA primary transcript
    #    [1] start position of a mature miRNA in the primary transcript
    #    [2] end position of a mature miRNA in the primary transcript
    #    [3] mapped sequence with the read count
    #    [4] miRNA name
    #    [5] 5' offset from the annotation (annotation_start = mapped_start - offset)
    #    [6] 3' offset from the annotation (annotation_end = mapped_end - offset)
    ## NOTE: the file contains multimappers! Remove reduncancy or apportion those in the further process
    $ bedtools intersect -a ${BED_OUTPUT} -b hsa_mature_ext5nt.bed -wo | \
      awk 'BEGIN{OFS="\t"}{ d=$2-$8; e=$3-$9; if(d<0){d*=-1}
                            if(d<=5 && e<=5){print $1,$2,$3,$4,$10,$2-$8,$3-$9} }' | \
      sort -k4,4 -k5,5 | python mir-flt.py - | sort -k1,1 -k2,2n > map_result.bed
```

#### miRBase v21 details

* 1881 annotated hairpins in human
     >> 932 have 5p/3p mature miRNA pair annotations
     >> 949 do not have the annotations, i.e. a mature miRNA per hairpin

     Some pairs have multiple copies in the genomic regions.

          Example 1: Exact copies
              There are two copies of a hairpin that generates hsa-let-7a-5p and hsa-let-7a-3p.

          Example 2: Hairpins that generates the same sequence but sligntly different pair
              There are two hairpins that generates hsa-let-7f-5p (exactly the same sequence).
              However the sequences of the pairs are sligntly different:
                  hsa-let-7f-1-3p: CUAUACAAUCUAUUGCCUUCCC
                  hsa-let-7f-2-3p: CUAUACAGUCUACUGUCUUUCC

     So checking the strand selection, it would be better to check the expression of miRNAs per hairpin.

* 2813 mature miRNAs including redundant miRNA copies in genomic loci


# tRNA pipeline

## Annotation and reference sequencese
First of all, we will generate annotation and reference sequences of
mature and precursor tRNAs with the following python scripts. In this
method, we will use **40nt** upstream and downstream of tRNAs to
define **5' leader** and **3' trailer** sequences.
* **`trf-bed.py`** converts a tRNAscan-SE output into BED containing
  genomic tRNA loci.
* **`trf-ref.py`** generates a set of annotation and sequence files for
   mature and precursor tRNAs. For the detail of the difference
   between mature and precursor tRNAs, see [tRNA
   biology](https://github.com/weng-lab/junko/tree/master/bioinfo_toolbox/smallrna#trna-biology)
   section below. The details of `trf-ref.py` outputs are as follows:

|  Format  |      File Name       |            Description                   |
|:--------:|:--------------------:|------------------------------------------|
|  FASTA   | `precursor_trna.fa`  | Precursor tRNA sequences (non-redundant) |
|  FASTA   | `mature_trna.fa`     | Mature tRNA sequences (non-redundant)    |
|   BED    | `precursor_trna.bed` | Precursor tRNA annotation                |
|   BED    | `mature_trna.bed`    | Mature tRNA annotation                   |
|   BED    | `${PREFIX}_clean.bed`| tRNA annotation on genomic loci          |

In the 4th column of `*_trna.bed`, you will find annotations of tRNA (sub)sequences:

|      Flag       |      Description                   |
|:---------------:|------------------------------------|
| `5_leader`      | 5' leader sequence                 |
| `3_trailer`     | 3' trailer sequence                |
| `intron`        | Intron sequence                    |
| `mature_trna`   | mature tRNA                        |
| `exon_junction` | exon-exon junction of spliced tRNA |

One of the outputs, `${PREFIX}_clean.bed`, is very similar to the tRNA
annotation generated by `trf-bed.py`. The main difference is the 4th
column containing unique identifiers assigned as non-redundant tRNA
sequences. The identifiers are consistent with the sequence names in
FASTA.

```shell
## Run tRNAscan-SE to get position and structure outputs
## For human genome, exclude mitochondrial genome because GENCODE already has the annotation.
#     -o: position output
#     -f: structure output
#     -q: quiet mode
#     -b: no column headers for position output
$ tRNAscan-SE -o ${tRNA_POS} -f ${tRNA_SS} -q -b ${GENOME_FASTA}

## Convert tRNAscan-SE position output to BED
## Remove undefined anticodon tRNAs
$ sort -V -k1,1 -k3,3n ${tRNA_POS} | python trf-bed.py - | grep -v "???" > genomic_trna.bed

## Add mitochondrial tRNAs from GENCODE GTF annotation file if it is human dataset
$ grep Mt_tRNA ${GENCODE_GTF} | \
  awk 'BEGIN{OFS="\t"; FS="\t"}
       { if($3=="transcript"){
             split($9,tmp,"gene_name "); split(tmp[2],name,";");
             print $1,$4-1,$5,name[1],"0",$7,"no_intron"} }' | \
  sed 's/"//g' >> genomic_trna.bed

## Define precursor tRNA sequences (two parameters for leader and trailer sequence length)
#    5' leader sequence length = 40 bp
#    3' trailer sequence length = 40 bp
$ awk 'BEGIN{OFS="\t"}{ $4 = $4"|"$1":"$2"-"$3"|"$7"|leader=40|trailer=40";
                        $2 -= 40; $3 += 40; print }' genomic_trna.bed | \
  bedtools getfasta -name -s -bed stdin -fi ${GENOME_FASTA} -fo stdout |
  python trf-ref.py genomic_trna.bed -
  mv genomic_trna_clean.bed genomic_trna.bed
```

###### tRNA biology
tRNA precursors will experience several post-processing after
transcription. 5' leader and 3' trailer sequences and introns in tRNA
precursors will be removed, and *CCA* nucleotides will be added to the
3' ends. Interestingly for some histidine tRNAs, *G* nucleotide will
be added to the 5' ends.

## Mapping tRNA derived small RNAs
Before mapping, make sure the 3' adapter sequences are removed, and
the reads are collapsed to non-redundant reads in FASTA format. If you
don't know your 3' adapters, try [DNApi](https://github.com/jnktsj/DNApi).

To apply the entire pipeline discribed in this page, the reads should
be collapsed with the read count delimited with `_` in a FASTA header:

    >AAAAGCTGGGTTGAGAGGGCG_28
     AAAAGCTGGGTTGAGAGGGCG

To obtain reliable tRNA mapping reads, I recommend to discard reads
mapped to mature/precursor miRNAs and rRNAs.

```shell
## If you want to remove rRNA and miRNA mapping reads, add this step
$ bowtie-build ${rRNA_FASTA},${miRNA_HAIRPIN_FASTA} ${BLACKLIST}
$ bowtie -f -S -v 0 -k 1 -p 8 --un ${CLEAN_FASTA} ${BLACKLIST} ${FASTA} &> /dev/null

## Map reads in FASTA to mature and precursor tRNAs
$ bowtie-build mature_trna.fa,precursor_trna.fa tRNA
$ bowtie -f -S -v 0 -a -p 8 --best --strata tRNA ${CLEAN_FASTA} 2> /dev/null | \
  samtools view -F20 -b - | bedtools bamtobed -i stdin > ${BED_OUTPUT}
```

## Filtering mapped reads
It is common that you find low-complexity reads mapped to too many
genomic loci. Some tRNA precursors have poly(A) 5' leader and 3'
trailer sequences, however the poly(A) reads mapped to the regions
might not come from the regions (because there are too many candidate
loci to allocate!).

Many people use how many times a read is aligned to different loci in
a reference genome for filtering out ambiguously mapped reads, however
this criterion discards a significant portion of reads in small RNA
sequencing. For example, a read come from a highly expressed small RNA
and sequenced 1 million times will be discarded 1 million times when
each read is mapped to multiple loci more than the cutoff.

I recommend the following cutoff based on apportioned read counts by
mapped times on the genome:

    count = {read count of a collapsed read} / {mapped time on the genome}
    if (count < 1) discard the read

Empirically, a read should be discarded if the apportioned read count
drops below 1. Mapped times should be determined with **mapped times
on the reference genome**. You need to map reads to the genome in
addition to mapping to tRNAs.

For filtering ambiguous reads, the following shell script snippet
would be useful.

```shell
## If you haven't mapped reads to the genome, try this step
## It will take looong time, so maybe grab some coffee
$ bowtie-build ${GENOME_FASTA} ${GENOME_INDEX}
$ bowtie -f -S -v 0 -a -p 8 --best --strata ${GENOME_INDEX} ${CLEAN_FASTA} 2> /dev/null | \
  samtools view -F4 | cut -f1 | sort | uniq -c | awk '{print $2"\t"$1}' > read_maptime.tab

## Extract ambiguous reads
$ awk '{ split($1,a,"_"); read_count=a[2];
         if(read_count/$2 < 1){print "\t"$1"\t"} }' > blacklist.tab

## Filter out the blacklisted reads
$ grep -v -F -f blacklist.tab ${BED_INPUT} > ${CLEAN_BED}
$ rm read_maptime.tab blacklist.tab
```

## Assigning mapped reads to tRNA annotation
For annotating reads, we will use the following python script in this
step.

* **`trf-flt.py`** filters, assign, and annotates mapped reads. When a
    read is mapped to a mature tRNA, a precursor tRNA, or a tRNA
    pseudogene, the read will be assigned to an annotation with higher
    priority (mature > precursor > pseudogene).

Each criterion for read annotation is:

| tRNA fragment type        | Flag | Criterion |
|---------------------------|:----:|-----------|
| 5'-end derived            | `5`  | read starts at 0~5nt from the 5' end of a mature tRNA |
| 3'-end derived            | `3`  | read ends at -7nt (CCA + -5nt) ~ 0nt from the 3' end of a mature tRNA |
| 5'-leader drived          | `l`  | read starts with any position in a 5' leader sequence and ends at +/-5nt from the 5' end of the mature tRNA |
| 3'-trailer derived        | `t`  | read starts at +/-5nt from the 3'end of a mature tRNA and ends with any position in the 3' leader sequence |
| intron derived            | `i`  | read overlaps at least 5nt with an intron |
| derived from other region | `o`  | read mapped to a mature tRNA, that doesn't fall in any criteria above |

Each field in the output BED file is:

1. Mature or precursor tRNA name
2. 0-based leftmost mapping position on mature or precursor tRNA
3. 1-based rightmost mapping position on mature or precursor tRNA
4. Read with the read count delimited with `_`
5. Read information delimited with `:` (tRNA fragment type
   {`5`,`3`,`l`,`t`,`i`,`o`}, fragment length, and additional
   nucleotide overap on mature tRNAs {`G` for 5' His-tRNA,
   `CCA`,`CC`,`C` for 3' end} if any)
6. Offset from 5'/3' end for {`5`,`3`,`l`,`t`}, intron overlap length
   for {`i`}, `.` for {`o`}

**Note:** The output BED file contains a read mapped to multiple
 tRNAs, i.e. multimappers. We will apportion those for profiling the
 expression in the further process.

```shell
## Intersect mapped reads to obtain relative mapping position of mature tRNAs
## and precursor tRNAs (mainly to intron, 5'leader, and 3'trailer sequences)
$ bedtools intersect -a ${BED_OUTPUT} -b mature_trna.bed -wo > mature_reads.bed
$ bedtools intersect -a ${BED_OUTPUT} -b precursor_trna.bed -wo > precursor_reads.bed

## Annotate reads to tRNAs
$ python trf-flt.py mature_reads.bed precursor_reads.bed | sort -k1,1 -k2,2n > trna_fragment.bed

## If you want to obtain reads mapped to exon-exon junctions, try:
$ bedtools intersect -a trna_fragment.bed -b mature_trna.bed -wo | \
  grep "exon_junction" > junction_reads.bed
```

## Read length distribution

After mapping and cleaning reads, most people would be interested in
the length distribution of tRNA fragments.

* **`trf-len.py`** tallies read counts according to the read
    lengths. It also output the read counts by length for each tRNA
    fragment type (`5`, `3`, `l`, `t`, `i`, `o`).

```shell
$ python trf-len.py -m ${MIN_LENGTH} -x ${MAX_LENGTH} trna_fragment.bed > trna_fragment.len

## Example output:
# length        total   5       3       l       t       i       o
16      33115   2686.0  21721.0 76.0    118.0   0       8514.0
17      31201   1584.0  22298.0 151.0   116.0   0       7052.0
18      26735   3652.0  9637.0  140.0   33.0    0       13273.0
19      14329   4326.0  4207.0  43.0    39.0    5.0     5709.0
20      18051   8072.0  3775.0  40.0    128.0   5.0     6031.0
21      12671   2075.0  3820.0  10.0    505.0   8.0     6253.0
22      13525   2411.0  5599.0  17.0    311.0   9.0     5178.0
23      6773    1200.0  2913.0  4.0     165.0   2.0     2489.0
24      6766    860.0   3123.0  21.5    249.5   3.0     2509.0
25      6563    1226.0  2607.0  2.0     98.0    2.0     2628.0
26      13202   5401.0  3387.0  9.0     70.0    3.0     4332.0
27      12341   6299.0  4364.0  4.0     179.0   1.0     1494.0
28      22419   15433.0 4855.0  6.0     212.0   3.0     1910.0
29      24744   13445.0 10066.0 14.0    131.0   1.0     1087.0
30      66902   45566.0 20317.0 25.5    88.5    3.0     902.0
31      172477  145627.0        26108.0 19.0    101.0   5.0     617.0
32      1668588 1636919.0       29933.0 21.0    535.0   0       1180.0
33      2540096 2497018.0       42104.0 58.0    68.0    2.0     846.0
34      593126  577081.0        14940.0 107.0   190.0   3.0     805.0
35      132894  118723.0        13181.0 69.0    114.0   0       807.0
36      60300   48316.0 10206.0 88.0    70.0    2.0     1618.0
37      19729   5439.0  13308.0 86.0    430.0   6.0     460.0
38      15796   1918.0  11294.0 25.0    2128.0  9.0     422.0
39      15180   5740.0  8855.0  5.0     181.0   11.0    388.0
40      6453    560.0   5292.0  13.0    133.0   28.0    427.0
```

Each column in the output is:

1. Read length
2. Total reads in the length
3. 5'end derived tRNA fragments
4. 3'end derived tRNA fragments
5. 5'leader derived tRNA fragments
6. 3'trailer derived tRNA fragments
7. Intron derived tRNA fragments
8. tRNA fragments derived from other regions

## Profiling tRNA fragments
We will profile tRNA fragments using the following scripts:

* **`trf-cnt.py`** generates a read count matrix for each tRNA species,
    isoacceptor, and isotype for each tRNA fragment type by grouping
    the reads by the lenths. It generates **18 files** (6 tRNA
    fragment types * {tRNA species, isoacceptor, isotype}). The
    nomenclature of the outputs is
    `${PREFIX}.${tRNA_FRAGMENT_TYPE}.[species|isoaccp|isotype].tab`.
* **`trf-spp.py`** adds missing tRNAs (i.e. tRNAs that do not have any
    reads and were not detected) to the data matrix.
* **`trf-col.py`** extracts tRNA fragments in specific length and make
    an expression matrix.

```shell
## Generate read counts mapped to tRNAs for each tRNA fragment type
$ python trf-cnt.py -m ${MIN_LENGTH} -x ${MAX_LENGTH} trna_fragment.bed

## Add tRNAs that are not detected in the dataset
$ cat mature_trna.bed precursor_trna.bed | cut -f1 | sort -u > trna_names
$ python trf-spp.py trna_names ${OUTPUT_LENGTH_MATRIX}

## Extract read counts in a specific range of lengths
$ python trf-col.py -m 30 -x 35 ${FILE_1} ${FILE_2} ...
```

For the input of `trf-col.py`, it is easy to specify files as `ls *.5.species.tab`.