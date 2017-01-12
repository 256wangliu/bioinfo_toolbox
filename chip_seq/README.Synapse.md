# Downloading datasets from Synapse
For HD study, I decided to use some H3K4me3 libraries on psychENCODE
Synapse portal. This is a note for me who is forgetful to remember how
I downloaded data.

## Querying with `synapseclient`
### Installation
First of all, we need to know how to query what datasets are
available. To do so, we need to download and install
[`synapseclient`](https://pypi.python.org/pypi/synapseclient). It
seems the R package is also available, but I chose python this time
because shell command `synapse` will be also installed with the python
package.

Just install it locally if you run it on a server:
```shell
$ python setup.py build
$ python setup.py install --user
```

### Querying the data
Well, I want to avoid to read whole documentation of synapseclient...
What we need to get for downloading datasets of interest is just a set
of Synapse IDs. I found this website is very useful to learn how to
query quickly: [Link to the part of psychENCODE website on
Synapse](https://www.synapse.org/#!Synapse:syn4921369/wiki/392338).

See **Annotation Dictionary** in the page to know the keywords for
querying.


Since I wanted to download all NeuN+/- H3K4me3 ChIP-seq datasets
sequenced from (dorsolateral) prefrontal cortex, I ran:

```shell
$ synapse query 'SELECT id,disease,tissueType,assayTarget,cellType FROM file \
                    WHERE projectId=="syn4921369" and dataType=="histone modification" and \
                          organism=="Homo sapiens" and fileType=="bam" and study=="EpiMap"' | \
     grep Prefrontal | grep NeuN | grep H3K4me3 | cut -f1 > ID_LIST
```

It may look pretty primitive, but now we have all the Synapse IDs in
the `ID_LIST` file.

### Downloading BAMs
Now, time to use `synapseclient`. To skip typing username and password,
prepare `~/.synapseConfig` under the home directory. The format is as
follows:

```
[authentication]
username: USERNAME
password: PASSWORD
```

I made a simple python script that uses `synapseclient`. If you want
to check how the module works, check out `pec_bluk.py` (although it's
super simple):

```shell
$ python pec_bulk.py ${downlaodPath} ID_LIST
```

If you are impatient and want to download files in parallel, the
following command may work:

```shell
$ for i in `cat ID_LIST | cut -f1`; do nohup echo $i > $i.txt; python pec_bulk.py ${downloadPath} $i.txt > /dev/null 2> log & done
```

And, delete `*.txt` after downloading.

For more documentations, this is the documantation page:
[Link](http://python-docs.synapse.org/index.html).

### Restoring FASTQs from BAMs
Since I wanted to map reads to hg38, not hg19 used in psychENCODE, I
restored FASTQ reads using `Picard`.

```shell
$ java -Xmx10g -jar ~/soft/picard-tools-1.130/picard.jar SamToFastq I=${BAMFILE} FASTQ=${PREFIX}_R1.fastq SECOND_END_FASTQ=${PREFIX}_R2.fastq
```
