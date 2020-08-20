## Raw sequencing data

All raw reads in FASTQ format are currently pending acceptation at the NCBI Gene Expression Omnibus (GEO).
Once they are accepted, the raw FASTQ files will be available under accession [GSE145172](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE145172).

## Read mapping and transcript quantification

Reads were pseudoaligned and quantified with [kallisto|bustools](https://github.com/pachterlab/kb_python).
We pseudoaligned the reads to a modified mm10 mouse genome, where we manually added entries for our transgenes of interest to the genome FASTA and GTF files.
First, we build the kallisto index as follows:

```bash
kb ref -i index.idx -g t2g.txt \
  -f1 cdna.fa -f2 intron.fa \
  -c1 cdna_t2c.txt -c2 intron_t2c.txt \
  --lamanno mm10_with_transgenes.fasta mm10_with_transgenes.gtf
```

Note that we built a special kind of index that can furthermore distinguish between spliced and unspliced mRNAs (--lamanno).
This option is not really required to reproduce the results reported in the paper, though.

We then pseudoalign and quantify our reads with kallisto|bustools as follows:

```bash
# sample #1
kb count --h5ad --verbose \
  -t 24 -m 110G \
  -i index.idx -g t2g.txt \
  -x 10xv3 -o sample_1 \
  -c1 cdna_t2c.txt -c2 intron_t2c.txt \
  --lamanno --filter bustools \
  AS-438105-LR-48144_R1.fastq.gz AS-438105-LR-48144_R2.fastq.gz \
  AS-438107-LR-48144_R1.fastq.gz AS-438107-LR-48144_R2.fastq.gz \
  AS-438109-LR-48144_R1.fastq.gz AS-438109-LR-48144_R2.fastq.gz \
  AS-438111-LR-48144_R1.fastq.gz AS-438111-LR-48144_R2.fastq.gz
  
# sample #2
kb count --h5ad --verbose \
  -t 24 -m 110G \
  -i index.idx -g t2g.txt \
  -x 10xv3 -o sample_2 \
  -c1 cdna_t2c.txt -c2 intron_t2c.txt \
  --lamanno --filter bustools \
  AS-438113-LR-48144_R1.fastq.gz AS-438113-LR-48144_R2.fastq.gz \
  AS-438115-LR-48144_R1.fastq.gz AS-438115-LR-48144_R2.fastq.gz \
  AS-438117-LR-48144_R1.fastq.gz AS-438117-LR-48144_R2.fastq.gz \
  AS-438119-LR-48144_R1.fastq.gz AS-438119-LR-48144_R2.fastq.gz
  ```
  
This produces two unfiltered count matrices in h5ad-format here:
"sample_1/counts_unfiltered/adata.h5ad" and "sample_2/counts_unfiltered/adata.h5ad".
These files can be opened in Scanpy for further filtering and pre-processing of our single cell data, "01_scRNA-seq-preprocessing.html".
