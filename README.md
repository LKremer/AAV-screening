# Part I: Screening of AAV serotypes with RNA sequencing

## Data
#### Raw sequencing results (FASTQ)
These files are deposited at GEO (accession [GSE145172](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE145172))

## Results
#### All plots of barcode / serotype proportions are here:
[02_analysis/plots/](02_analysis/plots/)

#### All supplementary tables of barcode / serotype proportions are here:
[02_analysis/tables/supplement/](02_analysis/tables/supplement/)

## Code
#### The code to count the barcodes in our FASTQ files is here:
[01_counting-barcodes/AAV-lib1/Snakefile](01_counting-barcodes/AAV-lib1/Snakefile)  
[01_counting-barcodes/AAV-lib3/Snakefile](01_counting-barcodes/AAV-lib3/Snakefile)  
These scripts can be run with [Snakemake](https://snakemake.readthedocs.io).

#### The R code for calculating normalized barcode proportions and for plotting is here:
[02_analysis/Lib1-analysis.R](02_analysis/Lib1-analysis.R)  
[02_analysis/Lib3-analysis.R](02_analysis/Lib3-analysis.R)  


# Part II: Single cell RNA-seq of cells labeled with AAV1_P5
## Data
#### Our filtered single cell RNA-seq count matrix, including row- and column labels and metadata
[03_single-cell-RNA-seq/count_matrix/](03_single-cell-RNA-seq/count_matrix/)
## Code + Results
The Python and R code including intermediate results and plots can be found in these HTML files (download and then open in a web browser).
#### The Python code used to get from a raw unfiltered count matrix to our UMAP and clustering
[03_single-cell-RNA-seq/01_scRNA-seq-preprocessing.html](03_single-cell-RNA-seq/01_scRNA-seq-preprocessing.html)
(or view the [raw Python code](03_single-cell-RNA-seq/01_scRNA-seq-preprocessing.py))
#### The R code used for the analyses depicted in Figures 3 and S4h-l, including these plots
[03_single-cell-RNA-seq/02_scRNA-seq-analysis.html](03_single-cell-RNA-seq/02_scRNA-seq-analysis.html)
(or view the [raw R Markdown code](03_single-cell-RNA-seq/02_scRNA-seq-analysis.Rmd))



## Computational Methods of Part I (AAV screening):

### Counting barcode occurrences in FASTQ files
NGS-samples were sequenced and demultiplexed by the DKFZ Genomics and Proteomics Core Facility using `bcl2fastq` 2.19.0.316.
This resulted in two (paired-end) FASTQ-files per sample.
Each FASTQ consists of reads resulting from the targeted barcode amplification and up to 50% PhiX DNA that was spiked in to increase library complexity.  
Each AAV variant is associated with a unique 15-mer barcode sequence.
To quantify the most successful AAV, we simply counted how often each barcode occurred in each FASTQ file, bearing in mind the following pitfalls:
1. barcode-sequences might occur outside of the amplicon by chance, e.g. in the PhiX genome
2. barcodes might have sequencing errors
3. barcodes occur on the forward and reverse strand  

To circumvent issues 1 and 2, we opted for a strategy where we only count barcodes matching the expected amplicon structure. 
This was achieved with the following regex (regular expression; defines a text search pattern):
`(?<=[NGCAT]{33}TGCTC)[NGCAT]{15}(?=CAGGG[NGCAT]{45})`.
Variable 15-mers `[NGCAT]{15}` are only counted if they are flanked by the expected regions `TGCTC` and `CAGGG`.
Furthermore, we enforce a minimum of 33 upstream nucleotides and 15 downstream nucleotides, in addition to the flanking regions, to only count 15-mers at the expected position.
15-mers matching this regex were extracted and counted with the standard GNU command-line tools `grep`, `sort` and `uniq`.
15-mers sequenced from the reverse strand were counted with an equivalent reverse complement regex and added to the forward counts.
### Assigning barcodes to AAV capsids
Raw 15-mer counts were further processed in R.
Most observed 15-mers matched a known barcode exactly (library #1: 74%, library #3: 87%), which allowed us to assign them to a unique AAV variant.
The remaining 15-mer counts were added to the counts of the closest known barcode, allowing for a maximum of two mismatches.
### Normalization
Each sequenced sample corresponds to one tube with up to 500 FAC-sorted cells.
To downweigh samples with lower cell numbers, barcode counts were scaled by the respective number of FACS events (usually 500, Supplementary Table 2).
Barcode counts of the same cell type and biological replicate (set) were then summed.
The AAV libraries used for transduction contain slightly unequal proportions of AAV variants, which means that some AAV variants may have an advantage due to increased starting concentration.
To remedy this problem, barcode counts were further scaled by their abundance in the transduction library, so that barcode counts corresponding to more frequent AAV capsids were decreased and vice versa.
### Identification of candidate AAVs with high transduction efficiency 
For each sample, normalized barcode proportions were calculated by dividing the normalized barcode counts by the total number of valid barcodes.
The mean of these proportions was then used to rank AAVs within and across celltypes (Figure 1d-i).
Two AAVs that performed consistently well across replicates and in both experiments were chosen for further validation.
All scripts used in the analysis are available in this GitHub repo here.
