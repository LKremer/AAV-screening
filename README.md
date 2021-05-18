## *In vivo* high-throughput screening of novel adeno-associated viral capsids identifies variants for transduction of adult neural stem cells within the subventricular zone

This repository contains all code and data that was generated as part of the study above by Dehler and Kremer *et al*.

# Part I: Screening of AAV serotypes with RNA sequencing

## Data
#### Raw sequencing results:
The FASTQ files are deposited at GEO (accession [GSE145172](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE145172), AAV_library_1 and AAV_library_3).

#### Information about input libraries #1 and #3
The table [02_analysis/tables/input_library_info.csv](02_analysis/tables/input_library_info.csv) lists the relative abundance of each AAV variant in the input library and the barcode that was assigned to each variant.

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
#### Raw sequencing results:
The FASTQ files are deposited at GEO (accession [GSE145172](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE145172), single-cell RNA-seq_sample #1 and #2).

#### Our filtered single cell RNA-seq count matrix, including row- and column labels and metadata:
[03_single-cell-RNA-seq/count_matrix/](03_single-cell-RNA-seq/count_matrix/)

## Code + Results

#### The shell code used for read mapping and transcript quantification:
[03_single-cell-RNA-seq/README.md](03_single-cell-RNA-seq/README.md)  

#### The Python code used to get from a raw unfiltered count matrix to our UMAP and clustering:
[03_single-cell-RNA-seq/01_scRNA-seq-preprocessing.html](03_single-cell-RNA-seq/01_scRNA-seq-preprocessing.html)
(or view the [raw Python code](03_single-cell-RNA-seq/01_scRNA-seq-preprocessing.py))

#### The R code used for the analyses depicted in Figures 3 and S4h-l, including these plots:
[03_single-cell-RNA-seq/02_scRNA-seq-analysis.html](03_single-cell-RNA-seq/02_scRNA-seq-analysis.html)
(or view the [raw R Markdown code](03_single-cell-RNA-seq/02_scRNA-seq-analysis.Rmd))  
If you want to view the HTML files above, please download them and open them in a web browser.
