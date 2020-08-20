## raw_counts_cellXgene.mtx.gz
This file is a raw count matrix in gzipped MatrixMarket (.mtx) format.
Rows correspond to cells and columns correspond to genes.
The respective row names (= cell barcodes) are given in __rowNames_cellBarcodes.tsv.gz__.
The respective column names (= Ensembl gene IDs and corresponding gene symbols) are given in __colNames_genes.tsv.gz__.

## metadata.tsv.gz
This is a gzipped tab-separated values (.tsv) file with additional information for every cell.
Each row in this file corresponds to one cell in the above count matrix.
Here is a brief explanation of the columns:

column name | explanation
------------ | -------------
barcode | the cell barcode (16bp)
sample | the biological replicate (either "sample #1" or "sample #2")
n_UMIs | the number of observed unique molecular identifiers (UMIs)
n_genes | the number of observed genes (i.e. genes with at least one observed mRNA molecule)
percent_mito | the % of UMIs that correspond to mitochondrial genes (indicates apoptosis)
leiden | the [Leiden cluster](https://www.nature.com/articles/s41598-019-41695-z) that the cell was assigned to
UMAP1 | the first [UMAP](https://arxiv.org/abs/1802.03426) coordinate
UMAP2 | the second [UMAP](https://arxiv.org/abs/1802.03426) coordinate
celltype | the cell type, this is a more informative re-labeling of the leiden column
in_lineage | a boolean indicating whether the cell is part of the main NSC lineage or not, based on the leiden column
eYFP_counts | the number of observed eYFP UMIs
NeoR_counts |  the number of observed NeoR UMIs
Cre_counts |  the number of observed Cre recombinase UMIs
transgene_exp | a brief summary of transgene expression, e.g. "eYFP+ NeoR-"
transduction | * whether the cell was transduced or not, based on the column above
leiden_raw | the Leiden cluster, but without regressing out the effects of the cell cycle
UMAP1_raw | the first UMAP coordinate, but without regressing out the effects of the cell cycle
UMAP2_raw | the second UMAP coordinate, but without regressing out the effects of the cell cycle
batch | same column as "sample" but with "WT1" instead of "sample #1"

\* the value in the "transduction" column is either "transduced" (for "eYFP+ NeoR-" or "eYFP+ NeoR+"), "not transduced" (for "eYFP- NeoR+") or "unknown" (for "eYFP- NeoR-").
