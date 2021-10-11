# cerebellum-atlas-analysis
Code and scripts used in the analysis of a comprehensive snRNAseq mouse cerebellar cortex atlas, as described in our recent
[publication](https://www.nature.com/articles/s41586-021-03220-z). 

Example scripts used for clustering and cell type annotation are available in `annotation`. 

## Figures

Markdown files to reproduce the main figures are available in `main_figures` and below:

[Figure 1: Cerebellar cell type annotation](https://raw.githack.com/MacoskoLab/cerebellum-atlas-analysis/master/main_figures/cb_figure1.html)

[Figure 2: Characterization of spatial variation](https://raw.githack.com/MacoskoLab/cerebellum-atlas-analysis/master/main_figures/cb_figure2.html)

[Figure 3: Cross-cluster continuity of neuronal populations](https://raw.githack.com/MacoskoLab/cerebellum-atlas-analysis/master/main_figures/cb_figure3.html)

[Figure 4: Characterization of MLI2  cells](https://raw.githack.com/MacoskoLab/cerebellum-atlas-analysis/master/main_figures/cb_figure4.html)

[Figure 4: Development of MLI2  cells](https://raw.githack.com/MacoskoLab/cerebellum-atlas-analysis/master/main_figures/cb_figure4_DE.html)

Some custom functions for Seurat/LIGER object manipulation and analysis (not yet available in LIGER) are included in `cb_custom_functions.R` and `src`. 

## Data

Processed count matrices can be downloaded from [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE165371) and [Single Cell Portal](https://singlecell.broadinstitute.org/single_cell/study/SCP795). Additional metadata and smFISH
data are available in `data`. Original individual cell type analysis objects (Seurat/LIGER) can be made available upon request; UMAP and t-SNE coordinates corresponding to many individual analysis objects are available in `data`. 

First round annotations for all recovered nuclei profiles (780,553 + 11,812 FACS sorted putative Purkinje profiles) are available in `data/full_cb_metadata.csv`. 

