# cerebellum-atlas-analysis
Code and scripts used in analysis of snRNAseq mouse cerebellum atlas, as described in our recent
[preprint](https://www.biorxiv.org/content/10.1101/2020.03.04.976407v1.abstract). 

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
data are available in `data`.

