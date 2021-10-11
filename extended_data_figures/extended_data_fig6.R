# Reproducing Extended Data Fig 6 panels

## Load utility functions and data

library(Seurat)
# v3.1.5
library(monocle3)
library(ggplot2)
library(cowplot)

# Additional plotting functions
source('../src/fig4_monocle.R')

# Define an extra utility function
Matrix.column_norm <- function(A){
  if (class(A)[1] == "dgTMatrix") {
    temp = summary(A)
    col.names = colnames(A)
    row.names = rownames(A)
    A = sparseMatrix(i=temp[,1],j=temp[,2],x=temp[,3])
    rownames(A) = row.names
    colnames(A) = col.names
  }
  A@x <- A@x / rep.int(Matrix::colSums(A), diff(A@p))
  return(A)
}


cb_dev_s = readRDS('../data/cb_dev_annotated.RDS')

# Create monocle object
cb_dev_mon = new_cell_data_set(expression_data = cb_dev_s@assays$RNA@counts)

# add batch information 
cb_dev_mon@colData[['sample']] = sapply(colnames(cb_dev_mon), function(x) {
  splits = strsplit(x, '_')[[1]]
  # if from the dev datasets
  if (splits[2] %in% c(1, 2)) {
    paste0(splits[1], '_', splits[2])
  } else {
    paste0('P60_', splits[1], '_', splits[2])
  }
})
# set gene row name 
rowData(cb_dev_mon)[['gene_short_name']] = row.names(rowData(cb_dev_mon))

cb_dev_mon@colData[['age']] = sapply(cb_dev_mon@colData$sample, function(x) {
  strsplit(x, '_')[[1]][1]
})

# add liger-based clusters (cell types)
cb_dev_mon@colData[['liger_clusters']] = cb_dev_s@meta.data$cluster

# basic preprocessing
cb_dev_mon <- preprocess_cds(cb_dev_mon, num_dim = 30, use_genes = cb_dev_s@assays$RNA@var.features)
# set UMAP
reducedDims(cb_dev_mon)$UMAP = cb_dev_s@reductions$umap@cell.embeddings

# generate partitions/cluster from UMAP
cb_dev_mon <- cluster_cells(cb_dev_mon, resolution = 5e-4, random_seed = 1)

cb_dev_mon <- learn_graph(cb_dev_mon)


# Compute pseudotime
# Select top node as root along with rightmost node of 
# non-MLI partition (if running in interactive)

# these nodes are 
root_pr_nodes = c('Y_62', 'Y_243')

cb_dev_mon = order_cells(cb_dev_mon, root_pr_nodes = root_pr_nodes)

## Panel A

cb_dev_mon@colData$age = factor(cb_dev_mon@colData$age, 
                                levels = c('E18', 'P0', 'P4', 'P8', 'P12', 'P16'), ordered = T)

# compute normalization as in liger
norm_genes = Matrix.column_norm(cb_dev_s@assays$RNA@counts)

pdf('extended_data_fig6_panelA.pdf', width = 7.5, height = 5, useDingbats = F)
plot_list = list()
for (gene in c( 'Tfap2b', 'Ascl1',  'Neurog1', 'Neurog2', 'Pax2', 'Klhl1')) {
  plot_list[[gene]] = plot_cells_colors(cb_dev_mon,
                                     genes = c(gene),
                                     label_cell_groups=FALSE,
                                     label_leaves=FALSE,
                                     label_branch_points=FALSE,
                                     label_roots = F,
                                     graph_label_size=3,
                                     trajectory_graph_segment_size = 0.5,
                                     trajectory_graph_color = 'grey60',
                                     cell_size = 0.5,
                                     normalized_values = norm_genes) + theme(legend.position = 'none') +
    # remove axes 
    NoAxes()
}
plot_grid(plotlist = plot_list, nrow = 2)
dev.off()

## Panel C 

pdf('extended_data_fig6_panelC.pdf', width = 4, height = 3, useDingbats = F)
plot_cells(cb_dev_mon,
           color_cells_by = 'liger_clusters',
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           label_roots = F,
           graph_label_size=3,
           trajectory_graph_segment_size = 0.5,
           trajectory_graph_color = 'grey60',
           cell_size = 0.5)
dev.off()


## Panel D

# can we get some genes that are in the shared part of the trajectory and also in the MLI2 cluster
traj_genes_dev <- graph_test(cb_dev_mon, neighbor_graph="principal_graph", cores=6)
pr_dev_ids <- row.names(subset(traj_genes_dev, q_value < 1e-3))

# find gene modules 
gene_modules_df <- find_gene_modules(cb_dev_mon[pr_dev_ids,], resolution = 1e-1, random_seed = 1)

# 34 modules
table(gene_modules_df$module)

# select modules associated with three MLI related clusters (not other interneuron populations)
gene_modules_df_sub = gene_modules_df[gene_modules_df$module %in% c(5, 14, 24, 31, 23, 7, 10, 1, 6, 2, 26, 22,
                                                                 16, 8, 20),]

gene_modules_df_sub$module = factor(droplevels(gene_modules_df_sub$module), levels = c(5, 14, 24, 31, 23, 7, 10, 1, 6, 2, 26, 22,
                                                                                     16, 8, 20),
                                   ordered = T)

# note that we are calculating the matrix using the full dataset (including the other int branch)
agg_mat_master <- aggregate_gene_expression(cb_dev_mon, gene_group_df = gene_modules_df_sub)
row.names(agg_mat_master) = stringr::str_c("Module ", row.names(agg_mat_master))

# use Seurat functionality for this split line heatmap

# first remove corrupted dimreduc slots 
cb_dev_s@reductions$tsne@assay.used = "RNA"
cb_dev_s@reductions$inmf@assay.used = "RNA"
cb_dev_s@reductions$umap@assay.used = "RNA"

cb_dev_s_sub = subset(cb_dev_s, idents= c('Early_MLI', 'MLI1', 'MLI2'))

cellname_ordered = cb_dev_s_sub@meta.data %>% 
  tibble::rownames_to_column('names') %>%
  group_by(cluster) %>% arrange(Pseudotime, .by_group = TRUE)

# artificially add these to "scale.data" slot 
cb_dev_s_sub@assays$RNA@scale.data = agg_mat_master[,row.names(cb_dev_s_sub@meta.data)]

pdf('extended_data_fig6_panelD.pdf', width = 8, height = 5, useDingbats = F)
DoHeatmap(cb_dev_s_sub, assay = 'RNA', features = row.names(agg_mat_master), cells = cellname_ordered$names,
          group.by = 'ident', group.bar = T, draw.lines = T, slot = 'scale.data', label = F,
          lines.width = 50)
dev.off()

## Panel E

# plot relevant module genes -- pick only 6
# original width was 9.5, height 6
pdf('extended_data_fig6_panelE.pdf', width = 7.5, height = 5, useDingbats = F)
plot_list = list()
for (gene in c( 'Serpine2',  'St18', 'Palmd', 'Ar', 'Fam135b', 'Gjd2')) {
  plot_list[[gene]] = plot_cells_colors(cb_dev_mon,
                                     genes = c(gene),
                                     label_cell_groups=FALSE,
                                     label_leaves=FALSE,
                                     label_branch_points=FALSE,
                                     label_roots = F,
                                     graph_label_size=3,
                                     trajectory_graph_segment_size = 0.5,
                                     trajectory_graph_color = 'grey60',
                                     cell_size = 0.5,
                                     normalized_values = norm_genes) + theme(legend.position = 'none') +
    # remove axes 
    NoAxes()
}
plot_grid(plotlist = plot_list, nrow = 2)
# modules here are 5, 7, 6, 2, 20, 26
dev.off()

# Package versions for panels above  
sessionInfo()
# R version 3.6.3 (2020-02-29)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Debian GNU/Linux 10 (buster)
# 
# Matrix products: default
# BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.8.0
# LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.8.0
# 
# locale:
#   [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8        LC_COLLATE=C.UTF-8    
# [5] LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8    LC_PAPER=C.UTF-8       LC_NAME=C             
# [9] LC_ADDRESS=C           LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
# 
# attached base packages:
#   [1] stats4    parallel  stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] pheatmap_1.0.12             dplyr_0.8.5                 ggplot2_3.3.0               liger_0.5.0                
# [5] patchwork_1.0.0             Matrix_1.2-18               cowplot_1.0.0               Seurat_3.1.5               
# [9] shiny_1.5.0                 monocle3_0.2.2              SingleCellExperiment_1.8.0  SummarizedExperiment_1.16.1
# [13] DelayedArray_0.12.3         BiocParallel_1.20.1         matrixStats_0.56.0          GenomicRanges_1.38.0       
# [17] GenomeInfoDb_1.22.1         IRanges_2.20.2              S4Vectors_0.24.4            Biobase_2.46.0             
# [21] BiocGenerics_0.32.0        
# 
# loaded via a namespace (and not attached):
#   [1] Rtsne_0.15             colorspace_1.4-1       grr_0.9.5              ellipsis_0.3.1        
# [5] ggridges_0.5.2         mclust_5.4.6           XVector_0.26.0         rstudioapi_0.11       
# [9] farver_2.0.3           leiden_0.3.3           listenv_0.8.0          ggrepel_0.8.2         
# [13] riverplot_0.6          fansi_0.4.1            codetools_0.2-16       splines_3.6.3         
# [17] knitr_1.28             jsonlite_1.6.1         ica_1.0-2              cluster_2.1.0         
# [21] png_0.1-7              uwot_0.1.8             sctransform_0.2.1      compiler_3.6.3        
# [25] httr_1.4.1             assertthat_0.2.1       fastmap_1.0.1          lazyeval_0.2.2        
# [29] cli_2.0.2              later_1.0.0            htmltools_0.5.0        tools_3.6.3           
# [33] rsvd_1.0.3             igraph_1.2.5           gtable_0.3.0           glue_1.4.1            
# [37] GenomeInfoDbData_1.2.2 RANN_2.6.1             reshape2_1.4.4         rappdirs_0.3.1        
# [41] Rcpp_1.0.4.6           vctrs_0.3.0            ape_5.3                nlme_3.1-144          
# [45] iterators_1.0.12       lmtest_0.9-37          xfun_0.24              stringr_1.4.0         
# [49] globals_0.12.5         mime_0.9               lifecycle_0.2.0        irlba_2.3.3           
# [53] future_1.17.0          zlibbioc_1.32.0        MASS_7.3-51.5          zoo_1.8-8             
# [57] scales_1.1.1           doSNOW_1.0.18          promises_1.1.0         RColorBrewer_1.1-2    
# [61] reticulate_1.16        pbapply_1.4-2          gridExtra_2.3          Matrix.utils_0.9.8    
# [65] stringi_1.4.6          foreach_1.5.0          rlang_0.4.6            pkgconfig_2.0.3       
# [69] bitops_1.0-6           lattice_0.20-40        ROCR_1.0-11            purrr_0.3.4           
# [73] labeling_0.3           htmlwidgets_1.5.1      tidyselect_1.1.0       RcppAnnoy_0.0.16      
# [77] plyr_1.8.6             magrittr_1.5           R6_2.4.1               snow_0.4-3            
# [81] pillar_1.4.4           withr_2.2.0            fitdistrplus_1.1-1     survival_3.1-8        
# [85] RCurl_1.98-1.2         tibble_3.0.1           future.apply_1.5.0     tsne_0.1-3            
# [89] crayon_1.3.4           utf8_1.1.4             KernSmooth_2.23-16     plotly_4.9.2.1        
# [93] viridis_0.5.1          grid_3.6.3             data.table_1.12.8      FNN_1.1.3             
# [97] digest_0.6.25          xtable_1.8-4           tidyr_1.1.0            httpuv_1.5.4          
# [101] munsell_0.5.0          viridisLite_0.3.0    

#################################################################################################################
# Panel F generated on another machine 
# note that this setup includes Seurat v2 -- different syntax

# first convert Seurat object to v2 
# cb_dev_s

cb_dev_s_mli = SubsetData(cb_dev_s, ident.use = c('MLI1', "MLI2"), subset.raw = T)
# then limit again to p12-p16 cells
# goes from 1622 to 1370 cells
cb_dev_s_mli = SubsetData(cb_dev_s_mli, subset.name = 'Age', accept.value = c('P12', 'P16'),
                             subset.raw = T)
cb_dev_s_mli@meta.data[['age']] = 'P12-16'
cb_dev_s_mli@meta.data[['joint_clusters']] = cb_dev_s_mli@ident

# generate ints_s object as in Figure 4 walkthrough
# ints_s
mli_s = SubsetData(ints_s, subset.name = 'cluster', accept.value = c('MLI1', 'MLI2'), subset.raw = T)
mli_s@meta.data[['age']] = 'P60'
mli_s@meta.data[['joint_clusters']] = mli_s@meta.data$cluster

mli_merged = MergeSeurat(mli_s, cb_dev_s_mli)
mli_merged = NormalizeData(mli_merged)
mli_merged = SetIdent(mli_merged, ident.use = mli_merged@meta.data$joint_clusters)

# Dotplot of genes with different expression patterns in early and adult interneurons
pdf('extended_data_fig6_panelF.pdf', width = 7, height = 4, useDingbats = F)
SplitDotPlotGG(mli_merged, grouping.var = 'age', genes.plot =  rev(c( 'Fos', 'Junb', 'Fosl2',  'Fam135b',
                                                                      'Sorcs2')), 
               cols.use = c('red', 'blue'), dot.scale = 10, plot.legend = T)
dev.off()


# package versions here
sessionInfo()

# R version 3.5.3 (2019-03-11)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Debian GNU/Linux 9 (stretch)
# 
# Matrix products: default
# BLAS: /usr/lib/libblas/libblas.so.3.7.0
# LAPACK: /usr/lib/lapack/liblapack.so.3.7.0
# 
# locale:
#   [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
# [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                 
# [9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] moments_0.14       png_0.1-7          Matrix.utils_0.9.7 pastecs_1.3.21     Seurat_2.3.4       liger_0.4.1       
# [7] patchwork_0.0.1    Matrix_1.2-17      cowplot_0.9.4      ggplot2_3.1.1     
# 
# loaded via a namespace (and not attached):
#   [1] Rtsne_0.15          colorspace_1.4-1    grr_0.9.5           class_7.3-15        modeltools_0.2-22  
# [6] ggridges_0.5.1      mclust_5.4.3        htmlTable_1.13.1    base64enc_0.1-3     rstudioapi_0.10    
# [11] proxy_0.4-23        npsurv_0.4-0        ggrepel_0.8.0       flexmix_2.3-15      riverplot_0.6      
# [16] bit64_0.9-7         mvtnorm_1.0-10      codetools_0.2-16    splines_3.5.3       R.methodsS3_1.7.1  
# [21] lsei_1.2-0          robustbase_0.93-4   knitr_1.22          Formula_1.2-3       jsonlite_1.6       
# [26] ica_1.0-2           cluster_2.0.7-1     kernlab_0.9-27      R.oo_1.22.0         compiler_3.5.3     
# [31] httr_1.4.0          backports_1.1.4     assertthat_0.2.1    lazyeval_0.2.2      lars_1.2           
# [36] acepack_1.4.1       htmltools_0.3.6     tools_3.5.3         igraph_1.2.4        gtable_0.3.0       
# [41] glue_1.3.1          RANN_2.6.1          reshape2_1.4.3      dplyr_0.8.0.1       Rcpp_1.0.1         
# [46] trimcluster_0.1-2.1 gdata_2.18.0        ape_5.3             nlme_3.1-137        iterators_1.0.10   
# [51] fpc_2.1-11.1        gbRd_0.4-11         lmtest_0.9-36       xfun_0.6            stringr_1.4.0      
# [56] irlba_2.3.3         gtools_3.8.1        DEoptimR_1.0-8      MASS_7.3-51.1       zoo_1.8-5          
# [61] scales_1.0.0        doSNOW_1.0.16       parallel_3.5.3      RColorBrewer_1.1-2  reticulate_1.11.1  
# [66] pbapply_1.4-0       gridExtra_2.3       rpart_4.1-13        segmented_0.5-3.0   latticeExtra_0.6-28
# [71] stringi_1.4.3       foreach_1.4.4       checkmate_1.9.1     caTools_1.17.1.2    boot_1.3-20        
# [76] bibtex_0.4.2        Rdpack_0.10-1       SDMTools_1.1-221    rlang_0.3.4         pkgconfig_2.0.2    
# [81] dtw_1.20-1          prabclus_2.2-7      bitops_1.0-6        lattice_0.20-38     RANN.L1_2.5.2      
# [86] ROCR_1.0-7          purrr_0.3.2         htmlwidgets_1.3     bit_1.1-14          tidyselect_0.2.5   
# [91] plyr_1.8.4          magrittr_1.5        R6_2.4.0            snow_0.4-3          gplots_3.0.1.1     
# [96] Hmisc_4.2-0         pillar_1.4.1        foreign_0.8-71      withr_2.1.2         fitdistrplus_1.0-14
# [101] mixtools_1.1.0      survival_2.43-3     nnet_7.3-12         tsne_0.1-3          tibble_2.1.2       
# [106] crayon_1.3.4        hdf5r_1.1.1         KernSmooth_2.23-15  grid_3.5.3          data.table_1.12.2  
# [111] FNN_1.1.3           metap_1.1           digest_0.6.19       diptest_0.75-7      tidyr_0.8.3        
# [116] R.utils_2.8.0       stats4_3.5.3        munsell_0.5.0      





