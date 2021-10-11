# Reproducing Extended Data Fig 3 panels

# Some utility functions
library(ggplot2)
source('../src/fig2_utils.R')

# function to generate pseudo replicates across regions 
split_rep1_rep2_meta = function(seurat_meta) {
  seurat_meta[['individual']] = sapply(row.names(seurat_meta), function(x) {
    strsplit(x, '_')[[1]][2]
  })
  cross = table(seurat_meta$individual, seurat_meta$region)
  
  # rep1 is the max
  rep1_inds = c()
  # rep 2 is the second highest
  rep2_inds = c()
  
  for (reg in colnames(cross)) {
    sorted = sort(cross[,reg], decreasing = T)
    if (length(which(sorted == 0)) > 4) {
      print(paste0('Region ', reg, ' has representation from fewer than 2 replicates.'))
    }
    rep1_inds = c(rep1_inds, names(sorted)[1])
    rep2_inds = c(rep2_inds, names(sorted)[2])
  }
  names(rep1_inds) = names(rep2_inds) = colnames(cross)
  
  seurat_meta[['rep_desig']] = 'NA'
  for (lobe in names(rep1_inds)) {
    rep1_ind_val = rep1_inds[lobe]
    seurat_meta[seurat_meta$individual == rep1_ind_val & seurat_meta$region == lobe, 'rep_desig'] = 'rep1'
  } 
  for (lobe in names(rep2_inds)) {
    rep2_ind_val = rep2_inds[lobe]
    seurat_meta[seurat_meta$individual == rep2_ind_val & seurat_meta$region == lobe, 'rep_desig'] = 'rep2'
  } 
  
  return(seurat_meta)
}

# Read in metadata corresponding to object with granule cells downsampled to 60k (as in Figure 2)
cluster_metadata_ds = read.csv('../data/cluster_metadata_ds.csv', header = T, row.names = 1,
                               stringsAsFactors = F)

# split up by analyses 
cluster_metadata_ds[['analysis_cluster']] = sapply(cluster_metadata_ds$cluster, function(x) {
  if (x %in% c('OPC', 'ODC')) {
    return('OPC/ODC')
  } else if (x %in% c('PLI', 'MLI1', 'MLI2')) {
    return('MLI/PLI')
  } else if (x %in% c('Macrophages', 'Microglia')){
    return('Microglia')
  } else if (x %in% c('Choroid', 'Ependymal')) {
    return('Choroid')
  } else if (x %in% c('Dcn_Fibroblasts', 'Endothelial_mural', 'Endothelial_stalk')) {
    return('Endothelial/Fibroblasts')
  } else {
    return(x)
  }
})

cluster_metadata_ds[['lob_comp_cluster']] = sapply(cluster_metadata_ds$subcluster, function(x) {
  if (x %in% c('PLI_1', 'PLI_2', 'PLI_3')) {
    return('PLI')
  } else {
    return(x)
  }
})
cluster_metadata_ds$lob_comp_cluster = factor(cluster_metadata_ds$lob_comp_cluster)
cluster_metadata_ds$analysis_cluster = factor(cluster_metadata_ds$analysis_cluster)

## Panel A

celltypes_test = c('Purkinje', 'Granule', 'OPC/ODC', 'Golgi', 'UBC',
                   'MLI/PLI', 'Bergmann', 'Astrocyte' )

celltype_cors = c()
celltype_names = c()
for (celltype in celltypes_test) {
  celltype_sub = cluster_metadata_ds[cluster_metadata_ds$analysis_cluster == celltype,]
  celltype_sub$lob_comp_cluster = droplevels(celltype_sub$lob_comp_cluster)
  
  # add replicate designation
  celltype_sub2 = split_rep1_rep2_meta(celltype_sub)
  
  spat_enrich_cell = rep_spatial_enrichments_meta(celltype_sub2)
  celltype_rep_corrs = sapply(row.names(spat_enrich_cell[[1]]), function(x) {
    if (!(x %in% row.names(spat_enrich_cell$rep1) & x %in% row.names(spat_enrich_cell$rep2))) {
      print(paste0('Cluster ', x, ' only present in one replicate.'))
      NA
    } else {
      cor(spat_enrich_cell$rep1[x,], spat_enrich_cell$rep2[x,])
    }
  })
  celltype_cors = c(celltype_cors, celltype_rep_corrs)
  celltype_names = c(celltype_names, rep(celltype, length(celltype_rep_corrs)))
}


full_rep_corrs_df = data.frame(pearson = celltype_cors,
                               subtype = names(celltype_cors),
                               celltype = celltype_names)

full_rep_corrs_df = full_rep_corrs_df[order(full_rep_corrs_df$pearson),]
full_rep_corrs_df$subtype = factor(full_rep_corrs_df$subtype, levels = full_rep_corrs_df$subtype)
full_rep_corrs_df$top_stat = 'Not'

# pdf('ExtendedDataFig3_panelA.pdf', useDingbats = F,
#     width = 7, height = 5)
ggplot(full_rep_corrs_df, aes(x = subtype, y = pearson, col = celltype)) + geom_point() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 0.95)) + 
  geom_hline(yintercept = 0.85, linetype = 'dashed')
# dev.off()


# Panel D 

# read in full cluster metadata (available in data directory)
cluster_metadata = read.csv('../data/cluster_metadata.csv', header = T, row.names = 1)

# Prepare cerebellar cartoon polygons for plotting later panels
# Lobule coordinates available in data directory
coords_full = c('coordsI_half.csv',  'coordsII_half.csv', 'coordsIII_half.csv', 
                'coordscul_half.csv', 'coordsVI_half.csv', 'coordsVII_half.csv',
                'coordsVIII_half.csv', 'coordsIX_half.csv', 'coordsX_half.csv',
                'coordsan12.csv', 'coordsan22.csv', 'coordsprm2.csv', 
                'coordssim2.csv', 'coordscop2.csv', 'coordsf2.csv', 
                'coordspf2.csv')

coords_full = paste0('../data/coords/', coords_full)

coords_names = c('I', 'II', 'III', 'CUL', 'VI', 'VII', 'VIII', 'IX', 'X',
                 'AN1', 'AN2', 'PRM', 'SIM', 'COP', 'F', 'PF')

SPs = generate_SPs(coords_full, coords_names)

# pdf('ExtendedDataFig3_panelD.pdf', useDingbats = F)
spatial_enrichment_map_meta(cluster_metadata, cluster = 'UBC', sps = SPs, do.print = F, color_divergent = T)
# dev.off()









