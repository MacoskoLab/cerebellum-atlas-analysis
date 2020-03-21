# final liger script for joint astrocyte/bergmann
source('../cb_custom_functions.R')
library(liger)
astro_berg_seurat = readRDS('AstroBerg_seurat.rds')
astro_berg_seurat@meta.data[['orig.animal']] = sapply(rownames(astro_berg_seurat@meta.data), function(x) {
  strsplit(x, split = '_')[[1]][2]
})
astro_berg_seurat@meta.data[['sex']] = sapply(astro_berg_seurat@meta.data$orig.animal, function(x) {
  substr(x, start = 1, stop = 1)
})

# check whether any from sorted populations ended up here -- F0S1 F004
table(astro_berg_seurat@meta.data$orig.animal)
# these are clean, but if not, could remove with this step
astro_berg_seurat = SetIdent(astro_berg_seurat, ident.use = astro_berg_seurat@meta.data$orig.animal)
astro_berg_seurat = SubsetData(astro_berg_seurat, ident.remove = c('F0S1'), subset.raw = T)

ab_list = SplitObject(astro_berg_seurat, attribute.1 = 'sex', do.clean = T)

ab_liger = createLiger(raw.data = list(male = ab_list$M@raw.data,
                                       female = ab_list$F@raw.data))
rm(ab_list)
gc() 

ab_liger = getCBAnimal(ab_liger)
ab_liger = getCBRegion(ab_liger)
ab_liger = getPercentMito(ab_liger)

# pct.mito cleaning step
high_mito = which(ab_liger@cell.data$percent_mito > 0.05)
low_mito = rownames(ab_liger@cell.data)[-high_mito]
ab_liger = subsetLiger(ab_liger, cells.use = low_mito)

ab_liger = normalize(ab_liger)

##### gene selection
# var thresh can vary here
ab_liger = selectGenes(ab_liger, var.thresh = 0.12)
# 2456 genes selected 
spatial_genes = spatial_variable_select_clean_fast(ab_liger, cell.var = 'region')
# keep track of how many spatial genes were selected
# store in all purpose agg.data slot 
length(spatial_genes)
#1383
ab_liger@agg.data[['spatial_genes1']] = spatial_genes
ab_liger@var.genes = union(spatial_genes, ab_liger@var.genes)
length(ab_liger@var.genes)
# 2875

ab_liger = scaleNotCenter(ab_liger)
# k can vary here, lambda should stay 5
ab_liger = optimizeALS(ab_liger, k = 25)
ab_liger@H.norm = do.call(rbind, ab_liger@H)
# or run UMAP here
ab_liger = runTSNE(ab_liger)
# let's keep k.param 30 unless absolutely necessary
# resolution can vary 
ab_liger = clusterLouvainJaccard2(ab_liger, k.param = 30, resolution = 2)
plotByDatasetAndCluster(ab_liger)
# note that this first round is mostly to separate astrocytes and bergmanns effectively and 
# remove any doublets between them while ensuring they're not too similar to either population

plotFeature2(ab_liger, 'percent_mito', by.dataset = F)
plotGene(ab_liger, 'Rbfox3', use.raw = T, plot.by = 'none', max.clip = 30, clip.absolute = T)
plotGene(ab_liger, 'Ttr', use.raw = T, plot.by = 'none', max.clip = 30, clip.absolute = T)
plotGene(ab_liger, 'Gdf10', use.raw = T, plot.by = 'none', max.clip = 30, clip.absolute = T)
plotGene(ab_liger, 'Bcan', use.raw = T, plot.by = 'none', max.clip = 30, clip.absolute = T)
plotGene(ab_liger, 'Shroom3', use.raw = T, plot.by = 'none', max.clip = 30, clip.absolute = T)
plotGene(ab_liger, 'Dcn', use.raw = T, plot.by = 'none', max.clip = 30, clip.absolute = T)
plotGene(ab_liger, 'Acta2', use.raw = T, plot.by = 'none', max.clip = 30, clip.absolute = T)
plotGene(ab_liger, 'Gad1', use.raw = T, plot.by = 'none', max.clip = 30, clip.absolute = T)
plotGene(ab_liger, 'Calb1', use.raw = T, plot.by = 'none', max.clip = 30, clip.absolute = T)
plotGene(ab_liger, 'Parp14', use.raw = T, plot.by = 'none', max.clip = 30, clip.absolute = T)

pdf('cluster_hold.pdf')
plotByDatasetAndCluster(ab_liger)
dev.off()
# clusters to remove: 7, 14, 16 (all 3 doublets), 19 (Purk contamination) 21 (astro/berg doublet), 23 (shroom3), 30 (doublet),
# 31 (doublet), 26 (Purk contamination), 25 (Robo2/Meis2 contamination)
ab_seurat = ligerToSeurat(ab_liger)
DimPlot(ab_seurat, reduction.use = 'tsne', cells.highlight = WhichCells(ab_seurat, ident = 19))
DimPlot(ab_seurat, reduction.use = 'tsne', cells.highlight = WhichCells(ab_seurat, ident = 20))
DimPlot(ab_seurat, reduction.use = 'tsne', cells.highlight = WhichCells(ab_seurat, ident = 21))
DimPlot(ab_seurat, reduction.use = 'tsne', cells.highlight = WhichCells(ab_seurat, ident = 15))
markers27 = FindMarkers(ab_seurat, ident.1 = 27, logfc.threshold = 0.3, max.cells.per.ident = 1000)
# clusters to keep for astrocyte object 
# 0, 2, 3, 5, 6, 11, 18, 22, 24, 29
# clusters to keep for bergmann object
# 1, 4, 8, 9, 10, 12, 13, 15, 17, 20, 27, 28, 32

astro_liger = subsetLiger(ab_liger, clusters.use = c(0, 2, 3, 5, 6, 11, 18, 22, 24, 29))
bergmann_liger = subsetLiger(ab_liger, clusters.use = c(1, 4, 8, 9, 10, 12, 13, 15, 17, 20, 27, 28, 32))
saveRDS(astro_liger, '../bergmann_astro/astro_liger.RDS')
saveRDS(bergmann_liger, '../bergmann_astro/bergmann_liger.RDS')