# Golgi analysis and annotation

source('../cb_custom_functions.R')
library(liger)
library(Seurat)

golgi_obj = readRDS('golgis_initial.RDS')

golgi_obj@meta.data[['orig.animal']] = sapply(rownames(golgi_obj@meta.data), function(x) {
  strsplit(x, split = '_')[[1]][2]
})
golgi_obj@meta.data[['sex']] = sapply(golgi_obj@meta.data$orig.animal, function(x) {
  substr(x, start = 1, stop = 1)
})

table(golgi_obj@meta.data$orig.animal)
# remove sorted populations if necessary
golgi_obj = SetIdent(golgi_obj, ident.use = golgi_obj@meta.data$orig.animal)
golgi_obj = SubsetData(golgi_obj, ident.remove = c('F0S1', 'F004'), subset.raw = T)

##### convert to liger 
golgi_list = SplitObject(golgi_obj, attribute.1 = 'sex', do.clean = T)

golgi_liger = createLiger(raw.data = list(male = golgi_list$M@raw.data,
                                          female = golgi_list$F@raw.data))
rm(golgi_list)
gc()

golgi_liger = getCBAnimal(golgi_liger)
golgi_liger = getCBRegion(golgi_liger)
golgi_liger = getPercentMito(golgi_liger)

# pct.mito cleaning step -- don't need to remove any in this case
high_mito = which(golgi_liger@cell.data$percent_mito > 0.05)
low_mito = rownames(golgi_liger@cell.data)[-high_mito]
golgi_liger = subsetLiger(golgi_liger, cells.use = low_mito)

golgi_liger = normalize(golgi_liger)

##### gene selection
golgi_liger = selectGenes(golgi_liger, var.thresh = 0.25)
# 2610 genes selected 
spatial_genes = spatial_variable_select_clean_fast(golgi_liger, cell.var = 'region')
# keep track of how many spatial genes were selected
golgi_liger@agg.data[['spatial_genes1']] = spatial_genes
golgi_liger@var.genes = union(spatial_genes, golgi_liger@var.genes)
length(golgi_liger@var.genes)
# 2639

golgi_liger = scaleNotCenter(golgi_liger)

golgi_liger = optimizeALS(golgi_liger, k = 20)
golgi_liger@H.norm = do.call(rbind, golgi_liger@H)
# or run UMAP here
golgi_liger = runTSNE(golgi_liger)

golgi_liger = makeSNN(golgi_liger, k.param = 30, force.recalc = T)
# high resolution for doublet removal
golgi_liger = calcClusters(golgi_liger, resolution = 2)
plotByDatasetAndCluster(golgi_liger)

# Check markers for first round clusters
golgi_seurat = ligerToSeurat2(golgi_liger)
golgi_markers = FindAllMarkers(golgi_seurat, logfc.threshold = 0.3)

##################################################
# round 2
golgi_liger2 = subsetLiger(golgi_liger, clusters.use = setdiff(0:20, c(4, 14, 16, 17, 19, 20, 12, 15, 6)))

golgi_liger2 = normalize(golgi_liger2)

##### gene selection
# var thresh can vary here
golgi_liger2 = selectGenes(golgi_liger2, var.thresh = 0.2)
length(golgi_liger2@var.genes)
# 1757
spatial_genes = spatial_variable_select_clean_fast(golgi_liger2, cell.var = 'region')
# store in all purpose agg.data slot 
length(spatial_genes)
#344
golgi_liger2@agg.data[['spatial_genes1']] = spatial_genes
golgi_liger2@var.genes = union(spatial_genes, golgi_liger2@var.genes)
length(golgi_liger2@var.genes)
# 1780

golgi_liger2 = scaleNotCenter(golgi_liger2)
# k can vary here, lambda should stay 5
# smaller k here
golgi_liger2 = optimizeALS(golgi_liger2, k = 10)
golgi_liger2@H.norm = do.call(rbind, golgi_liger2@H)
golgi_liger2 = runTSNE(golgi_liger2)
golgi_liger2 = makeSNN(golgi_liger2, k.param = 30, force.recalc = T)
golgi_liger2 = calcClusters(golgi_liger2, resolution = 0.1, force.recalc = T)
plotByDatasetAndCluster(golgi_liger2)

golgi_liger2 = runUMAP(golgi_liger2)
plotByDatasetAndCluster(golgi_liger2)

# check final markers
golgi_seurat = ligerToSeurat2(golgi_liger2)
golgi_markers2 = FindAllMarkers(golgi_seurat, logfc.threshold = 0.3)

# save final object
saveRDS(golgi_liger2, 'golgi_liger.RDS')
saveRDS(golgi_seurat, 'golgi_seurat.RDS')
