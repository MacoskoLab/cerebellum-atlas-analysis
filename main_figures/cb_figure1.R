#########################################
# Figure 1 

library(Seurat)  # v2.3.4
library(Matrix)

# Read in full dataset (as mtx or RDS file)
# Raw counts can be downloaded from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE165371
# (cb_adult_mouse)
counts_cb = readMM('cb_adult_mouse.mtx')
barcodes_cb = readLines('cb_adult_mouse_barcodes.txt')
genes_cb = readLines('cb_adult_mouse_genes.txt')

row.names(counts_cb) = genes_cb
colnames(counts_cb) = barcodes_cb

# full_cb = readRDS('full_cb_mouse.RDS')

# Create Seurat object if reading in raw data
full_cb = CreateSeuratObject(raw.data = counts_cb)

full_cb = NormalizeData(full_cb)

# basic dimensionality reduction for course cell type visualization
full_cb = FindVariableGenes(full_cb, do.plot = F)  # top 2000 highly variable genes
full_cb = ScaleData(full_cb, genes.use = full_cb@var.genes)

full_cb = RunPCA(full_cb, pcs.compute = 50, do.print = F)
full_cb = RunUMAP(full_cb, dims.use = 1:25)

# read in cluster metadata (available in data directory)
cluster_metadata = read.csv('cluster_metadata.csv', header = T, row.names = 1)
full_cb@meta.data = cluster_metadata


# Figure 1B
# UMAP by cluster
png("Figure1_panelB.png", units="in", width=6, height=5, res=400)
DimPlot(full_cb, reduction.use = 'umap', do.label = F, group.by = "cluster",
        no.legend = T, no.axes = T, pt.size = 0.5)
dev.off()


# Figure 1C

library(Matrix.utils)
library(dendsort)

# read in metadata corresponding to object with GCs downsampled to 60k
cluster_metadata_ds = read.csv('cluster_metadata_ds.csv', header = T, row.names = 1)
# read in union of variable genes determined in cluster subanalyses (used for dendrogram)
var_gene_union = readLines('var_genes_union.txt')

ds_gran_cb = SubsetData(full_cb, cells.use = row.names(cluster_metadata_ds), subset.raw = T)
ds_gran_cb = SetIdent(ds_gran_cb, ident.use = ds_gran_cb@meta.data$subcluster)

# generate dendrogram
dge.transpose = Matrix::t(ds_gran_cb@raw.data)
metacells = aggregate.Matrix(dge.transpose, groupings = ds_gran_cb@ident, fun="sum")
metacells.norm = Matrix::t(Matrix.column_norm(Matrix::t(metacells)))
metacells.scaled = scale(as.matrix(metacells.norm))
colnames(metacells.scaled) = row.names(ds_gran_cb@raw.data)

# subset to just relevant genes
m = metacells.scaled[, var_gene_union]

dist.mat = dist(m)

hier = dendsort(hclust(dist.mat, method = 'mcquitty'))

# Figure 1c dendrogram only
pdf("Figure1_panelC_dendrogram.pdf",useDingbats = F, width = 5, height = 4)
plot(hier, labels=F, hang=-1)
dev.off()

# Genes for dotplot
genes_plot = rev(c('Ppp1r17', 'Gabra6', 'Eomes', 'Lypd6', 'Acvr1c', 'Klhl1', 'Lgi2', 'Gdf10', 'Aqp4', 
                   'Mobp', 'Ppfibp1', 'Dcn', 'Kcnj8', 'Ttr', 'Mrc1', 'C1qa', 'Flt1',
                   'Foxj1'))

ordered_ident = factor(ds_gran_cb@meta.data$subcluster, levels = row.names(m)[hier$order])
ds_gran_cb@ident = ordered_ident
ds_gran_cb = NormalizeData(ds_gran_cb)
dotplot = Seurat::DotPlot(ds_gran_cb, genes.plot = genes_plot, dot.min = 0.08, scale.by = 'radius',
                          x.lab.rot = T, dot.scale = 6, plot.legend = T, do.return = T)
dotplot = dotplot + theme(axis.text.x=element_text(angle=45,hjust = 0.95, vjust=0.95))


# Figure 1C dotplot 
pdf('Figure1_panelC_dotplot.pdf', useDingbats = F,
    width = 9, height = 8)
dotplot
dev.off()


################################
# Extra functions required 

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


