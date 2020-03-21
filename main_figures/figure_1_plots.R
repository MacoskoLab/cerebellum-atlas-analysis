# FIGURE 1 -- DENDROGRAM AND DOTPLOT

full_seurat_gran_ds = readRDS('/home/vkozarev/final_redo/seurat/full_seurat_gran_ds_v3.RDS')
var_gene_union = readRDS('/home/vkozarev/final_redo/var_gene_union_v3.RDS')

library(Matrix)
library(Matrix.utils)
dge.transpose = Matrix::t(full_seurat_gran_ds@raw.data)
metacells<-aggregate.Matrix(dge.transpose,groupings = full_seurat_gran_ds@ident,fun="sum")
metacells.norm = Matrix::t(Matrix.column_norm(Matrix::t(metacells)))
metacells.scaled = scale(as.matrix(metacells.norm))
colnames(metacells.scaled) = rownames(full_seurat_gran_ds@raw.data)

m = metacells.scaled[,var_gene_union]

dist.mat = dist(m)
library(dendsort)
# this corresponds to weighted pair group clustering (wpgma)
h3 = dendsort(hclust(dist.mat,method = 'mcquitty'))

# dotplot genes
dp_genes_final_v2 = rev(c('Ppp1r17', 'Gabra6', 'Eomes', 'Lypd6', 'Prkcd', 'Klhl1', 'Lgi2', 'Gdf10', 
                          'Aqp4', 
                          'Mobp', 'Ppfibp1', 'Dcn', 'Kcnj8', 'Ttr', 'Mrc1', 'C1qa', 'Flt1',
                          'Foxj1'))

ordered_ident = factor(full_seurat_gran_ds@meta.data$liger_ident, levels = row.names(m)[h3$order])
full_seurat_gran_ds@ident = ordered_ident
full_seurat_gran_ds = NormalizeData(full_seurat_gran_ds)
dotplot = Seurat::DotPlot(full_seurat_gran_ds, genes.plot = dp_genes_final_v2, dot.min = 0.08, 
                          scale.by = 'radius',
                          x.lab.rot = T, dot.scale = 6, plot.legend = T, do.return = T)
dotplot = dotplot+ theme(axis.text.x=element_text(angle=45,hjust = 0.95, vjust=0.95))

# generate dotplot
pdf('/home/vkozarev/final_redo/full_dotplot_coarse_v4.pdf', useDingbats = F,
    width = 9, height = 8)
dotplot
dev.off()

# generate dendrogram without leaf labels
pdf("/home/vkozarev/final_redo/full_Dendrogram_var_gene_union_redo3.pdf",useDingbats = F,
    width = 5, height = 4)
plot(h3,labels=F,hang=-1)
dev.off()