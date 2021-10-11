# Reproducing Extended Data Fig 2 panels
library(ggplot2)
library(Seurat)
# v2.3.4
library(liger)
# v0.4.1
source('../cb_custom_functions.R')
source('../src/fig2_utils.R')

# Requires use of individual analysis objects or tSNE/UMAP coordinates

purk_s = readRDS('purkinje_seurat.RDS')
purk_plot = TSNEPlot(purk_s, no.axes = T, 
                     no.legend = T, pt.size = 0.05, do.label = F, do.return = T)

ubc_s = readRDS('ubc_seurat.RDS')
ubc_plot = TSNEPlot(ubc_s, no.axes = T, 
                    no.legend = T, pt.size = 0.1, do.label = F, do.return = T)

ints_s = readRDS('interneurons_seurat.RDS')
ints_plot = TSNEPlot(ints_s, no.axes = T, 
                     no.legend = T, pt.size = 0.05, do.label = F, do.return = T)

gran_s = readRDS('granule_seurat.RDS')
# note that this is actually a UMAP -- just called tSNE here for convenience
gran_plot = TSNEPlot(gran_s, no.axes = T, 
                     no.legend = T, pt.size = 0.05, do.label = F, do.return = T)

# png("ExtendedDataFig2_row1.png", units="in", width=16, height=4, res=300)
plot_grid(purk_plot, gran_plot, ubc_plot, ints_plot, nrow = 1)
# dev.off()

golgi_s = readRDS('golgi_seurat.RDS')
golgi_plot = TSNEPlot(golgi_s, no.axes = T, 
                      no.legend = T, pt.size = 0.05, do.label = F, do.return = T)

bergmann_s = readRDS('bergmann_seurat.RDS')
bergmann_plot = TSNEPlot(bergmann_s, no.axes = T, 
                         no.legend = T, pt.size = 0.05, do.label = F, do.return = T)
astro_s = readRDS('astrocyte_seurat.RDS')
astro_plot = TSNEPlot(astro_s, no.axes = T, 
                      no.legend = T, pt.size = 0.05, do.label = F, do.return = T)

oligo_s = readRDS('oligo_seurat.RDS')
oligo_plot = TSNEPlot(oligo_s, no.axes = T, 
                      no.legend = T, pt.size = 0.05, do.label = F, do.return = T)

# png("ExtendedDataFig2_row2.png", units="in", width=16, height=4, res=300)
plot_grid(golgi_plot, bergmann_plot, astro_plot, oligo_plot, nrow = 1)
# dev.off()

endo_s = readRDS('endothelial_seurat.RDS')
mural_stalk_endo = TSNEPlot(endo_s, no.axes = T, 
                       no.legend = T, pt.size = 0.05, do.label = F, do.return = T)

# choroid
choroid_s = readRDS('choroid_seurat.RDS')
choroid_plot = TSNEPlot(choroid_s, no.axes = T, 
                        no.legend = T, pt.size = 0.2, do.label = F, do.return = T)

micro_s = readRDS('microglia_seurat.RDS')
micro_plot = TSNEPlot(micro_s, no.axes = T, 
                      no.legend = T, pt.size = 0.2, do.label = F, do.return = T)

# png("ExtendedDataFig2_row3.png", units="in", width=12, height=4, res=300)
plot_grid(mural_stalk_endo, choroid_plot, micro_plot, nrow = 1)
# dev.off()

# Dotplots for later panels 
oligo_s = NormalizeData(oligo_s)
# using order from dendrogram in Fig 1 -- generate hier using that script
oligo_s@ident = factor(oligo_s@ident, levels = row.names(m)[hier$order][30:39])
# pdf('ExtendedDataFig2_panelL.pdf', useDingbats = F, 
#     height = 5, width = 8)
# replaced Il33 with Palm2
DotPlot2(oligo_s, (rev(c('Hcn2', 'Il33', 'Synpr', 'Klk6', 'Sepp1', 'Ptma', 'Abr', 'Gpr17', 'Pdgfra', 'Top2a'))), 
         dot.scale = 10,scale.min = 0.1,plot.legend = T, scale.by = 'size',
         x.lab.rot = 45)
# dev.off()

golgi_s = NormalizeData(golgi_s)
# pdf('ExtendedDataFig2_panelM.pdf', useDingbats = F, 
#     height = 2, width = 6)

DotPlot2(golgi_s, ((c('Sst', 'Nxph1', 'Gjd2', 'Sorcs3', 'Grm2'))), 
         dot.scale = 10,scale.min = 0.1,plot.legend = T, scale.by = 'size',
         x.lab.rot = 45)
# dev.off()


astro_s = NormalizeData(astro_s)
# pdf('ExtendedDataFig2_panelN.pdf', useDingbats = F, 
#     height = 2, width = 4)

DotPlot2(astro_s, (rev(c('Gpc5', 'Gfap'))), 
         dot.scale = 10,scale.min = 0.1,plot.legend = T, scale.by = 'size',
         x.lab.rot = 45)
# dev.off()

