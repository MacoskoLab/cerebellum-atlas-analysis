# Reproducing Extended Data Fig 4 panels 

library(ggplot2)
library(Seurat)
# v2.3.4

# Additional plotting functions 

library(dplyr)
library(reshape2)
# Based on DoHeatmap from Seurat with additional color functionality 
DoHeatmap2 = function (object, data.use = NULL, use.scaled = TRUE, cells.use = NULL, 
                       genes.use = NULL, disp.min = -2.5, disp.max = 2.5, group.by = "ident", 
                       group.order = NULL, draw.line = TRUE, col.low = "#FF00FF", 
                       col.mid = "#000000", col.high = "#FFFF00", slim.col.label = FALSE, 
                       remove.key = FALSE, rotate.key = FALSE, title = NULL, cex.col = 10, 
                       cex.row = 10, group.label.loc = "bottom", group.label.rot = FALSE, 
                       group.cex = 15, group.spacing = 0.15, assay.type = "RNA", 
                       do.plot = TRUE) 
{
  if (is.null(x = data.use)) {
    if (use.scaled) {
      data.use <- GetAssayData(object, assay.type = assay.type, 
                               slot = "scale.data")
    }
    else {
      data.use <- GetAssayData(object, assay.type = assay.type, 
                               slot = "data")
    }
  }
  cells.use <- SetIfNull(x = cells.use, default = object@cell.names)
  cells.use <- intersect(x = cells.use, y = colnames(x = data.use))
  if (length(x = cells.use) == 0) {
    stop("No cells given to cells.use present in object")
  }
  genes.use <- SetIfNull(x = genes.use, default = rownames(x = data.use))
  genes.use <- intersect(x = genes.use, y = rownames(x = data.use))
  if (length(x = genes.use) == 0) {
    stop("No genes given to genes.use present in object")
  }
  if (is.null(x = group.by) || group.by == "ident") {
    cells.ident <- object@ident[cells.use]
  }
  else {
    cells.ident <- factor(x = FetchData(object = object, 
                                        cells.use = cells.use, vars.all = group.by)[, 1])
    names(x = cells.ident) <- cells.use
  }
  cells.ident <- factor(x = cells.ident, labels = intersect(x = levels(x = cells.ident), 
                                                            y = cells.ident))
  data.use <- data.use[genes.use, cells.use, drop = FALSE]
  if ((!use.scaled)) {
    data.use = as.matrix(x = data.use)
    if (disp.max == 2.5) 
      disp.max = 10
  }
  data.use <- MinMax(data = data.use, min = disp.min, max = disp.max)
  data.use <- as.data.frame(x = t(x = data.use))
  data.use$cell <- rownames(x = data.use)
  colnames(x = data.use) <- make.unique(names = colnames(x = data.use))
  data.use <- data.use %>% melt(id.vars = "cell")
  names(x = data.use)[names(x = data.use) == "variable"] <- "gene"
  names(x = data.use)[names(x = data.use) == "value"] <- "expression"
  data.use$ident <- cells.ident[data.use$cell]
  if (!is.null(group.order)) {
    if (length(group.order) == length(levels(data.use$ident)) && 
        all(group.order %in% levels(data.use$ident))) {
      data.use$ident <- factor(data.use$ident, levels = group.order)
    }
    else {
      stop("Invalid group.order")
    }
  }
  breaks <- seq(from = min(data.use$expression), to = max(data.use$expression), 
                length = length(x = PurpleAndYellow()) + 1)
  data.use$gene <- with(data = data.use, expr = factor(x = gene, 
                                                       levels = rev(x = unique(x = data.use$gene))))
  data.use$cell <- with(data = data.use, expr = factor(x = cell, 
                                                       levels = cells.use))
  if (rotate.key) {
    key.direction <- "horizontal"
    key.title.pos <- "top"
  }
  else {
    key.direction <- "vertical"
    key.title.pos <- "left"
  }
  heatmap <- ggplot(data = data.use, mapping = aes(x = cell, 
                                                   y = gene, fill = expression)) + 
    geom_tile() + viridis::scale_fill_viridis(name = "Expression", discrete = F, option = 'C',
                                              guide = guide_colorbar(direction = key.direction, title.position = key.title.pos)) + 
    scale_y_discrete(position = "right", labels = rev(genes.use)) + 
    theme(axis.line = element_blank(), axis.title.y = element_blank(), 
          axis.ticks.y = element_blank(), strip.text.x = element_text(size = group.cex), 
          axis.text.y = element_text(size = cex.row), axis.text.x = element_text(size = cex.col), 
          axis.title.x = element_blank())
  if (slim.col.label) {
    heatmap <- heatmap + theme(axis.title.x = element_blank(), 
                               axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
                               axis.line = element_blank(), axis.title.y = element_blank(), 
                               axis.ticks.y = element_blank())
  }
  else {
    heatmap <- heatmap + theme(axis.text.x = element_text(angle = 90))
  }
  if (!is.null(x = group.by)) {
    if (group.label.loc == "top") {
      switch <- NULL
    }
    else {
      switch <- "x"
    }
    heatmap <- heatmap + facet_grid(facets = ~ident, drop = TRUE, 
                                    space = "free", scales = "free", switch = switch, 
    ) + scale_x_discrete(expand = c(0, 0), drop = TRUE)
    if (draw.line) {
      panel.spacing <- unit(x = group.spacing, units = "lines")
    }
    else {
      panel.spacing <- unit(x = 0, units = "lines")
    }
    heatmap <- heatmap + theme(strip.background = element_blank(), 
                               panel.spacing = panel.spacing)
    if (group.label.rot) {
      heatmap <- heatmap + theme(strip.text.x = element_text(angle = 90))
    }
  }
  if (remove.key) {
    heatmap <- heatmap + theme(legend.position = "none")
  }
  if (!is.null(x = title)) {
    heatmap <- heatmap + labs(title = title)
  }
  if (do.plot) {
    heatmap
  }
  return(heatmap)
}

##################################################################################
library(tidyr)
# Based on Seurat's SplitDotPlotGG function, allows for plotting of pct.exp for non-gene features
# (columns in meta.data)
SplitDotPlot2 = function (object, grouping.var, genes.plot = NULL, gene.groups, features.plot = NULL,
                          cols.use = c("blue", 
                                       "red"), col.min = -2.5, col.max = 2.5, dot.min = 0, dot.scale = 6, 
                          group.by, plot.legend = FALSE, do.return = FALSE, x.lab.rot = FALSE) 
{
  if (!missing(x = group.by)) {
    object <- SetAllIdent(object = object, id = group.by)
  }
  grouping.data <- FetchData(object = object, vars.all = grouping.var)[names(x = object@ident), 
                                                                       1]
  ncolor <- length(x = cols.use)
  ngroups <- length(x = unique(x = grouping.data))
  if (ncolor < ngroups) {
    stop(paste("Not enough colors supplied for number of grouping variables. Need", 
               ngroups, "got", ncolor, "colors"))
  }
  else if (ncolor > ngroups) {
    cols.use <- cols.use[1:ngroups]
  }
  idents.old <- levels(x = object@ident)
  idents.new <- paste(object@ident, grouping.data, sep = "_")
  colorlist <- cols.use
  names(x = colorlist) <- levels(x = grouping.data)
  object@ident <- factor(x = idents.new, levels = unlist(x = lapply(X = idents.old, 
                                                                    FUN = function(x) {
                                                                      lvls <- list()
                                                                      for (i in seq_along(along.with = levels(x = grouping.data))) {
                                                                        lvls[[i]] <- paste(x, levels(x = grouping.data)[i], 
                                                                                           sep = "_")
                                                                      }
                                                                      return(unlist(x = lvls))
                                                                    })), ordered = TRUE)
  if (!is.null(genes.plot)) {
    data.to.plot <- data.frame(FetchData(object = object, vars.all = genes.plot))
  } else {
    genes.plot = c()
  }
  
  if (!is.null(features.plot)) {
    # for some reason, a fetch data call here screws it up
    features_data = object@meta.data[, features.plot]
    if (!exists('data.to.plot')) {
      data.to.plot = features_data
    } else {
      data.to.plot = cbind(data.to.plot, features_data)
    }
    genes.plot = c(features.plot, genes.plot)
  }
  
  data.to.plot$cell <- rownames(x = data.to.plot)
  data.to.plot$id <- object@ident
  data.to.plot <- data.to.plot %>% gather(key = genes.plot, 
                                          value = expression, -c(cell, id))
  data.to.plot <- data.to.plot %>% group_by(id, genes.plot) %>% 
    summarize(avg.exp = ExpMean(x = expression), pct.exp = Seurat:::PercentAbove(x = expression, 
                                                                                 threshold = 0))
  data.to.plot <- data.to.plot %>% ungroup() %>% group_by(genes.plot) %>% 
    mutate(avg.exp = scale(x = avg.exp)) %>% mutate(avg.exp.scale = as.numeric(x = cut(x = MinMax(data = avg.exp, 
                                                                                                  max = col.max, min = col.min), breaks = 20)))
  
  
  data.to.plot <- data.to.plot %>% separate(col = id, into = c("ident1", 
                                                               "ident2"), sep = "_") %>% rowwise() %>% mutate(palette.use = colorlist[[ident2]], 
                                                                                                              ptcolor = colorRampPalette(colors = c("grey", palette.use))(20)[avg.exp.scale]) %>% 
    unite("id", c("ident1", "ident2"), sep = "_")
  data.to.plot$genes.plot <- factor(x = data.to.plot$genes.plot, 
                                    levels = rev(x = sub(pattern = "-", replacement = ".", 
                                                         x = genes.plot)))
  data.to.plot$pct.exp[data.to.plot$pct.exp < dot.min] <- NA
  data.to.plot$id <- factor(x = data.to.plot$id, levels = levels(object@ident))
  palette.use <- unique(x = data.to.plot$palette.use)
  if (!missing(x = gene.groups)) {
    names(x = gene.groups) <- genes.plot
    data.to.plot <- data.to.plot %>% mutate(gene.groups = gene.groups[genes.plot])
  }
  data.to.plot$pct.exp <- data.to.plot$pct.exp * 100
  p <- ggplot(data = data.to.plot, mapping = aes(x = genes.plot, 
                                                 y = id)) + geom_point(mapping = aes(size = pct.exp, color = ptcolor)) + 
    scale_radius(range = c(0, dot.scale)) + scale_color_identity() + 
    theme(axis.title.x = element_blank(), axis.title.y = element_blank())
  if (!missing(x = gene.groups)) {
    p <- p + facet_grid(facets = ~gene.groups, scales = "free_x", 
                        space = "free_x", switch = "y") + theme(panel.spacing = unit(x = 1, 
                                                                                     units = "lines"), strip.background = element_blank(), 
                                                                strip.placement = "outside")
  }
  if (x.lab.rot) {
    p <- p + theme(axis.text.x = element_text(angle = 90, 
                                              vjust = 0.5))
  }
  if (!plot.legend) {
    p <- p + theme(legend.position = "none")
  }
  else if (plot.legend) {
    plot.legend <- cowplot::get_legend(plot = p)
    palettes <- list()
    for (i in seq_along(along.with = colorlist)) {
      palettes[[names(colorlist[i])]] <- colorRampPalette(colors = c("grey", 
                                                                     colorlist[[i]]))(20)
    }
    gradient.legends <- mapply(FUN = Seurat:::GetGradientLegend, palette = palettes, 
                               group = names(x = palettes), SIMPLIFY = FALSE, USE.NAMES = FALSE)
    p <- p + theme(legend.position = "none")
    legends <- cowplot::plot_grid(plotlist = gradient.legends, 
                                  plot.legend, ncol = 1, rel_heights = c(1, rep.int(x = 0.5, 
                                                                                    times = length(x = gradient.legends))), scale = rep(0.5, 
                                                                                                                                        length(gradient.legends)), align = "hv")
    p <- cowplot::plot_grid(p, legends, ncol = 2, rel_widths = c(1, 
                                                                 0.3), scale = c(1, 0.8))
  }
  suppressWarnings(print(p))
  if (do.return) {
    return(p)
  }
}
################################################################################
# Panel A 
# Use joint object from SCP, download to data folder

ubc_joint = readRDS('../data/mouse_human_joint_ubc_obj.RDS')

# note labeling error in caption of figure -- all lower dimensional representations are tSNE
# except for granule, which is UMAP

# png('ExtendedDataFig4_panelA', units="in", width=5, height=3.5, res=300)
DimPlot(ubc_joint, reduction.use = 'tsne', group.by = 'species', no.axes = T, pt.size=0.1)
# dev.off()


# Panel B
ubc_mouse = SubsetData(ubc_joint, subset.name = 'species', accept.value = 'mouse', 
                       subset.raw = T, do.scale = T, do.center = F)
ubc_hum = SubsetData(ubc_joint, subset.name = 'species', accept.value = 'human', 
                     subset.raw = T, do.scale = T, do.center = F)

# gene markers of traditionally on/off UBCs and those loading highly on relevant iNMF factors
genes_ubc = c('Dgkg', 'Dgkb', 'Trpc3', 'Slc24a2', 'Esrrg','Grm1',  'Plcb4',  
              'Auts2','Kcnj6','Plcb1', 'Grm2',  'Calb2', 'Slc44a5')

# order cells according to factor loadings for "off" factor
mouse_order = ubc_mouse@meta.data[order(ubc_mouse@meta.data$off_factor_loading, decreasing = F),]
mouse_heatmap_ubc = DoHeatmap2(ubc_mouse, use.scaled = T,  genes.use = genes_ubc, 
                               cells.use = row.names(mouse_order), group.by = 'species',
                               slim.col.label = T, group.label.rot = F, draw.line = F,
                               disp.max = 2.5, remove.key = T)

mouse_heatmap_ubc = mouse_heatmap_ubc + theme(axis.title.y=element_blank(),
                                              axis.text.y=element_blank(),
                                              axis.ticks.y=element_blank())

human_order = ubc_hum@meta.data[order(ubc_hum@meta.data$off_factor_loading, decreasing = F),]
human_heatmap_ubc = DoHeatmap2(ubc_hum, use.scaled = T,  genes.use = genes_ubc, 
                               cells.use = row.names(human_order), draw.line = F,
                               slim.col.label = T, group.label.rot = F, group.by = 'species',
                               disp.max = 2.5)

# png('ExtendedDataFig4_panelB.png', units="in", width=9, height=8, res=400)
plot_grid(mouse_heatmap_ubc, human_heatmap_ubc, rel_widths = c(1, 1.9))
# dev.off()

# Panel C

ints_joint = readRDS('../data/mouse_human_joint_mli_pli_obj.RDS')
ints_joint = SetIdent(ints_joint, ident.use = ints_joint@meta.data$joint_cluster)

species_plot = DimPlot(ints_joint, reduction.use = 'tsne', group.by = 'species',
                       do.return = T, no.axes = T, pt.size=0.1) 

ident_plot = DimPlot(ints_joint, reduction.use = 'tsne', group.by = 'joint_cluster',
                       do.return = T, no.axes = T, pt.size=0.1, cols.use = c('#C13EB2', '#3EC14D', 'blue'),
                     no.legend = T, do.label = T) 

# png('ExtendedDataFig4_panelC.png', units="in", width=5, height=7, res=300)
plot_grid(species_plot, ident_plot, ncol = 1)
# dev.off()

# Panel D

ints_mouse = SubsetData(ints_joint, subset.name = 'species', accept.value = 'mouse', 
                       subset.raw = T, do.scale = T, do.center = F)
ints_hum = SubsetData(ints_joint, subset.name = 'species', accept.value = 'human', 
                     subset.raw = T, do.scale = T, do.center = F)

# marker genes for interneuron populations
genes_ints_shared = c('Ptprk','Lhfpl3','Grm8','Lypd6','Lrrn2','Tmeff2','Arl4a','Mb21d2',
                      'Megf11', 'Sorcs3',
                      'Nxph1','Cdh22','Lmo4','Itga8','Utrn','Dlgap3',
                      'Tspan9', 'Klhl1','Chrm2','Gprin3','Grik4','Grm5')

mouse_heatmap = DoHeatmap2(ints_mouse, use.scaled = T,  genes.use = genes_ints_shared, 
                           slim.col.label = T, group.label.rot = T,
                           disp.max = 3.5, title = 'mouse', remove.key = T)
# remove gene labels on right for mouse plot 
# note slight differences may be due to different ordering of cells
mouse_heatmap = mouse_heatmap + theme(axis.title.y=element_blank(),
                                      axis.text.y=element_blank(),
                                      axis.ticks.y=element_blank())

# try to shuffle cells to break up blocks 
cells_human = unlist(lapply(c('MLI1', 'MLI2', 'PLI'), function(x) {
  set.seed(1)
  sample(names(ints_hum@ident)[which(ints_hum@ident == x)])
}), use.names = F)

human_heatmap = DoHeatmap2(ints_hum, use.scaled = T,  genes.use = genes_ints_shared, 
                           cells.use = cells_human,
                           slim.col.label = T, group.label.rot = T,
                           disp.max = 3.5, title = 'human')

# png('ExtendedDataFig4_panelD.png', units="in", width=10, height=9, res=400)
plot_grid(mouse_heatmap, human_heatmap, rel_widths = c(1, 1.4))
# dev.off()

# Panel E

golgi_joint = readRDS('../data/mouse_human_joint_golgi_obj.RDS')

species_plot = DimPlot(golgi_joint, reduction.use = 'tsne', group.by = 'species',
                       do.return = T, no.axes = T, pt.size=0.1) 

ident_plot = DimPlot(golgi_joint, reduction.use = 'tsne', group.by = 'joint_cluster',
                     do.return = T, no.axes = T, pt.size=0.1, cols.use = c('#F8766D', 'blue'),
                     no.legend = T, do.label = T) 

# plotting species and cluster plots
# png('ExtendedDataFig4_panelE.png', units="in", width=5, height=7, res=300)
plot_grid(species_plot, ident_plot, ncol = 1)
# dev.off()

# Panel F
golgi_mouse = SubsetData(golgi_joint, subset.name = 'species', accept.value = 'mouse', 
                        subset.raw = T, do.scale = T, do.center = F)
golgi_hum = SubsetData(golgi_joint, subset.name = 'species', accept.value = 'human', 
                      subset.raw = T, do.scale = T, do.center = F)

genes_golgi = c( 'Vwc2',  'Slc6a5','Grm2', 'Sorcs3', 'Gjd2', 'Sst', 'Nxph1',  'Ebf3', 
                 'Gabrg1')

mouse_heatmap_golgi = DoHeatmap2(golgi_mouse, use.scaled = T,  genes.use = genes_golgi, 
                                 slim.col.label = T, group.label.rot = F, group.spacing = 0.3,
                                 disp.max = 3, title = 'mouse', remove.key = T)
# remove gene labels on right for mouse plot 
mouse_heatmap_golgi = mouse_heatmap_golgi + theme(axis.title.y=element_blank(),
                                                  axis.text.y=element_blank(),
                                                  axis.ticks.y=element_blank())


human_heatmap_golgi = DoHeatmap2(golgi_hum, use.scaled = T,  genes.use = genes_golgi, 
                                 slim.col.label = T, group.label.rot = F, group.spacing = 0.3,
                                 disp.max = 3, title = 'human')

# png('ExtendedDataFig4_panelF.png', width = 9, height = 8, units = 'in', res = 400)
plot_grid(mouse_heatmap_golgi, human_heatmap_golgi, rel_widths = c(1, 1.5))
# dev.off()

# Panel G

granule_joint = readRDS('../data/mouse_human_joint_granule_obj.RDS')
granule_joint = SetIdent(granule_joint, ident.use = granule_joint@meta.data$joint_subcluster)

# Note that DimPlot plots points according to grouping cluster, with smaller cluster on top
# point order can be shuffled if granule_joint object is first converted to liger object and 
# and plotted with liger function (plotFeature)
species_plot = DimPlot(granule_joint, reduction.use = 'umap', group.by = 'species',
                       do.return = T, no.axes = T, pt.size=0.001) 

ident_plot = DimPlot(granule_joint, reduction.use = 'umap', group.by = 'joint_subcluster',
                     do.return = T, no.axes = T, pt.size=0.1,
                     no.legend = T, do.label = T) 

# png('ExtendedDataFig4_panelG.png', units="in", width=5, height=7, res=300)
plot_grid(species_plot, ident_plot, ncol = 1)
# dev.off()

# Panel H
genes_gran = c('Galntl6', 'Kcnq5','Plcl1',  
               'Nr4a2', 'Nr4a3', 'Fos',
               'Chrm3',  'Ldb2', 'Shank2')

# underscores in cluster names causes issues for the plotting function
levels(granule_joint@ident) = c('Granule-1-2', 'Granule-3', 'Granule-4', 'Granule-5')

genes_plot = SplitDotPlot2(granule_joint, grouping.var = 'species', genes.plot =  rev(genes_gran),
                      cols.use = c('red', 'blue'), dot.scale = 9,  x.lab.rot = T, 
                      plot.legend = T, do.return = T)

# first generate region columns
for (lobule in unique(granule_joint@meta.data$region)) {
  print(lobule)
  granule_joint@meta.data[[paste0('Lobule', lobule)]] = sapply(granule_joint@meta.data$region, function(x) {
    as.numeric(x == lobule)
  })
}

lobules = SplitDotPlot2(granule_joint, grouping.var = 'species',
                        features.plot = rev(c('LobuleII', 'LobuleVII', 'LobuleVIII', 'LobuleIX', 'LobuleX')),
                        cols.use = c('orangered', 'cyan'), dot.scale = 9,  x.lab.rot = T,
                        plot.legend = T, do.return = T)

lobules = lobules + theme(axis.text.y = element_blank(),
                          axis.ticks.y = element_blank(),
                          axis.line.y = element_blank())

# pdf('ExtendedDataFig4_panelH.pdf', width = 7, height = 4, useDingbats = F)
plot_grid(genes_plot, lobules, nrow = 1, align = "h", axis = "b", rel_widths = c(1.5, 1))
# dev.off()