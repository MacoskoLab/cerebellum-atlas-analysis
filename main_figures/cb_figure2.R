# Figure 2

# read in full dataset as in panel a -- object full_cb
full_cb = SetIdent(full_cb, ident.use = full_cb@meta.data$subcluster)
purk_s = SubsetData(full_cb, subset.name = 'cluster', accept.value = 'Purkinje', 
                    subset.raw = T)

#  remember to normalize data
purk_s = NormalizeData(purk_s)
genes_purkinje = rev(c("Aldoc",'Kctd16', 'Gpr176', 'Drd3', "Ephb2","Tshz2","Fam19a2",
                       "Alcam","Tox2","2900055J20Rik","Vmn2r30"))
pdf("Figure2_panelB", useDingbats = F, width = 5, height = 3.5)
# use slightly customized dotplot function
DotPlot2(purk_s, genes_purkinje,
         scale.by = 'radius',dot.min = 0.2,dot.scale = 6,
         scale.min = 0.1,plot.legend = T, x.lab.rot = 45)
dev.off()



# panel c

# identifing good Aldoc pos and negative genes (for initial subsplit)
purk_s@meta.data[['Aldoc_stat']] = sapply(as.character(purk_s@ident), function(x) {
  split = strsplit(x, '-')[[1]]
  if (split[1] == 'Anti') {
    return('Aldoc_neg')
  } else {
    return('Aldoc_pos')
  }
})
purk_s = SetIdent(purk_s, ident.use = purk_s@meta.data$Aldoc_stat)

# plot spatial distributions
pdf('Figure2_PanelC_left.pdf', useDingbats = F)
spatial_enrichment_map(purk_s, cluster = 'Aldoc_pos', do.print = F, color_divergent = T)
dev.off()

pdf('Figure2_PanelC_right.pdf', useDingbats = F)
spatial_enrichment_map(purk_s, cluster = 'Aldoc_neg', do.print = F, color_divergent = T)
dev.off()

# panel d 
# map cluster names to numbered cluster names 
# reset the idents and reorder the levels
purk_s = SetIdent(purk_s, ident.use = purk_s@meta.data$subcluster)

# spatial distributions for all purkinje clusters 
for (clust in levels(purk_s@ident)) {
  pdf(paste0('Figure2_panelD_', clust, '.pdf'), useDingbats = F)
  print(spatial_enrichment_map(purk_s, cluster = clust, do.print = F, color_divergent = T))
  dev.off()
}

# panel e 
# remember that full size object (>400k cells)
gran_s = SubsetData(full_cb, subset.name = 'cluster', accept.value = 'Granule', 
                    subset.raw = T)

gran_s = NormalizeData(gran_s)

pdf("Figure2_panelE.pdf",useDingbats = F, width = 5, height = 3.5)
DotPlot2(gran_s,(c("Chrm3","Nr4a2","Galntl6", 'Ebf1', "Gprin3", "Rasgrf1")),dot.scale = 6,
         scale.min = 0.1,plot.legend = T,
         x.lab.rot = 45)
dev.off()


# panel f
# plot spatial distributions of clusters
for (clust_i in 1:3) {
  clust = levels(gran_s@ident)[i]
  pdf(paste0('Figure2_panelF_', clust, '.pdf'), useDingbats = F)
  print(spatial_enrichment_map(gran_s, cluster = clust, do.print = F, color_divergent = T))
  dev.off()
}


# panel g 
# Bergmann glia 
berg_s = SubsetData(full_cb, subset.name = 'cluster', accept.value = 'Bergmann', 
                    subset.raw = T)

pdf('Figure2_panelG', useDingbats = F, 
    height = 3, width = 4)
DotPlot2(berg_s, c('Wif1', 'Mybpc1'), dot.scale = 8,
         scale.min = 0.1,plot.legend = T, scale.by = 'size',
         x.lab.rot = 45)
dev.off()


# panel h
# spatial distribution of cluster
pdf('Figure1_panelH.pdf', useDingbats = F)
spatial_enrichment_map(bergmann_s, cluster = 'Bergmann_2', do.print = F, color_divergent = T)
dev.off()


########################################################################################
# customized functions

library(dplyr)
library(tidyr)
PercentAbove <- function(x, threshold){
  return(length(x = x[x > threshold]) / length(x = x))
}

# modify to allow for label rotation automatically
DotPlot2 = function (object, genes.plot, cols.use = c("lightgrey", "blue"),
                     col.min = -2.5, col.max = 2.5, dot.min = 0, dot.scale = 6,
                     scale.by = "radius", scale.min = NA, scale.max = NA, group.by,
                     plot.legend = FALSE, do.return = FALSE, x.lab.rot = NULL) {
  scale.func <- switch(EXPR = scale.by, size = scale_size,
                       radius = scale_radius, stop("'scale.by' must be either 'size' or 'radius'"))
  if (!missing(x = group.by)) {
    object <- SetAllIdent(object = object, id = group.by)
  }
  data.to.plot <- data.frame(FetchData(object = object, vars.all = genes.plot))
  colnames(x = data.to.plot) <- genes.plot
  data.to.plot$cell <- rownames(x = data.to.plot)
  data.to.plot$id <- object@ident
  #  data.to.plot$dataset <- object@meta.data$orig.ident
  data.to.plot <- data.to.plot %>% gather(key = genes.plot,
                                          value = expression, -c(cell, id))
  data.to.plot <- data.to.plot %>% group_by(id, genes.plot) %>%
    summarize(avg.exp = mean(expm1(x = expression)), pct.exp = PercentAbove(x = expression,
                                                                            threshold = 0))
  data.to.plot <- data.to.plot %>% ungroup() %>% group_by(genes.plot) %>%
    mutate(avg.exp.scale = scale(x = avg.exp)) %>% mutate(avg.exp.scale = MinMax(data = avg.exp.scale,
                                                                                 max = col.max, min = col.min))
  data.to.plot$genes.plot <- factor(x = data.to.plot$genes.plot,
                                    levels = rev(x = genes.plot))
  data.to.plot$pct.exp[data.to.plot$pct.exp < dot.min] <- NA
  p <- ggplot(data = data.to.plot, mapping = aes(x = genes.plot,
                                                 y = id)) + geom_point(shape = 20,mapping = aes(size = pct.exp, color = avg.exp.scale)) +
    scale.func(range = c(0, dot.scale), limits = c(scale.min,
                                                   scale.max)) + theme(axis.title.x = element_blank(),
                                                                       axis.title.y = element_blank())
  if (length(x = cols.use) == 1) {
    p <- p + scale_color_distiller(palette = cols.use)
  }
  else {
    p <- p + scale_color_gradient(low = cols.use[1], high = cols.use[2])
  }
  if (!plot.legend) {
    p <- p + theme(legend.position = "none")
  }
  if (!is.null(x.lab.rot)) {
    p <- p + theme(axis.text.x = element_text(angle = x.lab.rot, hjust = 0.95, vjust=0.95))
  }
  suppressWarnings(print(p))
  if (do.return) {
    return(p)
  }
}

# spatial enrichment map function
# make sure to have coordinates preloaded 
spatial_enrichment_map = function(seurat, gene, cluster = NULL, quantile.p = 0.5, use.pos.expr = T, use.raw = T,
                                  order_regions = NULL, do.print = T, color_divergent = F, return_df = F) {
  # right now this assumes all regions are represented in all cell types -- fix later
  region_prop = table(seurat@meta.data$region) / ncol(seurat@raw.data)
  # watch out for partial matching!
  levels_include = levels(factor(seurat@meta.data$region))
  if (is.null(order_regions)) {
    order_regions = c("I","II", "III", "CUL", "VI", "VII", "VIII",  "IX","X","AN1", "AN2",
                      "PRM", "SIM", "COP", "F",  "PF")
  }
  if (!is.null(cluster)) {
    high_express = names(seurat@ident)[which(seurat@ident == cluster)]
    title = paste0('Cluster ', cluster)
  } else {
    data.use = seurat@raw.data[gene, ]
    if (!use.raw) {
      data.use = seurat@scale.data[gene, ]
      use.scaled = T
    } else {
      use.scaled = F
    }
    if (use.pos.expr) {
      data.use = data.use[data.use > 0]
    }
    cutoff = quantile(data.use, probs = c(quantile.p))
    if (do.print) { print(cutoff) }
    high_express = WhichCells(seurat, subset.name = gene, accept.low = cutoff, 
                              use.raw = use.raw, use.scaled = use.scaled)
    title = paste0(gene, ", quantile = ", quantile.p)
  }
  if (do.print) { print(table(seurat@meta.data[high_express,]$region)) }
  gene_prop = table(factor(seurat@meta.data[high_express,]$region, levels = levels_include)) / 
    sum(table(factor(seurat@meta.data[high_express,]$region, levels = levels_include)))
  if (do.print) {
    barplot((gene_prop / region_prop)[order_regions], las = 2, main = title, 
            ylab = "Relative Proportion of high expressing cells")
    abline(h = 1.0)
  }
  values = data.frame((gene_prop / region_prop)[order_regions])
  row.names(values) = order_regions
  colnames(values) = c('region', 'avg')
  values[['log_avg']] = log(values$avg, base = 2)
  # print(values)
  
  if (return_df) {
    return(values)
  } else {
    # look in cerebellum analysis functions to see how these are made 
    Sps.df<-SpatialPolygonsDataFrame(SPs, values, match.ID = TRUE)
    
    # nn<-data.frame(Sps.df@data)
    # nn<-round(nn, 7)
    
    if (color_divergent) {
      col_names = RColorBrewer:::brewer.pal(11,"RdBu")
      palpos = colorRampPalette(c('white', 'red'), space = 'Lab')(floor((max(values$avg) - 1) * 50))
      palneg = colorRampPalette(c('white', 'blue'), space = 'Lab')(floor((1 - min(values$avg)) * 50))
      
      palette <- c(rev(palneg),palpos)
      
      # final attempt
      cols <- rev(colorRampPalette(brewer.pal(11, "RdBu"))(50))
      max_abs <- max(abs(values$avg - 1))
      # add an extra 0.1 plus because the top color can get lost
      brk <- lattice::do.breaks(c(0.99-max_abs, 1.01 + max_abs), 50)
      # print(brk)
      first_true <- which.max(brk > min(values$avg))
      last_true <- which.max(brk > max(values$avg))
      brk <- brk[(first_true -1):min(last_true, length(brk))]
      cols <- cols[(first_true -1):min(last_true, length(cols))]
      
      spplot(Sps.df, zcol="avg",
             main = title, col.regions= cols,
             at = brk)
      # colorkey = list(col = cols, 
      #                 at = brk))
    } else {
      spplot(Sps.df, zcol="avg",
             main = title, col.regions= viridis(50, option = "C"))
    }
  }
  
  # return((gene_prop / region_prop)[order_regions])
}