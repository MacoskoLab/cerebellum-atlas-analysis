# Additional functions needed for figure 2 and spatial divergence analysis
library(dplyr)
library(tidyr)
library(sp)


# Customized DotPlot function and utilities
# Allows label rotation
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
                                                 y = id)) + geom_point(shape = 20,
                                                                       mapping = aes(size = pct.exp, color = avg.exp.scale)) +
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

# Functions for cerebellar diagram plotting

generate_SPs = function(coords_files, region_names) {
  SPs = lapply(seq_along(coords_files), function(i) {
    x = coords_files[i]
    name = region_names[i]
    pol = Polygon(read.csv(x), hole = F)
    Polygons(list(pol), name)
  })
  SPs = SpatialPolygons(SPs)
  return(SPs)
}

# spatial shapes passed in through sps argument 
# Assumes number of sps matches number of "regions" in metadata object
spatial_enrichment_map_meta = function(seurat_meta, cluster, sps, reg_col = 'region',
                                       meta_col = 'cluster', quantile.p = 0.5, 
                                       use.pos.expr = T, use.raw = T,
                                       order_regions = NULL, do.print = T, 
                                       color_divergent = F, return_df = F) {
  # assumes all regions are represented in all cell types 
  region_prop = table(seurat_meta[[reg_col]]) / nrow(seurat_meta)
  levels_include = levels(factor(seurat_meta[[reg_col]]))
  if (is.null(order_regions)) {
    order_regions = c("I","II", "III", "CUL", "VI", "VII", "VIII",  "IX","X","AN1", "AN2",
                      "PRM", "SIM", "COP", "F",  "PF")
  }
  if (!is.null(cluster)) {
    high_express = row.names(seurat_meta)[which(seurat_meta[[meta_col]] == cluster)]
    title = paste0('Cluster ', cluster)
  } 
  gene_prop = table(factor(seurat_meta[high_express,][[reg_col]], levels = levels_include)) / 
    sum(table(factor(seurat_meta[high_express,][[reg_col]], levels = levels_include)))
  if (do.print) {
    barplot((gene_prop / region_prop)[order_regions], las = 2, main = title, 
            ylab = "Relative Proportion of high expressing cells")
    abline(h = 1.0)
  }
  values = data.frame((gene_prop / region_prop)[order_regions])
  row.names(values) = order_regions
  colnames(values) = c(reg_col, 'avg')
  values[['log_avg']] = log(values$avg, base = 2)
  
  if (return_df) {
    return(values)
  } else {
    # look in cerebellum analysis functions to see how these are made 
    Sps.df<-SpatialPolygonsDataFrame(sps, values, match.ID = TRUE)
    
    if (color_divergent) {
      col_names = RColorBrewer:::brewer.pal(11,"RdBu")

      palpos = function(x){return(colorRamp(col_names[1:6], space = 'Lab')(x) %>% rgb(maxColorValue = 255))}
      palneg = function(x){return(colorRamp(col_names[6:11], space = 'Lab')(x) %>% rgb(maxColorValue = 255))}

      minVal = min(values$avg)
      maxVal = max(values$avg)+0.1
      breaks = lattice::do.breaks(c(minVal, maxVal), 199)

      break_colors = breaks %>% lapply(function(break_val){
        if(break_val == 1){
          return(col_names[6])
        }else if(break_val < 1){
          # palneg from [0=blue, 1=white]
          return(palneg(1-break_val))
        }else{
          #palpos from [1=white, max=red]
          # -> [0=white, max-1=red]
          # -> [0=white, (max-1)/(max-1)=red]
          new_break = (break_val-1)/(maxVal-1)
          if(is.na(new_break)){
            browser()
          }
          if((1 - new_break) < 0) {
            subbreak = 0
          } else {
            subbreak = 1 - new_break
          }
          return(palpos(subbreak))
        }
      }) %>% unlist
      
      spplot(Sps.df, zcol="avg",
             main = title, col.regions= break_colors,
             at = breaks)

    } else {
      spplot(Sps.df, zcol="avg",
             main = title, col.regions= viridis(50, option = "C"))
    }
  }
}


# Tests for significant divergence from expected lobule distribution (chi square)
calcRegionSign_meta = function(seurat_meta, ident, ident_col = 'lob_comp_cluster', division='region',  
                               use_multinomial = F) {
  region_prop = table(seurat_meta[[division]]) / nrow(seurat_meta)
  df = length(region_prop) - 1
  levels_include = levels(factor(seurat_meta[[division]]))
  
  cluster_cells = row.names(seurat_meta)[which(seurat_meta[[ident_col]] == ident)]
  observed = table(factor(seurat_meta[cluster_cells,][[division]], levels = levels_include))
  
  if (use_multinomial) {
    # likelihood ratios based on multinomial test
    observed_ratios = observed / sum(observed)
    test_s = -2 * sum(observed * log(region_prop / observed_ratios))
    
  } else {
    # calc expected values for chi-squared pearsons
    expected = region_prop * table(seurat_meta[[ident_col]])[ident]
    
    test_s = sum((observed - expected)^2 / expected)
  }
  # chi square in both cases
  return(pchisq(test_s, df = df, lower.tail = F))
}

# Computes region divergence significance for all subclusters in a grouping
# Uses BH for multiple test correction
getScatterPlotData_meta = function(seurat_meta, celltype = 'None', ident_col = 'lob_comp_cluster',
                                   division = 'region', use_reps = F) {
  
  # print(levels(seurat_meta[[ident_col]]))
  if (use_reps) {
    split_objs = SplitObject(seurat, attribute.1 = 'rep_desig', subset.raw = T)
    split_dfs = list()
    
    for (rep in c('rep1', 'rep2')) {
      split_dfs[[rep]] = getScatterPlotData_meta(split_objs[[rep]], celltype = celltype, division = division,
                                                 use_reps = F)
      split_dfs[[rep]][['rep']] = rep
    }
    full_df = do.call(rbind, split_dfs)
    return(full_df)
  } else {
    p_values = sapply(levels(seurat_meta[[ident_col]]), function(x) {
      calcRegionSign_meta(seurat_meta, x, division = division)
    })

    names(p_values) = levels(seurat_meta[[ident_col]])
    
    fold_change_values = sapply(levels(seurat_meta[[ident_col]]), function(x) {
      df = spatial_enrichment_map_meta(seurat_meta, cluster = x, meta_col = ident_col,
                                       do.print = F, return_df = T)
      max(df$avg)
    })

    full_df = data.frame(p_values = p_values, fc_values = fold_change_values,
                         celltype = celltype)
    # set pi0 to 1 for Benj Hoch
    full_df[['q_values']] = p.adjust(p_values, method = 'BH')
    row.names(full_df) = levels(seurat_meta[[ident_col]])
    full_df[['subtype']] = row.names(full_df)
    return(full_df)
  }
}