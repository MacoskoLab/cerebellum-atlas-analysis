## Functions shared across cerebellum analyses
library(liger)
library(Seurat)

######################################################################################
## Supplementary liger functions (similar to SeuratExtraFunctions)

### this function plots a series of "binary" tsne plots for a particular feature column in 
# object@cell.data, plotting the presence of one feature, one plot at a time 
# example: using this to generate region enrichment plots for object
# library(cowplot)
# region_plots = plotFeatureSingle(astro_liger, feature = 'region', return.plots = T)
# plot_grid(region_plots[['I']], region_plots[['VI']], region_plots[['VII']], region_plots[['IX']], ncol = 2)
plotFeatureSingle = function(object, feature, by.dataset = F, title = NULL, 
                             subset.keep = NULL, alpha = 0.5,
                             pt.size = 0.3, text.size = 3, do.shuffle = T, rand.seed = 1, do.labels = F,
                             axis.labels = NULL, do.legend = T, legend.size = 5, option = 'plasma', 
                             zero.color = '#F5F5F5', return.plots = F) 
{
  if (is.null(subset.keep)) {
    subset.keep = unique(object@cell.data[[feature]])
  } 
  print(subset.keep)
  plot_list = list()
  for (type in subset.keep) {
    tmp_name = paste0('is_', type)
    object@cell.data[[tmp_name]] = as.character(sapply(object@cell.data[[feature]], function(x){
      x == type
    }))
    plots = plotFeature(object, feature = tmp_name, discrete = T, by.dataset = by.dataset, alpha = alpha,
                        return.plots = T)
    plots = plots + scale_color_manual(values=c(zero.color, '#F8766D'))
    plot_list[[type]] = plots
  }
  return(plot_list)
}

##### This fixes some errors in the original plotFeature function (latest version of liger has this 
# version as plotFeature) -- allows correct merging when by.dataset = F
# this fixes 
plotFeature2 <- function(object, feature, by.dataset = T, discrete = NULL, title = NULL, 
                         pt.size = 0.3, text.size = 3, do.shuffle = T, rand.seed = 1, do.labels = F,
                         alpha = 1, 
                         axis.labels = NULL, do.legend = T, legend.size = 5, option = 'plasma', 
                         zero.color = '#F5F5F5', return.plots = F) {
  dr_df <- data.frame(object@tsne.coords)
  colnames(dr_df) <- c("dr1", "dr2")
  if (!(feature %in% colnames(object@cell.data))) {
    stop('Please select existing feature in cell.data, or add it before calling.')
  }
  dr_df$feature <- object@cell.data[, feature]
  if (is.null(discrete)) {
    if (class(dr_df$feature) != "factor") {
      discrete <- FALSE
    } else {
      discrete <- TRUE
    }
  }
  if (!discrete){
    dr_df$feature[dr_df$feature == 0] <- NA
  }
  if (by.dataset) {
    dr_df$dataset <- object@cell.data$dataset
  } else {
    dr_df$dataset <- factor("single")
  }
  if (do.shuffle) {
    set.seed(rand.seed)
    idx <- sample(1:nrow(dr_df))
    dr_df <- dr_df[idx, ]
  }
  p_list <- list()
  for (sub_df in split(dr_df, f = dr_df$dataset)) {
    ggp <- ggplot(sub_df, aes(x = dr1, y = dr2, color = feature)) + geom_point(size = pt.size, alpha = alpha)
    
    # if data is not discrete
    if (discrete) {
      ggp <- ggp + guides(color = guide_legend(override.aes = list(size = legend.size))) +
        labs(col = feature)
      if (do.labels) {
        centers <- sub_df %>% group_by(feature) %>% summarize(
          dr1 = median(x = dr1),
          dr2 = median(x = dr2)
        )
        ggp <- ggp + geom_text(data = centers, mapping = aes(label = feature),
                               colour = "black", size = text.size)
      }
    } else {
      ggp <- ggp + scale_color_viridis_c(option = option,
                                         direction = -1,
                                         na.value = zero.color) + labs(col = feature)
    }
    
    if (by.dataset) {
      base <- as.character(sub_df$dataset[1])
    } else {
      base <- ""
    }
    if (!is.null(title)) {
      base <- paste(title, base)
    }
    ggp <- ggp + ggtitle(base)
    if (!is.null(axis.labels)) {
      ggp <- ggp + xlab(axis.labels[1]) + ylab(axis.labels[2])
    }
    if (!do.legend) {
      ggp <- ggp + theme(legend.position = "none")
    }
    p_list[[as.character(sub_df$dataset[1])]] <- ggp
  }
  if (by.dataset) {
    p_list <- p_list[names(object@raw.data)]
  }
  
  if (return.plots){
    if (length(p_list) == 1) {
      return(p_list[[1]])
    } else {
      return(p_list)
    }
  } else {
    for (plot in p_list) {
      print(plot)
    }
  }
}

#### an improvement to clusterLouvainJaccard -- allowing users to set more parameters and stores those
# parameters inside parameters slot of liger object (under list name 'Seurat')
clusterLouvainJaccard2 = function (object, resolution = 0.1, k.param = 30, n.iter = 10, 
                                   reduction.type = 'NMF', force.recalc = T, save.SNN = T, 
                                   n.start = 10, print.output = F, dims.use = 1:ncol(object@H.norm), ...) 
{
  temp.seurat = CreateSeuratObject(t(Reduce(rbind, object@scale.data)))
  temp.seurat@scale.data = t(Reduce(rbind, object@scale.data))
  rownames(object@H.norm) = colnames(temp.seurat@scale.data)
  temp.seurat@dr$NMF = new(Class = "dim.reduction", cell.embeddings = object@H.norm, 
                           key = "NMF")
  temp.seurat <- FindClusters(object = temp.seurat, reduction.type = reduction.type, 
                              dims.use = dims.use, force.recalc = force.recalc, save.SNN = save.SNN, 
                              resolution = resolution, k.param = k.param, n.iter = n.iter, 
                              n.start = n.start, print.output = print.output, ...)
  object@clusters = temp.seurat@ident
  if (!('Seurat' %in% names(object@parameters))) {
    object@parameters[['Seurat']] = list()
  }
  object@parameters[['Seurat']][['BuildSNN']] = temp.seurat@calc.params$BuildSNN
  res_params = paste0('FindClusters.res.', resolution)
  object@parameters[['Seurat']][[res_params]] = temp.seurat@calc.params[[res_params]]
  return(object)
}

### Function to add percent_mito column to cell.data slot in liger object 
getPercentMito <- function(object) {
  all.genes = Reduce(union, lapply(object@raw.data, rownames))
  mito.genes <- grep(pattern = "^mt-", x = all.genes, value = TRUE, ignore.case = T)
  percent_mito = unlist(lapply(object@raw.data, function(x) {
    colSums(x[mito.genes, ])/colSums(x)
  }), use.names = F)
  object@cell.data[['percent_mito']] <- percent_mito
  
  return(object)
}

### Function to add CB region column to cell.data slot in liger object 
# (using cell names like IXa_M003_TAGCACAAGCAAATCA)
getCBRegion <- function(object) {
  sup_region = as.character(sapply(rownames(object@cell.data), function(x){
    strsplit(x, '_')[[1]][1]
  }))
  regions = sapply(sup_region, function(x) {
    if (grepl('AN1', x = x)) {
      return('AN1')
    } else if (grepl('AN2', x = x)) {
      return('AN2')
    } else if (grepl('Calb', x) | grepl('Purk', x)) {
      return('Purk')
    } else if (grepl('COP', x = x)) {
      return('COP')
    } else if (grepl('VIII', x)) {
      return('VIII')
    } else if (grepl('VII', x)) {
      return('VII')
    } else if (grepl('VI', x) | grepl('Vl', x) | grepl('DEC', x)) {
      return('VI')
    } else if (grepl('IV', x) | grepl('CUL', x) | grepl('V', x)) {
      return('CUL')
    } else if (grepl('SIM', x)) {
      return('SIM')
    } else if (grepl('PF', x)) {
      return('PF')
    } else if (grepl('F', x)) {
      return('F')
    } else if (grepl('PRM', x)) {
      return('PRM')
    } else if (grepl('IX', x)) {
      return('IX')
    } else if (grepl('X', x)) {
      return('X')
    } else if (grepl('III', x)) {
      return('III')
    } else if (grepl('II', x)) {
      return('II')
    } else if (grepl('I', x)) {
      return('I')
    }
  })
  object@cell.data[['region']] = factor(regions)
  return(object)
}

### Function to add animal column to cell.data slot in liger object 
# (using cell names like IXa_M003_TAGCACAAGCAAATCA)
getCBAnimal <- function(object) {
  animal = as.character(sapply(rownames(object@cell.data), function(x){
    strsplit(x, '_')[[1]][2]
  }))
  object@cell.data[['animal']] = factor(animal)
  return(object)
}

##############################################################################
# plotGene with png
# plotByDatasetAndCluster with png

plotGene_png <- function(object, gene, use.raw = F, use.scaled = F, scale.by = 'dataset', 
                         log2scale = NULL, methylation.indices = NULL, plot.by = 'dataset', 
                         set.dr.lims = F, pt.size = 0.1, min.clip = NULL, max.clip = NULL, 
                         clip.absolute = F, points.only = F, option = 'plasma', cols.use = NULL, 
                         zero.color = '#F5F5F5', axis.labels = NULL, do.legend = T, return.plots = F,
                         use.png = T, width = 7, height = 6, dpi = 100) {
  if ((plot.by != scale.by) & (use.scaled)) {
    warning("Provided values for plot.by and scale.by do not match; results may not be very
            interpretable.")
  }
  if (use.raw) {
    if (is.null(log2scale)) {
      log2scale <- FALSE
    }
    # drop only outer level names
    gene_vals <- getGeneValues(object@raw.data, gene, log2scale = log2scale)
  } else {
    if (is.null(log2scale)) {
      log2scale <- TRUE
    }
    # rescale in case requested gene not highly variable
    if (use.scaled) {
      # check for feature 
      if (!(scale.by %in% colnames(object@cell.data)) & scale.by != 'none') {
        stop("Please select existing feature in cell.data to scale.by, or add it before calling.")
      }
      gene_vals <- getGeneValues(object@norm.data, gene)
      cellnames <- names(gene_vals)
      # set up dataframe with groups
      gene_df <- data.frame(gene = gene_vals)
      if (scale.by == 'none') {
        gene_df[['scaleby']] = 'none'
      } else {
        gene_df[['scaleby']] = factor(object@cell.data[[scale.by]])
      }
      gene_df1 <- gene_df %>%
        group_by(scaleby) %>%
        # scale by selected feature
        mutate_at(vars(-group_cols()), function(x) { scale(x, center = F)})
      gene_vals <- gene_df1$gene
      names(gene_vals) <- cellnames
      if (log2scale) {
        gene_vals <- log2(10000 * gene_vals + 1)
      }
    } else {
      # using normalized data
      # indicate methylation indices here 
      gene_vals <- getGeneValues(object@norm.data, gene, methylation.indices = methylation.indices,
                                 log2scale = log2scale)
    }
  }
  gene_vals[gene_vals == 0] <- NA
  dr_df <- data.frame(object@tsne.coords)
  rownames(dr_df) <- rownames(object@cell.data)
  dr_df$gene <- as.numeric(gene_vals[rownames(dr_df)])
  colnames(dr_df) <- c("dr1", "dr2", "gene")
  # get dr limits for later
  lim1 <- c(min(dr_df$dr1), max(dr_df$dr1))
  lim2 <- c(min(dr_df$dr2), max(dr_df$dr2))
  
  if (plot.by != 'none') {
    if (!(plot.by %in% colnames(object@cell.data))) {
      stop("Please select existing feature in cell.data to plot.by, or add it before calling.")
    }
    dr_df$plotby <- factor(object@cell.data[[plot.by]])
  } else {
    dr_df$plotby <- factor("none")
  }
  # expand clip values if only single provided
  num_levels <- length(levels(dr_df$plotby))
  if (length(min.clip) == 1) {
    min.clip <- rep(min.clip, num_levels)
    names(min.clip) <- levels(dr_df$plotby)
  }
  if (length(max.clip) == 1) {
    max.clip <- rep(max.clip, num_levels)
    names(max.clip) <- levels(dr_df$plotby)
  }
  if (!is.null(min.clip) & is.null(names(min.clip))) {
    if (num_levels > 1) {
      message("Adding names to min.clip according to levels in plot.by group; order may not be 
              preserved as intended if multiple clip values passed in. Pass in named vector to 
              prevent this.")
    }
    names(min.clip) <- levels(dr_df$plotby)
    }
  if (!is.null(max.clip) & is.null(names(max.clip))) {
    if (num_levels > 1) {
      message("Adding names to max.clip according to levels in plot.by group; order may not be 
              preserved as intended if multiple clip values passed in. Pass in named vector to 
              prevent this.")
    }
    names(max.clip) <- levels(dr_df$plotby)
    }
  p_list <- list()
  for (sub_df in split(dr_df, f = dr_df$plotby)) {
    # maybe do quantile cutoff here
    group_name <- as.character(sub_df$plotby[1])
    if (!clip.absolute) {
      max_v <- quantile(sub_df$gene, probs = max.clip[group_name], na.rm = T)
      min_v <- quantile(sub_df$gene, probs = min.clip[group_name], na.rm = T)
    } else {
      max_v <- max.clip[group_name]
      min_v <- min.clip[group_name]
    }
    sub_df$gene[sub_df$gene < min_v & !is.na(sub_df$gene)] <- min_v
    sub_df$gene[sub_df$gene > max_v & !is.na(sub_df$gene)] <- max_v
    
    ggp <- ggplot(sub_df, aes(x = dr1, y = dr2, color = gene)) + geom_point(size = pt.size) +
      labs(col = gene)
    
    if (!is.null(cols.use)) {
      ggp <- ggp + scale_color_gradientn(colors = cols.use,
                                         na.value = zero.color)
    } else {
      ggp <- ggp + scale_color_viridis_c(option = option,
                                         direction = -1,
                                         na.value = zero.color)
    }
    if (set.dr.lims) {
      ggp <- ggp + xlim(lim1) + ylim(lim2)
    }
    
    if (plot.by != 'none') {
      base <- as.character(sub_df$plotby[1])
    } else {
      base <- ""
    }
    ggp <- ggp + ggtitle(base)
    
    if (!is.null(axis.labels)) {
      ggp <- ggp + xlab(axis.labels[1]) + ylab(axis.labels[2])
    }
    if (!do.legend) {
      ggp <- ggp + theme(legend.position = "none")
    }
    if (points.only) {
      ggp <- ggp + theme(
        axis.line = element_blank(), axis.text.x = element_blank(),
        axis.text.y = element_blank(), axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(), legend.position = "none",
        panel.background = element_blank(), panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.background = element_blank(), plot.title = element_blank()
      )
    }
    if (use.png) {
      legend_s = get_legend(ggp)
      ggp <- AugmentPlot(ggp, width = width, height = height, dpi = dpi) + legend_s
    }
    p_list[[as.character(sub_df$plotby[1])]] <- ggp
  }
  if (plot.by == 'dataset') {
    p_list <- p_list[names(object@raw.data)]
  }
  
  if (return.plots){
    if (length(p_list) == 1) {
      return(p_list[[1]])
    } else {
      return(p_list)
    }
  } else {
    for (plot in p_list) {
      print(plot)
    }
  }
  }

plotByDatasetAndCluster_png <- function(object, clusters = NULL, title = NULL, pt.size = 0.3,
                                        text.size = 3, do.shuffle = T, rand.seed = 1,
                                        axis.labels = NULL, do.legend = T, legend.size = 5,
                                        return.plots = F, use.png = T, width = 7, height = 6, dpi = 100) {
  tsne_df <- data.frame(object@tsne.coords)
  colnames(tsne_df) <- c("tsne1", "tsne2")
  tsne_df$Dataset <- unlist(lapply(1:length(object@H), function(x) {
    rep(names(object@H)[x], nrow(object@H[[x]]))
  }))
  c_names <- names(object@clusters)
  if (is.null(clusters)) {
    # if clusters have not been set yet
    if (length(object@clusters) == 0) {
      clusters <- rep(1, nrow(object@tsne.coords))
      names(clusters) <- c_names <- rownames(object@tsne.coords)
    } else {
      clusters <- object@clusters
      c_names <- names(object@clusters)
    }
  }
  tsne_df$Cluster <- clusters[c_names]
  if (do.shuffle) {
    set.seed(rand.seed)
    idx <- sample(1:nrow(tsne_df))
    tsne_df <- tsne_df[idx, ]
  }
  
  p1 <- ggplot(tsne_df, aes(x = tsne1, y = tsne2, color = Dataset)) +
    geom_point(size = pt.size) +
    guides(color = guide_legend(override.aes = list(size = legend.size)))
  
  centers <- tsne_df %>% group_by(Cluster) %>% summarize(
    tsne1 = median(x = tsne1),
    tsne2 = median(x = tsne2)
  )
  p2 <- ggplot(tsne_df, aes(x = tsne1, y = tsne2, color = Cluster)) + geom_point(size = pt.size) +
    guides(color = guide_legend(override.aes = list(size = legend.size)))
  
  if (!use.png) {
    p2 <- p2 + geom_text(data = centers, mapping = aes(label = Cluster), 
                         colour = "black", size = text.size) 
  }
  
  if (!is.null(title)) {
    p1 <- p1 + ggtitle(title[1])
    p2 <- p2 + ggtitle(title[2])
  }
  if (!is.null(axis.labels)) {
    p1 <- p1 + xlab(axis.labels[1]) + ylab(axis.labels[2])
    p2 <- p2 + xlab(axis.labels[1]) + ylab(axis.labels[2])
  }
  if (!do.legend) {
    p1 <- p1 + theme(legend.position = "none")
    p2 <- p2 + theme(legend.position = "none")
  }
  if (use.png) {
    p1_legend = get_legend(p1)
    p2_legend = get_legend(p2)
    p1 <- AugmentPlot(p1, width = width, height = height, dpi = dpi) + p1_legend
    p2 <- AugmentPlot(p2, width = width, height = height, dpi = dpi) + 
      geom_text(data = centers, mapping = aes(tsne1, tsne2, label = Cluster), 
                colour = "black", size = text.size) + p2_legend
  }
  if (return.plots) {
    return(list(p1, p2))
  } else {
    print(p1)
    print(p2)
  }
}

###############################################################################
## Next few functions all have to do with plotFactor_new (with png option)
get_legend <- function(p) {
  tmp <- ggplot_gtable(ggplot_build(p))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}


# note that setting point size to -1 seems to make it disappear 
# (trick that I can't find in documentation, but got from Seurat)
# also tried Pethukov's ggrastr, but it caused too much distortion 
# also started taking way too long when playing with raster.width > 4

# still working on fixing legend for first page of plots 

plotFactors_new <- function(object, num.genes = 8, cells.highlight = NULL, plot.tsne = F,
                            return.plots = F, pt.size1 = 0.5, pt.size2 = 1, 
                            tsne.colors = c('lemonchiffon', 'red'), use.png = F,
                            width = 7, height = 6, dpi = 100) {
  k <- ncol(object@H.norm)
  pb <- txtProgressBar(min = 0, max = k, style = 3)
  
  W <- t(object@W)
  rownames(W) <- colnames(object@scale.data[[1]])
  Hs_norm <- object@H.norm
  H_raw = do.call(rbind, object@H)
  plot_list = list()
  tsne_list = list()
  
  for (i in 1:k) {
    # creates array of 2 by 1 plots 
    par(mfrow = c(2, 1))
    top_genes.W <- rownames(W)[order(W[, i], decreasing = T)[1:num.genes]]
    top_genes.W.string <- paste0(top_genes.W, collapse = ", ")
    factor_textstring <- paste0("Factor", i)
    
    plot_title1 <- paste(factor_textstring, "\n", top_genes.W.string, "\n")
    
    h_df = data.frame(x = 1:nrow(Hs_norm), h_norm = Hs_norm[, i],
                      h_raw = H_raw[, i], dataset = object@cell.data$dataset,
                      highlight = FALSE)
    
    top <- ggplot(h_df, aes(x = x, y=h_raw, col = dataset)) + geom_point(size = pt.size1) + 
      labs(x = 'Cell', y = 'Raw H Score') + ggtitle(plot_title1) + theme(legend.position = 'none')
    bottom <- ggplot(h_df, aes(x = x, y=h_norm, col = dataset)) + geom_point(size = pt.size1) + 
      labs(x = 'Cell', y = 'H_norm Score') + 
      theme(legend.position = 'top',
            legend.title = element_blank()) +
      guides(colour = guide_legend(override.aes = list(size = 2)))
    
    if (!is.null(cells.highlight)) {
      h_df[cells.highlight, 'highlight'] = TRUE
      top <- top + geom_point(data = subset(h_df, highlight == TRUE),
                              aes(x, h_raw),
                              col = "black",
                              size = pt.size1)
      bottom <- bottom + geom_point(data = subset(h_df, highlight == TRUE),
                                    aes(x, h_norm),
                                    col = "black",
                                    size = pt.size1)
    }
    
    if (use.png){
      top = AugmentPlot(top, width = width, height = height, dpi = dpi)
      top = top + labs(x = 'Cell', y = 'Raw H Score')
      # save legend before losing it -- might have to recreate manually
      bottom_legend = get_legend(bottom) 
      bottom = AugmentPlot(bottom, width = width, height = height, dpi = dpi)
      bottom = bottom + labs(x = 'Cell', y = 'H_norm Score') 
    }
    
    full = plot_grid(top, bottom, ncol = 1, rel_heights = c(1, 1))
    plot_list[[i]] = full
    if (plot.tsne) {
      tsne_df <- data.frame(Hs_norm[, i], object@tsne.coords)
      factorlab <- paste0("Factor", i)
      colnames(tsne_df) <- c(factorlab, "tSNE1", "tSNE2")
      p1 <- ggplot(tsne_df, aes_string(x = "tSNE1", y = "tSNE2", color = factorlab)) + 
        geom_point(size = pt.size2) +
        scale_color_gradient(low = tsne.colors[1], high = tsne.colors[2]) + ggtitle(label = paste('Factor', i)) + 
        theme(legend.position = 'none')
      if (use.png){
        p1 = AugmentPlot(p1, width = width, height = height, dpi = dpi)
      }
      tsne_list[[i]] = p1
    }
    setTxtProgressBar(pb, i)
  }
  if (return.plots) {
    return(list(plot_list, tsne_list))
  } else {
    for (i in 1:k) {
      print(plot_list[[i]])
      print(tsne_list[[i]])
    }
  }
}

library(png)
# AugmentPlot taken from Seurat v3

#' Augments ggplot2-based plot with a PNG image.
#'
#' Creates "vector-friendly" plots. Does this by saving a copy of the plot as a PNG file,
#' then adding the PNG image with \code{\link[ggplot2]{annotation_raster}} to a blank plot
#' of the same dimensions as \code{plot}. Please note: original legends and axes will be lost
#' during augmentation.
#'
#' @param plot A ggplot object
#' @param width,height Width and height of PNG version of plot
#' @param dpi Plot resolution
#'
#' @return A ggplot object
#'
#' @importFrom png readPNG
#' @importFrom ggplot2 ggplot_build ggsave ggplot aes_string geom_blank annotation_raster ggtitle
#'
#' @export
#'
#' @examples
#' \dontrun{
#' plot <- DimPlot(object = pbmc_small)
#' AugmentPlot(plot = plot)
#' }
#'
AugmentPlot <- function(plot, width = 10, height = 10, dpi = 100) {
  pbuild.params <- ggplot_build(plot = plot)$layout$panel_params[[1]]
  range.values <- c(
    pbuild.params$x.range,
    pbuild.params$y.range
  )
  xyparams <- GetXYAesthetics(
    plot = plot,
    geom = class(x = plot$layers[[1]]$geom)[1]
  )
  title <- plot$labels$title
  tmpfile <- tempfile(fileext = '.png')
  ggsave(
    filename = tmpfile,
    plot = plot + NoLegend() + NoAxes() + theme(plot.title = element_blank()),
    width = width,
    height = height,
    dpi = dpi
  )
  img <- readPNG(source = tmpfile)
  file.remove(tmpfile)
  blank <- ggplot(
    data = plot$data,
    mapping = aes_string(x = xyparams$x, y = xyparams$y)
  ) + geom_blank()
  blank <- blank + plot$theme + ggtitle(label = title)
  blank <- blank + annotation_raster(
    raster = img,
    xmin = range.values[1],
    xmax = range.values[2],
    ymin = range.values[3],
    ymax = range.values[4]
  )
  return(blank)
}

# other functions which are necessary -- GetXYAesthetics, NoLegend, NoAxes

# Get X and Y aesthetics from a plot for a certain geom
#
# @param plot A ggplot2 object
# @param geom Geom class to filter to
# @param plot.first Use plot-wide X/Y aesthetics before geom-specific aesthetics
#
# @return A named list with values 'x' for the name of the x aesthetic and 'y' for the y aesthetic
#
GetXYAesthetics <- function(plot, geom = 'GeomPoint', plot.first = TRUE) {
  geoms <- sapply(
    X = plot$layers,
    FUN = function(layer) {
      return(class(x = layer$geom)[1])
    }
  )
  geoms <- which(x = geoms == geom)
  if (length(x = geoms) == 0) {
    stop("Cannot find a geom of class ", geom)
  }
  geoms <- min(geoms)
  if (plot.first) {
    x <- as.character(x = plot$mapping$x %||% plot$layers[[geoms]]$mapping$x)[2]
    y <- as.character(x = plot$mapping$y %||% plot$layers[[geoms]]$mapping$y)[2]
  } else {
    x <- as.character(x = plot$layers[[geoms]]$mapping$x %||% plot$mapping$x)[2]
    y <- as.character(x = plot$layers[[geoms]]$mapping$y %||% plot$mapping$y)[2]
  }
  return(list('x' = x, 'y' = y))
}

NoLegend <- function(...) {
  no.legend.theme <- theme(
    # Remove the legend
    legend.position = 'none',
    # Validate the theme
    validate = TRUE,
    ...
  )
  return(no.legend.theme)
}

#' @inheritParams SeuratTheme
#' @param keep.text Keep axis text
#' @param keep.ticks Keep axis ticks
#'
#' @importFrom ggplot2 theme element_blank
#' @export
#'
#' @rdname SeuratTheme
#' @aliases NoAxes
#'
#' @examples
#' # Generate a plot with no axes
#' library(ggplot2)
#' df <- data.frame(x = rnorm(n = 100, mean = 20, sd = 2), y = rbinom(n = 100, size = 100, prob = 0.2))
#' p <- ggplot(data = df, mapping = aes(x = x, y = y)) + geom_point(mapping = aes(color = 'red'))
#' p + NoAxes()
#'
NoAxes <- function(..., keep.text = FALSE, keep.ticks = FALSE) {
  blank <- element_blank()
  no.axes.theme <- theme(
    # Remove the axis elements
    axis.line.x = blank,
    axis.line.y = blank,
    # Validate the theme
    validate = TRUE,
    ...
  )
  if (!keep.text) {
    no.axes.theme <- no.axes.theme + theme(
      axis.text.x = blank,
      axis.text.y = blank,
      axis.title.x = blank,
      axis.title.y = blank,
      validate = TRUE,
      ...
    )
  }
  if (!keep.ticks){
    no.axes.theme <- no.axes.theme + theme(
      axis.ticks.x = blank,
      axis.ticks.y = blank,
      validate = TRUE,
      ...
    )
  }
  return(no.axes.theme)
}

# also some hadley wickham functions 
`%||%` <- function(lhs, rhs) {
  if (!is.null(x = lhs)) {
    return(lhs)
  } else {
    return(rhs)
  }
}




###############################################################################
## Supplementary Seurat functions (original SeuratExtraFunctions further below)

# This function takes all possible pairs of clusters in a Seurat object and finds
# DE genes between them -- concatenates all dataframes into single result and returns
# it also downsamples by default to speed it up 
# ... should be additional parameters for FindMarkers
FindMarkersPairwise = function(object, max.cells.per.ident = 700, ...) {
  clusters <- levels(object@ident)
  pairwise <- combn(clusters, 2)
  
  results <- list()
  for(i in 1:ncol(pairwise)) {
    message('Finding DE genes between ', pairwise[1, i], ' and ', pairwise[2, i])
    markers_pair <- FindMarkers(object, ident.1 = pairwise[1, i], ident.2 = pairwise[2, i], 
                                max.cells.per.ident = max.cells.per.ident, ...)
    comparisons <- pairwise[, i]
    markers_pair$comparison <- paste(comparisons[1], comparisons[2], sep = '_')
    markers_pair$gene <- rownames(markers_pair)
    results[[i]] <- markers_pair
  }
  
  results.df <- do.call(rbind, results)
  return(results.df)
}

##############################################################################
# feature plot with raw data instead of scaled data
# also has necessary internal function dependencies below (all functions until next section
# are just dependencies)

FeaturePlot_raw = function (object, features.plot, min.cutoff = NA, max.cutoff = NA, 
                            dim.1 = 1, dim.2 = 2, cells.use = NULL, pt.size = 1, cols.use = c("yellow", 
                                                                                              "red"), 
                            pch.use = 16, overlay = FALSE, do.hover = FALSE, 
                            data.hover = "ident", do.identify = FALSE, reduction.use = "tsne", 
                            use.imputed = FALSE, nCol = NULL, no.axes = FALSE, no.legend = TRUE, 
                            coord.fixed = FALSE, dark.theme = FALSE, do.return = FALSE, 
                            vector.friendly = FALSE, png.file = NULL, png.arguments = c(10, 
                                                                                        10, 100)) 
{
  cells.use <- SetIfNull(x = cells.use, default = colnames(x = object@data))
  if (is.null(x = nCol)) {
    nCol <- 2
    if (length(x = features.plot) == 1) {
      nCol <- 1
    }
    if (length(x = features.plot) > 6) {
      nCol <- 3
    }
    if (length(x = features.plot) > 9) {
      nCol <- 4
    }
  }
  num.row <- floor(x = length(x = features.plot)/nCol - 1e-05) + 
    1
  if (overlay | do.hover) {
    num.row <- 1
    nCol <- 1
  }
  par(mfrow = c(num.row, nCol))
  dim.code <- GetDimReduction(object = object, reduction.type = reduction.use, 
                              slot = "key")
  dim.codes <- paste0(dim.code, c(dim.1, dim.2))
  data.plot <- as.data.frame(GetCellEmbeddings(object = object, 
                                               reduction.type = reduction.use, dims.use = c(dim.1, dim.2), 
                                               cells.use = cells.use))
  x1 <- paste0(dim.code, dim.1)
  x2 <- paste0(dim.code, dim.2)
  data.plot$x <- data.plot[, x1]
  data.plot$y <- data.plot[, x2]
  data.plot$pt.size <- pt.size
  names(x = data.plot) <- c("x", "y")
  data.use <- t(x = FetchData(object = object, vars.all = features.plot, 
                              cells.use = cells.use, use.imputed = use.imputed, use.raw = T))
  min.cutoff <- mapply(FUN = function(cutoff, feature) {
    ifelse(test = is.na(x = cutoff), yes = min(data.use[feature, 
                                                        ]), no = cutoff)
  }, cutoff = min.cutoff, feature = features.plot)
  max.cutoff <- mapply(FUN = function(cutoff, feature) {
    ifelse(test = is.na(x = cutoff), yes = max(data.use[feature, 
                                                        ]), no = cutoff)
  }, cutoff = max.cutoff, feature = features.plot)
  check_lengths = unique(x = vapply(X = list(features.plot, 
                                             min.cutoff, max.cutoff), FUN = length, FUN.VALUE = numeric(length = 1)))
  if (length(x = check_lengths) != 1) {
    stop("There must be the same number of minimum and maximum cuttoffs as there are features")
  }
  if (overlay) {
    pList <- list(BlendPlot(data.use = data.use, features.plot = features.plot, 
                            data.plot = data.plot, pt.size = pt.size, pch.use = pch.use, 
                            cols.use = cols.use, dim.codes = dim.codes, min.cutoff = min.cutoff, 
                            max.cutoff = max.cutoff, coord.fixed = coord.fixed, 
                            no.axes = no.axes, no.legend = no.legend, dark.theme = dark.theme))
  }
  else {
    pList <- mapply(FUN = SingleFeaturePlot, feature = features.plot, 
                    min.cutoff = min.cutoff, max.cutoff = max.cutoff, 
                    coord.fixed = coord.fixed, MoreArgs = list(data.use = data.use, 
                                                               data.plot = data.plot, pt.size = pt.size, pch.use = pch.use, 
                                                               cols.use = cols.use, dim.codes = dim.codes, no.axes = no.axes, 
                                                               no.legend = no.legend, dark.theme = dark.theme, 
                                                               vector.friendly = vector.friendly, png.file = png.file, 
                                                               png.arguments = png.arguments), SIMPLIFY = FALSE)
  }
  if (do.hover) {
    if (length(x = pList) != 1) {
      stop("'do.hover' only works on a single feature or an overlayed FeaturePlot")
    }
    if (is.null(x = data.hover)) {
      features.info <- NULL
    }
    else {
      features.info <- FetchData(object = object, vars.all = data.hover)
    }
    return(HoverLocator(plot = pList[[1]], data.plot = data.plot, 
                        features.info = features.info, dark.theme = dark.theme, 
                        title = features.plot))
  }
  else if (do.identify) {
    if (length(x = pList) != 1) {
      stop("'do.identify' only works on a single feature or an overlayed FeaturePlot")
    }
    return(FeatureLocator(plot = pList[[1]], data.plot = data.plot, 
                          dark.theme = dark.theme))
  }
  else {
    print(x = cowplot::plot_grid(plotlist = pList, ncol = nCol))
  }
  ResetPar()
  if (do.return) {
    return(pList)
  }
}


SetIfNull <- function(x, default) {
  if (is.null(x = x)) {
    return(default)
  } else {
    return(x)
  }
}

# globalVariables(names = c('x', 'y', 'gene'), package = 'Seurat', add = TRUE)
SingleFeaturePlot <- function(
  data.use,
  feature,
  data.plot,
  pt.size,
  pch.use,
  cols.use,
  dim.codes,
  min.cutoff,
  max.cutoff,
  coord.fixed,
  no.axes,
  no.title = FALSE,
  no.legend,
  dark.theme,
  vector.friendly = FALSE,
  png.file = NULL,
  png.arguments=c(10,10,100)
) {
  #first, consider vector friendly case
  if (vector.friendly) {
    previous_call <- blank_call <- png_call <- match.call()
    blank_call$pt.size <- -1
    blank_call$vector.friendly <- FALSE
    png_call$no.axes <- TRUE
    png_call$no.legend <- TRUE
    png_call$vector.friendly <- FALSE
    png_call$no.title <- TRUE
    blank_plot <- eval(blank_call, sys.frame(sys.parent()))
    png_plot <- eval(png_call, sys.frame(sys.parent()))
    png.file <- SetIfNull(x = png.file, default = paste0(tempfile(), ".png"))
    ggsave(filename = png.file, plot = png_plot,
           width = png.arguments[1],
           height = png.arguments[2],
           dpi = png.arguments[3])
    to_return <- AugmentPlot(blank_plot, png.file)
    file.remove(png.file)
    return(to_return)
  }
  data.gene <- na.omit(object = data.frame(data.use[feature, ]))
  #   Check for quantiles
  min.cutoff <- SetQuantile(cutoff = min.cutoff, data = data.gene)
  max.cutoff <- SetQuantile(cutoff = max.cutoff, data = data.gene)
  #   Mask any values below the minimum and above the maximum values
  data.gene <- sapply(
    X = data.gene,
    FUN = function(x) {
      return(ifelse(test = x < min.cutoff, yes = min.cutoff, no = x))
    }
  )
  data.gene <- sapply(
    X = data.gene,
    FUN = function(x) {
      return(ifelse(test = x > max.cutoff, yes = max.cutoff, no = x))
    }
  )
  data.plot$gene <- data.gene
  #   Stuff for break points
  if (length(x = cols.use) == 1) {
    brewer.gran <- brewer.pal.info[cols.use, ]$maxcolors
  } else {
    brewer.gran <- length(x = cols.use)
  }
  #   Cut points
  if (all(data.gene == 0)) {
    data.cut <- 0
  } else {
    data.cut <- as.numeric(x = as.factor(x = cut(
      x = as.numeric(x = data.gene),
      breaks = brewer.gran
    )))
  }
  data.plot$col <- as.factor(x = data.cut)
  #   Start plotting
  p <- ggplot(data = data.plot, mapping = aes(x = x, y = y))
  if (brewer.gran != 2) {
    if (length(x = cols.use) == 1) {
      p <- p + geom_point(
        mapping = aes(color = col),
        size = pt.size,
        shape = pch.use
      ) + scale_color_brewer(palette = cols.use)
    } else {
      p <- p + geom_point(
        mapping = aes(color = col),
        size = pt.size,
        shape = pch.use
      ) + scale_color_manual(values = cols.use)
    }
  } else {
    if (all(data.plot$gene == data.plot$gene[1])) {
      warning(paste0("All cells have the same value of ", feature, "."))
      p <- p + geom_point(color = cols.use[1], size = pt.size, shape = pch.use)
    } else {
      p <- p + geom_point(
        mapping = aes(color = gene),
        size = pt.size,
        shape = pch.use
      ) + scale_color_gradientn(
        colors = cols.use,
        guide = guide_colorbar(title = feature)
      )
    }
  }
  if (dark.theme) {
    p <- p + DarkTheme()
  }
  if (no.axes) {
    p <- p + theme(
      axis.line = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank()
    )
    if (!no.title) p <- p + labs(title = feature, x ="", y="")
    if (no.title) p <- p + labs(x ="", y="")
  } else {
    if (no.title) p <- p + labs(x = dim.codes[1], y = dim.codes[2])
    if (!(no.title)) p <- p + labs(title = feature, x = dim.codes[1], y = dim.codes[2])
  }
  if (no.legend) {
    p <- p + theme(legend.position = 'none')
  }
  if(coord.fixed) {
    p <- p + coord_fixed()
  }
  
  return(p)
}

library(stats)

# Find the quantile of a data
#
# Converts a quantile in character form to a number regarding some data
# String form for a quantile is represented as a number prefixed with 'q'
# For example, 10th quantile is 'q10' while 2nd quantile is 'q2'
#
# Will only take a quantile of non-zero data values
#
SetQuantile <- function(cutoff, data) {
  if (grepl(pattern = '^q[0-9]{1,2}$', x = as.character(x = cutoff), perl = TRUE)) {
    this.quantile <- as.numeric(x = sub(
      pattern = 'q',
      replacement = '',
      x = as.character(x = cutoff)
    )) / 100
    data <- unlist(x = data)
    data <- data[data > 0]
    cutoff <- quantile(x = data, probs = this.quantile)
  }
  return(as.numeric(x = cutoff))
}

# Reset Par
#
# Reset the graphing space to
# mfrow = c(1, 1)
#
ResetPar <- function(...) {
  par(mfrow = c(1, 1), ...)
}


# spatial gene selection 
###
# give option to do the averaging in raw gene space or scaled gene space 
# both are by default false
# cell.var is the column in cell.data to split by 
# see required utility function, ligerToSeurat_raw below 
library(pastecs)
spatial_variable_select_clean <-function(liger_object, cell.var = "regions", use.raw = F, 
                                         use.scale = F, alpha = 0.01){
  # convert to Seurat for averaging
  tmp.seurat = ligerToSeurat_raw(liger_object, renormalize = T, raw.only = T)
  tmp.seurat@meta.data[[cell.var]] = liger_object@cell.data[[cell.var]]
  tmp.seurat = SetIdent(tmp.seurat, ident.use = tmp.seurat@meta.data[[cell.var]])
  # get average expression for each region (returns dataframe)
  avg_lobe_expression_mat = AverageExpression(tmp.seurat, return.seurat = F,
                                              use.raw = use.raw, use.scale = use.scale)
  
  #Between lobe variability 
  mean<-rowMeans(avg_lobe_expression_mat)
  var<-RowVar(avg_lobe_expression_mat)
  #getting VMR
  avg_lobe_expression_mat[['VMR']] = var/mean
  avg_lobe_expression_mat[['logVMR']] = log(avg_lobe_expression_mat$VMR)
  
  # getting spatially significant gene list
  
  vmrs =as.matrix(avg_lobe_expression_mat[,"logVMR"])
  # convert to density (empirical dist)
  d = density(na.omit(vmrs))
  ts_y = ts(d$y)
  # inflection points of distribution
  tp = turnpoints(ts_y)
  idx = which.max(d$y[tp$tppos])
  # "mean" of distribution? 
  center=d$x[tp$tppos[idx]]
  belowcenter = vmrs[which(vmrs < center)]
  # building a null distribution?
  null = c(belowcenter,(belowcenter + 2*(center-belowcenter)))
  pvalsUpper=pnorm(vmrs,mean=center,sd= sd(null),lower.tail=F)
  # adjust p-values with Benjamini Hochberg correction 
  pvalsUpper=p.adjust(pvalsUpper,"fdr")
  sign_idx=which(pvalsUpper <= alpha)
  
  # get minimum significant vmr
  # cutoff_vmr = vmrs[sign_idx[which.min(vmrs[sign_idx])]]
  cutoff_vmr = min(vmrs[sign_idx])
  
  # plotting significant cutoff for vmr
  hist(vmrs)
  abline(v=cutoff_vmr)
  
  genelist<-rownames(avg_lobe_expression_mat)[which(avg_lobe_expression_mat$logVMR > cutoff_vmr)]
  #removing genes that arent in norm.data
  liger_gene_list = lapply(liger_object@norm.data, FUN = rownames)
  liger_gene_list[[length(liger_gene_list) + 1]] = genelist
  final_list = Reduce(intersect, liger_gene_list)
  return(final_list)
}

# define row variance function 
#find the row variance
# taken from stack overflow 
# https://stackoverflow.com/questions/25099825/row-wise-variance-of-a-matrix-in-r
RowVar <- function(x, ...) {
  rowSums((x - rowMeans(x, ...))^2, ...)/(dim(x)[2] - 1)
}

# keeping only basic values 
# this version of the function is required for spatial variable selection
ligerToSeurat_raw = function (object, nms = names(object@H), renormalize = T, use.liger.genes = T, 
                              raw.only = F, by.dataset = F) 
{
  if (!require("Seurat", quietly = TRUE)) {
    stop("Package \"Seurat\" needed for this function to work. Please install it.", 
         call. = FALSE)
  }
  maj_version <- packageVersion("Seurat")$major
  if (class(object@raw.data[[1]])[1] != "dgCMatrix") {
    mat <- as(x, "CsparseMatrix")
    object@raw.data <- lapply(object@raw.data, function(x) {
      as(x, "CsparseMatrix")
    })
  }
  raw.data <- MergeSparseDataAll(object@raw.data, nms)
  if (!raw.only) {
    scale.data <- do.call(rbind, object@scale.data)
    rownames(scale.data) <- colnames(raw.data)
    if (maj_version < 3) {
      var.genes <- object@var.genes
      inmf.obj <- new(Class = "dim.reduction", gene.loadings = t(object@W), 
                      cell.embeddings = object@H.norm, key = "iNMF_")
      rownames(inmf.obj@gene.loadings) <- var.genes
      tsne.obj <- new(Class = "dim.reduction", cell.embeddings = object@tsne.coords, 
                      key = "tSNE_")
    }
    else {
      var.genes <- object@var.genes
      if (any(grepl("_", var.genes))) {
        print("Warning: Seurat v3 genes cannot have underscores, replacing with dashes ('-')")
        var.genes <- gsub("_", replacement = "-", var.genes)
      }
      inmf.obj <- new(Class = "DimReduc", feature.loadings = t(object@W), 
                      cell.embeddings = object@H.norm, key = "iNMF_")
      rownames(inmf.obj@feature.loadings) <- var.genes
      tsne.obj <- new(Class = "DimReduc", cell.embeddings = object@tsne.coords, 
                      key = "tSNE_")
    }
    rownames(tsne.obj@cell.embeddings) <- rownames(scale.data)
    rownames(inmf.obj@cell.embeddings) <- rownames(scale.data)
    colnames(tsne.obj@cell.embeddings) <- paste0("tSNE_", 1:2)
  }
  new.seurat <- Seurat::CreateSeuratObject(raw.data)
  if (renormalize) {
    new.seurat <- Seurat::NormalizeData(new.seurat)
  }
  if (raw.only){
    return(new.seurat)
  }
  if (by.dataset) {
    ident.use <- as.character(unlist(lapply(1:length(object@raw.data), 
                                            function(i) {
                                              dataset.name <- names(object@raw.data)[i]
                                              paste0(dataset.name, as.character(object@clusters[colnames(object@raw.data[[i]])]))
                                            })))
  }
  else {
    ident.use <- as.character(object@clusters)
  }
  if (maj_version < 3) {
    if (use.liger.genes) {
      new.seurat@var.genes <- var.genes
    }
    new.seurat@scale.data <- t(scale.data)
    new.seurat@dr$tsne <- tsne.obj
    new.seurat@dr$inmf <- inmf.obj
    new.seurat <- SetIdent(new.seurat, ident.use = ident.use)
  }
  else {
    if (use.liger.genes) {
      new.seurat@assays$RNA@var.features <- var.genes
    }
    Seurat::SetAssayData(new.seurat, slot = "scale.data", 
                         t(scale.data), assay = "RNA")
    new.seurat@reductions$tsne <- tsne.obj
    new.seurat@reductions$inmf <- inmf.obj
    Idents(new.seurat) <- ident.use
  }
  return(new.seurat)
}

library(Matrix.utils)
spatial_variable_select_clean_fast<-function(liger_object, cell.var = "regions", use.raw = F, 
                                             use.scale = F, alpha = 0.01,density.bandwidth = "nrd0"){
  # convert to Seurat for averaging
  groupings = factor(liger_object@cell.data[[cell.var]])
  names(groupings) = rownames(liger_object@cell.data)
  # dge = MergeSparseDataAll(liger_object@raw.data)
  dge.norm = MergeSparseDataAll(liger_object@norm.data)
  x = Sparse_transpose(dge.norm)
  sums=aggregate.Matrix(x=x,groupings = groupings[intersect(names(groupings),colnames(dge.norm))])
  sums = as.matrix(sums)
  tab = table(groupings)
  avg_lobe_expression_mat = as.data.frame(t(sweep(sums,1,tab,"/")))
  rownames(avg_lobe_expression_mat) = rownames(dge.norm)
  #Between lobe variability 
  mean<-rowMeans(avg_lobe_expression_mat)
  var<-RowVar(avg_lobe_expression_mat)
  #getting VMR
  avg_lobe_expression_mat[['VMR']] = var/mean
  avg_lobe_expression_mat[['logVMR']] = log(avg_lobe_expression_mat$VMR)
  
  # getting spatially significant gene list
  
  vmrs =as.matrix(avg_lobe_expression_mat[,"logVMR"])
  # convert to density (empirical dist)
  d = density(na.omit(vmrs),bw = density.bandwidth)
  ts_y = ts(d$y)
  # inflection points of distribution
  tp = turnpoints(ts_y)
  idx = which.max(d$y[tp$tppos])
  # "mean" of distribution? 
  center=d$x[tp$tppos[idx]]
  belowcenter = vmrs[which(vmrs < center)]
  # building a null distribution?
  null = c(belowcenter,(belowcenter + 2*(center-belowcenter)))
  pvalsUpper=pnorm(vmrs,mean=center,sd= sd(null),lower.tail=F)
  # adjust p-values with Benjamini Hochberg correction 
  pvalsUpper=p.adjust(pvalsUpper,"fdr")
  sign_idx=which(pvalsUpper <= alpha)
  
  # get minimum significant vmr
  # cutoff_vmr = vmrs[sign_idx[which.min(vmrs[sign_idx])]]
  cutoff_vmr = min(vmrs[sign_idx])
  
  # plotting significant cutoff for vmr
  hist(vmrs)
  abline(v=cutoff_vmr)
  
  genelist<-rownames(avg_lobe_expression_mat)[which(avg_lobe_expression_mat$logVMR > cutoff_vmr)]
  # print(length(genelist))
  #removing genes that arent in norm.data
  liger_gene_list = lapply(liger_object@norm.data, FUN = rownames)
  liger_gene_list[[length(liger_gene_list) + 1]] = genelist
  final_list = Reduce(intersect, liger_gene_list)
  return(final_list)
}

########################################## new plotGene function (included in latest version of liger)
plotGene <- function(object, gene, use.raw = F, use.scaled = F, scale.by = 'dataset', 
                     methylation.indices = NULL, plot.by = 'dataset', set.dr.lims = F, 
                     pt.size = 0.1, min.clip = NULL, max.clip = NULL, clip.absolute = F, 
                     points.only = F, option = 'plasma', cols.use = NULL, zero.color = '#F5F5F5', 
                     axis.labels = NULL, do.legend = T, return.plots = F) {
  if ((plot.by != scale.by) & (use.scaled)) {
    warning("Provided values for plot.by and scale.by do not match; results may not be very
            interpretable.")
  }
  if (use.raw) {
    # drop only outer level names
    gene_vals <- getGeneValues(object@raw.data, gene)
  } else {
    if (use.scaled) {
      if (scale.by != 'dataset') {
        # check for feature 
        if (!(scale.by %in% colnames(object@cell.data)) & scale.by != 'none') {
          stop("Please select existing feature in cell.data to scale.by, or add it before calling.")
        }
        # have to rescale in this case 
        gene_vals <- getGeneValues(object@norm.data, gene)
        cellnames <- names(gene_vals)
        # set up dataframe with groups
        gene_df <- data.frame(gene = gene_vals)
        if (scale.by == 'none') {
          gene_df[['scaleby']] = 'none'
        } else {
          gene_df[['scaleby']] = factor(object@cell.data[[scale.by]])
        }
        # using dplyr
        gene_df1 <- gene_df %>%
          group_by(scaleby) %>%
          # scale by selected feature
          mutate_at(vars(-group_cols()), function(x) { scale(x, center = F)})
        gene_vals <- gene_df1$gene
        names(gene_vals) <- cellnames
      } else {
        # remember to use cols instead
        gene_vals <- getGeneValues(object@scale.data, gene, use.cols = T, log2scale = T)
      }
    } else {
      # using normalized data
      # indicate methylation indices here 
      gene_vals <- getGeneValues(object@norm.data, gene, methylation.indices = methylation.indices,
                                 log2scale = T)
    }
  }
  gene_vals[gene_vals == 0] <- NA
  dr_df <- data.frame(object@tsne.coords)
  rownames(dr_df) <- rownames(object@cell.data)
  dr_df$gene <- as.numeric(gene_vals[rownames(dr_df)])
  colnames(dr_df) <- c("dr1", "dr2", "gene")
  # get dr limits for later
  lim1 <- c(min(dr_df$dr1), max(dr_df$dr1))
  lim2 <- c(min(dr_df$dr2), max(dr_df$dr2))
  
  if (plot.by != 'none') {
    if (!(plot.by %in% colnames(object@cell.data))) {
      stop("Please select existing feature in cell.data to plot.by, or add it before calling.")
    }
    dr_df$plotby <- factor(object@cell.data[[plot.by]])
  } else {
    dr_df$plotby <- factor("none")
  }
  # expand clip values if only single provided
  num_levels <- length(levels(dr_df$plotby))
  if (length(min.clip) == 1) {
    min.clip <- rep(min.clip, num_levels)
    names(min.clip) <- levels(dr_df$plotby)
  }
  if (length(max.clip) == 1) {
    max.clip <- rep(max.clip, num_levels)
    names(max.clip) <- levels(dr_df$plotby)
  }
  if (!is.null(min.clip) & is.null(names(min.clip))) {
    if (num_levels > 1) {
      message("Adding names to min.clip according to levels in plot.by group; order may not be 
              preserved as intended if multiple clip values passed in. Pass in named vector to 
              prevent this.")
    }
    names(min.clip) <- levels(dr_df$plotby)
    }
  if (!is.null(max.clip) & is.null(names(max.clip))) {
    if (num_levels > 1) {
      message("Adding names to max.clip according to levels in plot.by group; order may not be 
              preserved as intended if multiple clip values passed in. Pass in named vector to 
              prevent this.")
    }
    names(max.clip) <- levels(dr_df$plotby)
    }
  p_list <- list()
  for (sub_df in split(dr_df, f = dr_df$plotby)) {
    # maybe do quantile cutoff here
    group_name <- as.character(sub_df$plotby[1])
    if (!clip.absolute) {
      max_v <- quantile(sub_df$gene, probs = max.clip[group_name], na.rm = T)
      min_v <- quantile(sub_df$gene, probs = min.clip[group_name], na.rm = T)
    } else {
      max_v <- max.clip[group_name]
      min_v <- min.clip[group_name]
    }
    sub_df$gene[sub_df$gene < min_v & !is.na(sub_df$gene)] <- min_v
    sub_df$gene[sub_df$gene > max_v & !is.na(sub_df$gene)] <- max_v
    
    ggp <- ggplot(sub_df, aes(x = dr1, y = dr2, color = gene)) + geom_point(size = pt.size) +
      labs(col = gene)
    
    if (!is.null(cols.use)) {
      ggp <- ggp + scale_color_gradientn(colors = cols.use,
                                         na.value = zero.color)
    } else {
      ggp <- ggp + scale_color_viridis_c(option = option,
                                         direction = -1,
                                         na.value = zero.color)
    }
    if (set.dr.lims) {
      ggp <- ggp + xlim(lim1) + ylim(lim2)
    }
    
    if (plot.by != 'none') {
      base <- as.character(sub_df$plotby[1])
    } else {
      base <- ""
    }
    ggp <- ggp + ggtitle(base)
    
    if (!is.null(axis.labels)) {
      ggp <- ggp + xlab(axis.labels[1]) + ylab(axis.labels[2])
    }
    if (!do.legend) {
      ggp <- ggp + theme(legend.position = "none")
    }
    if (points.only) {
      ggp <- ggp + theme(
        axis.line = element_blank(), axis.text.x = element_blank(),
        axis.text.y = element_blank(), axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(), legend.position = "none",
        panel.background = element_blank(), panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.background = element_blank(), plot.title = element_blank()
      )
    }
    p_list[[as.character(sub_df$plotby[1])]] <- ggp
  }
  if (plot.by == 'dataset') {
    p_list <- p_list[names(object@raw.data)]
  }
  
  if (return.plots){
    if (length(p_list) == 1) {
      return(p_list[[1]])
    } else {
      return(p_list)
    }
  } else {
    for (plot in p_list) {
      print(plot)
    }
  }
  }

getGeneValues <- function(list, gene, use.cols = F, methylation.indices = NULL, log2scale = F,
                          scale.factor = 10000) {
  gene_vals <- unlist(lapply(seq_along(list), function(i) {
    mtx <- unname(list)[[i]]
    if (use.cols) {
      mtx <- t(mtx)
    } 
    if (gene %in% rownames(mtx)) {
      gene_vals_int <- mtx[gene, ]
    } else {
      gene_vals_int <- rep(list(0), ncol(mtx))
      names(gene_vals_int) <- colnames(mtx)
    }
    if (log2scale & !(i %in% methylation.indices)) {
      gene_vals_int <- log2(scale.factor * gene_vals_int + 1)
    } 
    return(gene_vals_int)
  }),
  use.names = T)
  return(gene_vals)
}

# clustering functions 
library(Seurat)
makeSNN = function (object, dims.use = 1:ncol(object@H.norm), 
                    k.param = 30, prune.SNN = 1/15, print.output = TRUE, 
                    force.recalc = FALSE, filename = NULL, 
                    save.SNN = TRUE, nn.eps = 0) 
{
  data.use = object@H.norm[, dims.use]
  parameters.to.store <- as.list(environment(), all = TRUE)[names(formals("makeSNN"))]
  parameters.to.store$object <- NULL
  parameters.to.store$print.output <- NULL
  if (!is.null(object@parameters$makeSNN) && !force.recalc) {
    old.parameters <- object@parameters$makeSNN
    old.parameters$time <- NULL
    old.parameters$print.output <- NULL
    if (all(all.equal(old.parameters, parameters.to.store) == 
            TRUE)) {
      warning("Build parameters exactly match those of already computed and stored SNN. To force recalculation, set force.recalc to TRUE.")
      return(object)
    }
  }
  parameters.to.store[['time']] = Sys.time()
  object@parameters[['makeSNN']] = parameters.to.store
  n.cells <- nrow(x = data.use)
  if (n.cells < k.param) {
    warning("k.param set larger than number of cells. Setting k.param to number of cells - 1.")
    k.param <- n.cells - 1
  }
  
  if (print.output) {
    cat("Computing nearest neighbor graph\n", file = stderr())
  }
  my.knn <- RANN::nn2(data = data.use, k = k.param, searchtype = "standard", 
                      eps = nn.eps)
  nn.ranked <- my.knn$nn.idx
  
  if (print.output) {
    cat("Computing SNN\n", file = stderr())
  }
  if (save.SNN | is.null(filename)) {
    object@agg.data[['snn']] <- Seurat:::ComputeSNN(nn_ranked = nn.ranked, prune = prune.SNN)
    rownames(object@agg.data$snn) <- rownames(object@cell.data)
    colnames(object@agg.data$snn) <- rownames(object@cell.data)
    if (!is.null(filename)) {
      WriteEdgeFile(snn = object@snn, filename = filename, 
                    display_progress = print.output)
    }
  }
  else {
    DirectSNNToFile(nn_ranked = nn.ranked, prune = prune.SNN, 
                    display_progress = print.output, filename = filename)
  }
  return(object)
}

calcClusters = function(object, dims.use = 1:ncol(object@H.norm), 
                        k.param = 30, prune.SNN = 1/15, print.output = FALSE, 
                        save.SNN = TRUE, reuse.SNN = TRUE, 
                        force.recalc = FALSE, nn.eps = 0, modularity.fxn = 1, resolution = 0.8, 
                        algorithm = 1, n.start = 10, n.iter = 10, random.seed = 0, 
                        temp.file.location = NULL, edge.file.name = NULL) {
  snn.built <- FALSE
  if (!is.null(object@agg.data$snn)) {
    snn.built <- TRUE
    save.SNN <- TRUE
  }
  if ((missing(x = dims.use) && missing(x = k.param) && 
       missing(x = prune.SNN) && 
       snn.built) || reuse.SNN) {
    save.SNN <- TRUE
    if (reuse.SNN && !snn.built) {
      stop("No SNN stored to reuse.")
    }
    if (reuse.SNN && (!missing(x = dims.use) || 
                      !missing(x = k.param) || !missing(x = prune.SNN))) {
      warning("SNN was not be rebuilt with new parameters. Continued with stored\n               SNN. To suppress this warning, remove all SNN building parameters.")
    }
  }
  else {
    object <- makeSNN(object = object, dims.use = dims.use, 
                      k.param = k.param, prune.SNN = prune.SNN, 
                      print.output = print.output,  
                      force.recalc = force.recalc, filename = edge.file.name, 
                      save.SNN = save.SNN, nn.eps = nn.eps)
  }
  for (r in resolution) {
    parameters.to.store <- as.list(environment(), all = TRUE)[names(formals("calcClusters"))]
    parameters.to.store$resolution <- r
    if (!is.null(object@parameters[[paste0("calcClusters.res.", 
                                           r)]])  & force.recalc != TRUE) {
      parameters.to.store$object <- NULL
      parameters.to.store$print.output <- NULL
      old.parameters <- object@parameters[[paste0("calcClusters.res.", r)]]
      old.parameters$time <- NULL
      old.parameters$print.output <- NULL
      if (all(all.equal(old.parameters, parameters.to.store) == 
              TRUE)) {
        warning(paste0("Clustering parameters for resolution ", 
                       r, " exactly match those of already computed. \n  To force recalculation, set force.recalc to TRUE."))
        next
      }
    }
    parameters.to.store$object <- NULL
    parameters.to.store$print.output <- NULL
    parameters.to.store[['time']] = Sys.time()
    object@parameters[[paste0("calcClusters.res.", r)]] = parameters.to.store
    ids <- Seurat:::RunModularityClustering(SNN = object@agg.data$snn, modularity = modularity.fxn, 
                                            resolution = r, algorithm = algorithm, n.start = n.start, 
                                            n.iter = n.iter, random.seed = random.seed, print.output = print.output, 
                                            temp.file.location = temp.file.location, edge.file.name = edge.file.name)
    object@clusters = factor(ids)
    names(object@clusters) = rownames(object@cell.data)
    object <- groupSingletons(object = object, SNN = object@agg.data$snn)
    name <- paste0("res.", r)
    object@cell.data[[name]] = object@clusters
  }
  if (!save.SNN) {
    object@agg.data <- NULL
    object@parameters$makeSNN <- NULL
  }
  return(object)
}

groupSingletons = function (object, SNN) 
{
  singletons <- c()
  c_table = table(object@clusters)
  cluster_names <- unique(x = object@clusters)
  for (cluster in cluster_names) {
    if (c_table[as.character(cluster)] == 1) {
      singletons <- append(singletons, cluster)
    }
  }
  cluster_names <- setdiff(x = cluster_names, y = singletons)
  connectivity <- vector(mode = "numeric", length = length(x = cluster_names))
  names(x = connectivity) <- cluster_names
  for (i in singletons) {
    for (j in cluster_names) {
      subSNN = SNN[getClusterCells(object, cluster = i), 
                   match(getClusterCells(object, cluster = j), 
                         table = colnames(x = SNN))]
      if (is.object(x = subSNN)) {
        connectivity[j] <- sum(subSNN)/(nrow(x = subSNN) * 
                                          ncol(x = subSNN))
      }
      else {
        connectivity[j] <- mean(x = subSNN)
      }
    }
    m <- max(connectivity, na.rm = T)
    mi <- which(x = connectivity == m, arr.ind = TRUE)
    closest_cluster <- sample(x = names(x = connectivity[mi]), 
                              1)
    object@clusters[object@clusters == i] = closest_cluster
  }
  if (length(x = singletons) > 0) {
    message(paste(length(x = singletons), "singletons identified.", 
                  length(x = unique(object@clusters)), "final clusters."))
  }
  object@clusters = droplevels(object@clusters)
  return(object)
}
getClusterCells = function(object, cluster) {
  return(names(object@clusters)[which(object@clusters == cluster)])
}


# special plot gene to do both legend and points
# also change plotting order so that points with positive expression are plotted on top
plotGene_2020 <- function(object, gene, use.raw = F, use.scaled = F, scale.by = 'dataset', 
                          log2scale = NULL, methylation.indices = NULL, plot.by = 'dataset', 
                          set.dr.lims = F, pt.size = 0.1, min.clip = NULL, max.clip = NULL, 
                          clip.absolute = F, points.only = F, option = 'plasma', cols.use = NULL, 
                          zero.color = '#F5F5F5', axis.labels = NULL, do.legend = T, return.plots = F) {
  if ((plot.by != scale.by) & (use.scaled)) {
    warning("Provided values for plot.by and scale.by do not match; results may not be very
            interpretable.")
  }
  if (use.raw) {
    if (is.null(log2scale)) {
      log2scale <- FALSE
    }
    # drop only outer level names
    gene_vals <- getGeneValues(object@raw.data, gene, log2scale = log2scale)
  } else {
    if (is.null(log2scale)) {
      log2scale <- TRUE
    }
    # rescale in case requested gene not highly variable
    if (use.scaled) {
      # check for feature 
      if (!(scale.by %in% colnames(object@cell.data)) & scale.by != 'none') {
        stop("Please select existing feature in cell.data to scale.by, or add it before calling.")
      }
      gene_vals <- getGeneValues(object@norm.data, gene)
      cellnames <- names(gene_vals)
      # set up dataframe with groups
      gene_df <- data.frame(gene = gene_vals)
      if (scale.by == 'none') {
        gene_df[['scaleby']] = 'none'
      } else {
        gene_df[['scaleby']] = factor(object@cell.data[[scale.by]])
      }
      gene_df1 <- gene_df %>%
        group_by(scaleby) %>%
        # scale by selected feature
        mutate_at(vars(-group_cols()), function(x) { scale(x, center = F)})
      gene_vals <- gene_df1$gene
      names(gene_vals) <- cellnames
      if (log2scale) {
        gene_vals <- log2(10000 * gene_vals + 1)
      }
    } else {
      # using normalized data
      # indicate methylation indices here 
      gene_vals <- getGeneValues(object@norm.data, gene, methylation.indices = methylation.indices,
                                 log2scale = log2scale)
    }
  }
  gene_vals[gene_vals == 0] <- NA
  dr_df <- data.frame(object@tsne.coords)
  rownames(dr_df) <- rownames(object@cell.data)
  dr_df$gene <- as.numeric(gene_vals[rownames(dr_df)])
  colnames(dr_df) <- c("dr1", "dr2", "gene")
  
  # put points with no expression in the back
  dr_df[['back']] = is.na(dr_df$gene)
  dr_df = dr_df[order(dr_df$back, decreasing = T),]
  
  # get dr limits for later
  lim1 <- c(min(dr_df$dr1), max(dr_df$dr1))
  lim2 <- c(min(dr_df$dr2), max(dr_df$dr2))
  
  if (plot.by != 'none') {
    if (!(plot.by %in% colnames(object@cell.data))) {
      stop("Please select existing feature in cell.data to plot.by, or add it before calling.")
    }
    dr_df$plotby <- factor(object@cell.data[[plot.by]])
  } else {
    dr_df$plotby <- factor("none")
  }
  # expand clip values if only single provided
  num_levels <- length(levels(dr_df$plotby))
  if (length(min.clip) == 1) {
    min.clip <- rep(min.clip, num_levels)
    names(min.clip) <- levels(dr_df$plotby)
  }
  if (length(max.clip) == 1) {
    max.clip <- rep(max.clip, num_levels)
    names(max.clip) <- levels(dr_df$plotby)
  }
  if (!is.null(min.clip) & is.null(names(min.clip))) {
    if (num_levels > 1) {
      message("Adding names to min.clip according to levels in plot.by group; order may not be 
              preserved as intended if multiple clip values passed in. Pass in named vector to 
              prevent this.")
    }
    names(min.clip) <- levels(dr_df$plotby)
    }
  if (!is.null(max.clip) & is.null(names(max.clip))) {
    if (num_levels > 1) {
      message("Adding names to max.clip according to levels in plot.by group; order may not be 
              preserved as intended if multiple clip values passed in. Pass in named vector to 
              prevent this.")
    }
    names(max.clip) <- levels(dr_df$plotby)
    }
  p_list <- list()
  for (sub_df in split(dr_df, f = dr_df$plotby)) {
    # maybe do quantile cutoff here
    group_name <- as.character(sub_df$plotby[1])
    if (!clip.absolute) {
      max_v <- quantile(sub_df$gene, probs = max.clip[group_name], na.rm = T)
      min_v <- quantile(sub_df$gene, probs = min.clip[group_name], na.rm = T)
    } else {
      max_v <- max.clip[group_name]
      min_v <- min.clip[group_name]
    }
    sub_df$gene[sub_df$gene < min_v & !is.na(sub_df$gene)] <- min_v
    sub_df$gene[sub_df$gene > max_v & !is.na(sub_df$gene)] <- max_v
    
    ggp <- ggplot(sub_df, aes(x = dr1, y = dr2, color = gene)) + geom_point(size = pt.size) +
      labs(col = gene)
    
    if (!is.null(cols.use)) {
      ggp <- ggp + scale_color_gradientn(colors = cols.use,
                                         na.value = zero.color)
    } else {
      ggp <- ggp + scale_color_viridis_c(option = option,
                                         direction = -1,
                                         na.value = zero.color)
    }
    if (set.dr.lims) {
      ggp <- ggp + xlim(lim1) + ylim(lim2)
    }
    
    if (plot.by != 'none') {
      base <- as.character(sub_df$plotby[1])
    } else {
      base <- ""
    }
    ggp <- ggp + ggtitle(base)
    
    if (!is.null(axis.labels)) {
      ggp <- ggp + xlab(axis.labels[1]) + ylab(axis.labels[2])
    }
    if (!do.legend) {
      ggp <- ggp + theme(legend.position = "none")
    }
    if (points.only) {
      ggp <- ggp + theme(
        axis.line = element_blank(), axis.text.x = element_blank(),
        axis.text.y = element_blank(), axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.background = element_blank(), panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.background = element_blank(), plot.title = element_blank()
      )
    }
    p_list[[as.character(sub_df$plotby[1])]] <- ggp
  }
  if (plot.by == 'dataset') {
    p_list <- p_list[names(object@raw.data)]
  }
  
  if (return.plots){
    if (length(p_list) == 1) {
      return(p_list[[1]])
    } else {
      return(p_list)
    }
  } else {
    for (plot in p_list) {
      print(plot)
    }
  }
}

