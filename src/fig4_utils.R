# Additional functions needed for figure 4 and spatial divergence analysis
library(dplyr)
library(tidyr)
library(reshape2)


# Variation of liger's plotGene to do both legend and points
# also change plotting order so that points with positive expression are plotted on top
plotGene_ordered <- function(object, gene, use.raw = F, use.scaled = F, scale.by = 'dataset', 
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
  
  # put points with no expression in the back
  dr_df[['back']] = is.na(dr_df$gene)
  dr_df = dr_df[order(dr_df$back, decreasing = T),]
  
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