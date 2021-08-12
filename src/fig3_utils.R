# Additional functions needed for figure 3
library(pbapply)

generateObjectContinuityFits = function(object, ident1, ident2, use.raw.order = T,
                                        reverse_order = F, subset.idents = F, genes_test = 300, rand.seed = 1, 
                                        return_monocle = F) {
  
  if (length(levels(object@ident)) >2 & subset.idents) {
    print('Subsetting object to relevant identities')
    object = SubsetData(object, ident.use = c(ident1, ident2), subset.raw = T)
  }
  
  if (ncol(object@raw.data) > 5000) {
    print("Downsampling object")
    if (use.raw.order) {
      set.seed(rand.seed)
      cells.use = sample(colnames(object@raw.data), size = 5000)
      object = SubsetData(object, cells.use = cells.use, 
                          subset.raw = T)
    } else {
      set.seed(rand.seed)
      cells.use = sample(colnames(object@raw.data), size = 5000)
      object = SubsetData(object, cells.use = cells.use, 
                          subset.raw = T)
    }
  }
  
  object = NormalizeData(object)
  print("Calculating differential expression")
  # compare only two idents
  genes1 = FindMarkers(object, ident.1 = ident1, ident.2 = ident2)
  genes2 = FindMarkers(object, ident.1 = ident2, ident.2 = ident1)
  
  genes1_sig = row.names(genes1[genes1$p_val_adj < 0.05,])[1:300]
  genes2_sig = row.names(genes2[genes2$p_val_adj < 0.05,])[1:300]
  
  gene_union = union(genes1_sig, genes2_sig)
  # remove missing genes
  if (any(is.na(gene_union))) {
    gene_union = gene_union[-which(is.na(gene_union))]
  }
  
  # do monocle based ordering 
  print("Creating monocle object")
  monocle_obj = importCDS(object)
  monocle_obj <- estimateSizeFactors(monocle_obj)
  monocle_obj <- estimateDispersions(monocle_obj)
  
  # set DEGs as ordering genes
  monocle_obj <- setOrderingFilter(monocle_obj, gene_union)
  
  monocle_obj <- reduceDimension(monocle_obj, max_components = 2,
                                 method = 'DDRTree')
  # order the cells 
  monocle_obj <- orderCells(monocle_obj)
  
  pseudo_vals = monocle_obj@phenoData@data$Pseudotime
  names(pseudo_vals) = row.names(monocle_obj@phenoData@data)
  mon_pseudo_order = names(sort(pseudo_vals))
  if (reverse_order) {
    mon_pseudo_order = rev(mon_pseudo_order)
  }
  
  # switch norm and scale data back for log fitting
  object@data = liger::Matrix.column_norm(object@raw.data)
  object@scale.data = scale(t(object@data), scale = T, center = F)
  
  print("Calculating sigmoid fits")
  # allows negative and positive offset
  log_params_mono_ca = calcLogitParams(object, factor_order = mon_pseudo_order, 
                                       genes_test = gene_union,
                                       use_scale = T)
  # positive offset only
  log_params_mono_ca_pos = calcLogitParams(object, factor_order = mon_pseudo_order, 
                                           genes_test = gene_union,
                                           use_scale = T, version = 2)
  
  negv_results_ca = row.names(log_params_mono_ca[log_params_mono_ca$isConv & log_params_mono_ca$v < 0,])
  ca_combined = log_params_mono_ca
  ca_combined[negv_results_ca,] = log_params_mono_ca_pos[negv_results_ca,]
  
  if (return_monocle) {
    return(list(fits = ca_combined, pseudo_order = mon_pseudo_order, object = object,
                monocle = monocle_obj))
  }
  
  return(list(fits = ca_combined, pseudo_order = mon_pseudo_order, object = object))
}


calcLogitParams = function(object, factor_order, genes_test = NULL, start_vals = NULL, use_scale = F,
                           version = 1) {
  if (class(object) == 'liger') {
    norm_data = liger:::MergeSparseDataAll(object@norm.data)
  } else {
    norm_data = object@data
  }
  if (use_scale) {
    norm_data = t(object@scale.data)
  }
  if (is.null(genes_test)) {
    genes_test = row.names(norm_data)
  }
  norm_data = norm_data[, factor_order]
  # scale rank order to 1
  x = 1:length(factor_order) / length(factor_order)
  # remove missing
  genes_test = genes_test[Matrix::rowSums(norm_data[genes_test,]) > 0]
  
  gene_params = pblapply(seq_along(genes_test), FUN = function(i) {
    g = genes_test[i]
    y = norm_data[g,]
    # reverse expression if negative correlation
    cor_val = cor(x, y, method = 'spearman')
    if (cor_val < 0) {
      y = rev(y)
      reversed = 1
    } else {
      reversed = 0
    }
    fitmodel_list = fitLogitParams(x, y, start_vals = start_vals, offset_version = version)
    fitmodel = fitmodel_list[[1]]
    log_params = coef(fitmodel)
    if (length(log_params) == 4) {
      vert_shift = log_params[4]
    } else {
      vert_shift = 0
    }
    return(as.data.frame(list(log_params[1], log_params[2], log_params[3], vert_shift, reversed, cor_val, 
                              fitmodel$convInfo$isConv, fitmodel$convInfo$stopMessage,
                              fitmodel_list[[2]], fitmodel_list[[3]]), stringsAsFactors = F,
                         col.names = c('a', 'b', 'c', 'v', 'reversed', 'cor_val', 'isConv', 'stopMessage', 'SRS', 'model_used')))
  }) 
  full_set = do.call(rbind, gene_params)
  row.names(full_set) = genes_test
  return(full_set)
}

# Fit sigmoidal curve to gene expression

fitLogitParams = function(x, y, start_vals = NULL, offset_version = 1) {
  # default values to test based on initial experiments 
  # small grid of values
  if (is.null(start_vals)) {
    start_vals = list(
      list(a=2, b = 6, c = 0.3, v = 0),
      list(a = 20, b = 6, c = 1.3, v = 0),
      list(a = 2, b = 4, c = 2, v = 0),
      list(a = 2, b = 100, c = 0.75, v = 0)
    )
  } else {
    # pass custom grid of values
    start_vals = list(list(a = start_vals[1], b = start_vals[2], c = start_vals[3], v = start_vals[4]))
  }
  df = data.frame(index = x, gene = y + 1e-15)
  
  deviances = c()
  fitted_models = list()
  for (i in 1:length(start_vals)) {
    # first version allows negative or positive offset
    # version 2 is positive only
    if (offset_version == 1) {
      fitmodel1 <- nls(gene ~ a/(1 + exp(-b * (index-c))) + v, data=df, start = start_vals[[i]],
                       control = nls.control(maxiter = 250, tol = 1e-05, minFactor = 1e-10,
                                             printEval = FALSE, warnOnly = TRUE))
    } else {
      fitmodel1 <- nls(gene ~ a/(1 + exp(-b * (index-c))) + v^2, data=df, start = start_vals[[i]],
                       control = nls.control(maxiter = 250, tol = 1e-05, minFactor = 1e-10,
                                             printEval = FALSE, warnOnly = TRUE))
    }
    rss1 = fitmodel1$m$deviance()
    
    # same as model above but without any vertical shift
    fitmodel1_5 <- nls(gene ~ a/(1 + exp(-b * (index-c))), data=df, start = start_vals[[i]][1:3],
                       control = nls.control(maxiter = 250, tol = 1e-05, minFactor = 1e-10,
                                             printEval = FALSE, warnOnly = TRUE))
    rss1_5 = fitmodel1_5$m$deviance()
    
    # store the values -- index + half is naming convention
    fitted_models[[i*2 - 1]] = list(fitmodel1, rss1, i)
    fitted_models[[i*2]] = list(fitmodel1_5, rss1_5, i + 0.5)
    deviances = c(deviances, c(rss1, rss1_5))
  }
  
  # models with .5 are the ones without offsets 
  best_mod = which.min(deviances)
  return(fitted_models[[best_mod]])
}

# Aggregate continuity results, with options for how to rank gene results

makeCompareDF2 = function(list_param_dfs, names, n_cell_counts, spearman_quantile = 0.7, 
         spearman_cutoff = NULL, top_ranked_spearman = NULL, converged_only = F, positive_v = F,
         top_ranked_pval = NULL) {
  subset_list = lapply(seq_along(list_param_dfs), function(i) {
    x = list_param_dfs[[i]]
    quant = quantile(abs(x$cor_val), probs = spearman_quantile)
    if (!is.null(spearman_cutoff)) {
      quant = spearman_cutoff[i]
    }

    # original slope is a*b/4
    x[['mid_slope']] = x$b / 4
    if (converged_only) {
      x = x[x$isConv == TRUE,]
    }
    # exclude genes with negative vertical offset values 
    if (positive_v)  {
      x = x[x$v >= 0,]
    }
    if (!is.null(top_ranked_spearman)) {
      x[['abs_cor_val']] = abs(x$cor_val)
      subset1 = x[order(x$abs_cor_val, decreasing = T),][1:top_ranked_spearman,]
    } else if (!is.null(top_ranked_pval)) {
      # already ordered by decreasing pvalue
      subset1 = x[1:top_ranked_pval,]
    } else {
      subset1 = x[abs(x$cor_val) > quant,]
    }
    subset1[['dataset']] = names[i]
    subset1[['n_cell_counts']] = n_cell_counts[i]
    subset1[['gene']] = row.names(subset1)
    # c parameter corresponds to the midpoint x value 
    subset1[c('b', 'c', 'mid_slope', 'cor_val', 'gene', 'dataset', 'n_cell_counts')]
  })
  
  compare_df = do.call(rbind, subset_list)
  compare_df[['log_b']] = log10(compare_df$b)
  compare_df[['log_midslope']] = log10(compare_df$mid_slope)

  compare_df[['corrected_c']] = sapply(1:nrow(compare_df), function(i) {
    if (compare_df[i, 'cor_val'] < 0) {
      compare_df[i, 'n_cell_counts'] - compare_df[i, 'c']
    } else {
      compare_df[i, 'c']
    }
  })
  
  return(compare_df)
}

# Plotting functionality

sigmoid = function(params, x) {
  params[1] / (1 + exp(-params[2] * (x - params[3]))) + params[4]
}

plotFittedGene = function(object, gene, factor_order, fitparams, use_scale = F, single_params = FALSE,
                          cells_highlight = NULL, pt.size = 0.1, linesize = 1, color_by_cluster = T, 
                          continuous_color = NULL, use_gradient = F, low_col = 'grey', high_col = 'orangered',
                          do.legend = T, use_step = F, return.plot = F, return.data = F) {
  if (class(object) == 'liger') {
    norm_data = liger:::MergeSparseDataAll(object@norm.data)
  } else {
    norm_data = object@data
  }
  if (use_scale) {
    norm_data = t(object@scale.data)
  }
  norm_data = norm_data[, factor_order]
  # scale to same overall length
  x = 1:length(factor_order) / length(factor_order)
  
  fit_df = data.frame(index = x, gene = norm_data[gene,])
  # print(any(fit_df$gene < 0))
  fit_df[['cluster']] = object@ident[row.names(fit_df)]
  # must be column in meta.data slot
  if (!is.null(continuous_color)) {
    fit_df[['continuous']] = object@meta.data[row.names(fit_df), continuous_color]
  }
  
  if (single_params) {
    if (fitparams[5] == 1) {
      reversed = TRUE
      print('Order was reversed as negative correlation')
      fit_df$gene = rev(norm_data[gene,])
      row.names(fit_df) = rev(row.names(fit_df))
      fit_df$cluster = rev(fit_df$cluster)
      if (!is.null(continuous_color))
        fit_df$continuous = rev(fit_df$continuous)
      # print(head(fit_df))
    } else {
      reversed = FALSE
    }
    fit_df["fit_values"] = sigmoid(fitparams[1:4], x)
    if (use_step) {
      fit_df['fit_values'] = step_function(fitparams[1:4], x)
    }
    
  } else {
    if (fitparams[gene,][5] == 1) {
      print('Order was reversed as negative correlation')
      reversed = TRUE
      fit_df$gene = rev(norm_data[gene,])
      row.names(fit_df) = rev(row.names(fit_df))
      fit_df$cluster = rev(fit_df$cluster)
      if (!is.null(continuous_color))
        fit_df$continuous = rev(fit_df$continuous)
      # print(head(fit_df))
    } else {
      reversed = FALSE
    }
    
    fit_df["fit_values"] = sigmoid(as.numeric(fitparams[gene, 1:4]), x)
    # subtract vertical displacement and divide by range 
    fit_df["fit_values"] = (sigmoid(as.numeric(fitparams[gene, 1:4]), x) - as.numeric(fitparams[gene, 4])) /
      as.numeric(fitparams[gene, 1])
    if (use_step) {
      fit_df['fit_values'] = step_function(as.numeric(fitparams[gene, 1:4]), x)
    }
  }
  
  fit_df$gene = (fit_df$gene - as.numeric(fitparams[gene, 4])) / as.numeric(fitparams[gene, 1])
  
  if (!is.null(cells_highlight)) {
    fit_df[['highlight']] = FALSE
    fit_df[cells_highlight, 'highlight'] = TRUE
    # print(head(fit_df[fit_df$highlight == TRUE,]))
    plot1 = ggplot(fit_df, aes(x = index, y = gene, col = highlight)) + geom_point(size = pt.size) + 
      geom_line(aes(x = index, y = fit_values),col = 'red', size = linesize) + 
      labs(x= '', y = 'Normalized Gene Expression')+ ggtitle(paste0(gene, ' (factor ordered)'))
    plot1 = plot1 + scale_color_manual(values = c('black', 'orange')) + theme(legend.position = "none")
  } else {
    plot1 = ggplot(fit_df, aes(x = index, y = gene)) + geom_point(size = pt.size) + 
      geom_line(aes(x = index, y = fit_values),col = 'red', size = linesize) + 
      labs(x= '', y = 'Normalized Gene Expression')+ ggtitle(paste0(gene, ' (factor ordered)'))
  }
  if (color_by_cluster) {
    plot1 = ggplot(fit_df, aes(x = index, y = gene, col = cluster)) + geom_point(size = pt.size) + 
      geom_line(aes(x = index, y = fit_values), col = 'red', size = linesize) + 
      labs(x= '', y = 'Normalized Gene Expression')+ ggtitle(paste0(gene, ' (factor ordered)'))
  } else if (!is.null(continuous_color)) {
    plot1 = ggplot(fit_df, aes(x = index, y = gene, col = continuous)) + geom_point(size = pt.size) + 
      geom_line(aes(x = index, y = fit_values),col = 'red', size = linesize) + 
      scale_color_viridis_c(direction = -1, option = 'plasma') +
      labs(x= '', y = 'Normalized Gene Expression')+ ggtitle(paste0(gene, ' (factor ordered)'))
    if (use_gradient) {
      plot1 = plot1 + scale_color_gradient(low = low_col,
                                           high = high_col)
    }
  }
  if (!do.legend) {
    plot1 = plot1 + theme(legend.position = 'none')
  }
  if (reversed) {
    plot1 = plot1 + scale_x_reverse()
  }
  if (return.plot) {
    return(plot1)
  }
  if (return.data) {
    return(fit_df)
  }
  print(plot1)
  
}

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
