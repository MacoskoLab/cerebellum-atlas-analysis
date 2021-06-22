# Additional functions for plotting in monocle
library(dplyr)

# Variation of plot_cells which allows for different color maps 
# and gene scaling options 
# Required to be able to easily plot trajectory graph and scaled genes on single figure
plot_cells_colors = function (cds, x = 1, y = 2, reduction_method = c("UMAP", "tSNE", 
                                                                   "PCA", "LSI", "Aligned"), 
                           color_cells_by = "cluster", group_cells_by = c("cluster", 
                                                                          "partition"), genes = NULL, show_trajectory_graph = TRUE, 
                           trajectory_graph_color = "grey28", trajectory_graph_segment_size = 0.75, 
                           norm_method = c("log", "size_only"), label_cell_groups = TRUE, 
                           label_groups_by_cluster = TRUE, group_label_size = 2, labels_per_group = 1, 
                           label_branch_points = TRUE, label_roots = TRUE, label_leaves = TRUE, 
                           graph_label_size = 2, cell_size = 0.35, cell_stroke = I(cell_size/2), 
                           alpha = 1, min_expr = 0, rasterize = FALSE, scale_to_range = FALSE,
                           viridis_option = 'plasma', na.value = '#F5F5F5', direction = -1,
                           use_gradient = F, low_col = 'grey', high_col = 'blue',
                           normalized_values = NULL, scale.factor = 10000) 
{
  reduction_method <- match.arg(reduction_method)
  assertthat::assert_that(methods::is(cds, "cell_data_set"))
  assertthat::assert_that(!is.null(reducedDims(cds)[[reduction_method]]), 
                          msg = paste("No dimensionality reduction for", reduction_method, 
                                      "calculated.", "Please run reduce_dimension with", 
                                      "reduction_method =", reduction_method, "before attempting to plot."))
  low_dim_coords <- reducedDims(cds)[[reduction_method]]
  assertthat::assert_that(ncol(low_dim_coords) >= max(x, y), 
                          msg = paste("x and/or y is too large. x and y must", 
                                      "be dimensions in reduced dimension", "space."))
  if (!is.null(color_cells_by)) {
    assertthat::assert_that(color_cells_by %in% c("cluster", 
                                                  "partition", "pseudotime") | color_cells_by %in% 
                              names(colData(cds)), msg = paste("color_cells_by must one of", 
                                                               "'cluster', 'partition', 'pseudotime,", "or a column in the colData table."))
    if (color_cells_by == "pseudotime") {
      tryCatch({
        pseudotime(cds, reduction_method = reduction_method)
      }, error = function(x) {
        stop(paste("No pseudotime for", reduction_method, 
                   "calculated. Please run order_cells with", 
                   "reduction_method =", reduction_method, "before attempting to color by pseudotime."))
      })
    }
  }
  assertthat::assert_that(!is.null(color_cells_by) || !is.null(markers), 
                          msg = paste("Either color_cells_by or markers must", 
                                      "be NULL, cannot color by both!"))
  norm_method = match.arg(norm_method)
  group_cells_by = match.arg(group_cells_by)
  assertthat::assert_that(!is.null(color_cells_by) || !is.null(genes), 
                          msg = paste("Either color_cells_by or genes must be", 
                                      "NULL, cannot color by both!"))
  if (show_trajectory_graph && is.null(principal_graph(cds)[[reduction_method]])) {
    message("No trajectory to plot. Has learn_graph() been called yet?")
    show_trajectory_graph = FALSE
  }
  gene_short_name <- NA
  sample_name <- NA
  data_dim_1 <- NA
  data_dim_2 <- NA
  if (rasterize) {
    plotting_func <- ggrastr::geom_point_rast
  }
  else {
    plotting_func <- ggplot2::geom_point
  }
  S_matrix <- reducedDims(cds)[[reduction_method]]
  data_df <- data.frame(S_matrix[, c(x, y)])
  colnames(data_df) <- c("data_dim_1", "data_dim_2")
  data_df$sample_name <- row.names(data_df)
  data_df <- as.data.frame(cbind(data_df, colData(cds)))
  if (group_cells_by == "cluster") {
    data_df$cell_group <- tryCatch({
      clusters(cds, reduction_method = reduction_method)[data_df$sample_name]
    }, error = function(e) {
      NULL
    })
  }
  else if (group_cells_by == "partition") {
    data_df$cell_group <- tryCatch({
      partitions(cds, reduction_method = reduction_method)[data_df$sample_name]
    }, error = function(e) {
      NULL
    })
  }
  else {
    stop("Error: unrecognized way of grouping cells.")
  }
  if (color_cells_by == "cluster") {
    data_df$cell_color <- tryCatch({
      clusters(cds, reduction_method = reduction_method)[data_df$sample_name]
    }, error = function(e) {
      NULL
    })
  }
  else if (color_cells_by == "partition") {
    data_df$cell_color <- tryCatch({
      partitions(cds, reduction_method = reduction_method)[data_df$sample_name]
    }, error = function(e) {
      NULL
    })
  }
  else if (color_cells_by == "pseudotime") {
    data_df$cell_color <- tryCatch({
      pseudotime(cds, reduction_method = reduction_method)[data_df$sample_name]
    }, error = function(e) {
      NULL
    })
  }
  else {
    data_df$cell_color <- colData(cds)[data_df$sample_name, 
                                       color_cells_by]
  }
  if (show_trajectory_graph) {
    ica_space_df <- t(cds@principal_graph_aux[[reduction_method]]$dp_mst) %>% 
      as.data.frame() %>% dplyr::select_(prin_graph_dim_1 = x, 
                                         prin_graph_dim_2 = y) %>% dplyr::mutate(sample_name = rownames(.), 
                                                                                 sample_state = rownames(.))
    dp_mst <- cds@principal_graph[[reduction_method]]
    edge_df <- dp_mst %>% igraph::as_data_frame() %>% dplyr::select_(source = "from", 
                                                                     target = "to") %>% dplyr::left_join(ica_space_df %>% 
                                                                                                           dplyr::select_(source = "sample_name", source_prin_graph_dim_1 = "prin_graph_dim_1", 
                                                                                                                          source_prin_graph_dim_2 = "prin_graph_dim_2"), 
                                                                                                         by = "source") %>% dplyr::left_join(ica_space_df %>% 
                                                                                                                                               dplyr::select_(target = "sample_name", target_prin_graph_dim_1 = "prin_graph_dim_1", 
                                                                                                                                                              target_prin_graph_dim_2 = "prin_graph_dim_2"), 
                                                                                                                                             by = "target")
  }
  markers_exprs <- NULL
  expression_legend_label <- NULL
  if (!is.null(genes)) {
    if (!is.null(dim(genes)) && dim(genes) >= 2) {
      markers = unlist(genes[, 1], use.names = FALSE)
    }
    else {
      markers = genes
    }
    markers_rowData <- rowData(cds)[(rowData(cds)$gene_short_name %in% 
                                       markers) | (row.names(rowData(cds)) %in% markers), 
                                    , drop = FALSE]
    markers_rowData <- as.data.frame(markers_rowData)
    if (nrow(markers_rowData) == 0) {
      stop("None of the provided genes were found in the cds")
    }
    if (nrow(markers_rowData) >= 1) {
      cds_exprs <- SingleCellExperiment::counts(cds)[row.names(markers_rowData), 
                                                     , drop = FALSE]
      cds_exprs <- Matrix::t(Matrix::t(cds_exprs)/size_factors(cds))
      if (!is.null(dim(genes)) && dim(genes) >= 2) {
        genes = as.data.frame(genes)
        row.names(genes) = genes[, 1]
        genes = genes[row.names(cds_exprs), ]
        agg_mat = as.matrix(aggregate_gene_expression(cds, 
                                                      genes, norm_method = norm_method, scale_agg_values = FALSE))
        markers_exprs = agg_mat
        markers_exprs <- reshape2::melt(markers_exprs)
        colnames(markers_exprs)[1:2] <- c("feature_id", 
                                          "cell_id")
        if (is.factor(genes[, 2])) 
          markers_exprs$feature_id = factor(markers_exprs$feature_id, 
                                            levels = levels(genes[, 2]))
        markers_exprs$feature_label <- markers_exprs$feature_id
        norm_method = "size_only"
        expression_legend_label = "Expression score"
      }
      else {
        cds_exprs@x = round(10000 * cds_exprs@x)/10000
        markers_exprs = matrix(cds_exprs, nrow = nrow(markers_rowData))
        colnames(markers_exprs) = colnames(SingleCellExperiment::counts(cds))
        row.names(markers_exprs) = row.names(markers_rowData)
        # Added this in to emulate liger plotting 
        if (!is.null(normalized_values)) {
          markers_exprs = as.matrix(normalized_values[row.names(markers_exprs), colnames(markers_exprs),
                                                      drop = F])
        }
        markers_exprs <- reshape2::melt(markers_exprs)
        colnames(markers_exprs)[1:2] <- c("feature_id", 
                                          "cell_id")
        markers_exprs <- merge(markers_exprs, markers_rowData, 
                               by.x = "feature_id", by.y = "row.names")
        if (is.null(markers_exprs$gene_short_name)) {
          markers_exprs$feature_label <- as.character(markers_exprs$feature_id)
        }
        else {
          markers_exprs$feature_label <- as.character(markers_exprs$gene_short_name)
        }
        markers_exprs$feature_label <- ifelse(is.na(markers_exprs$feature_label) | 
                                                !as.character(markers_exprs$feature_label) %in% 
                                                markers, as.character(markers_exprs$feature_id), 
                                              as.character(markers_exprs$feature_label))
        markers_exprs$feature_label <- factor(markers_exprs$feature_label, 
                                              levels = markers)
        if (norm_method == "size_only") 
          expression_legend_label = "Expression"
        else if (!is.null(normalized_values)) {
          expression_legend_label = 'liger_scaled_vals'
        }
        else expression_legend_label = "log10(Expression)"
      }
      if (scale_to_range) {
        markers_exprs = dplyr::group_by(markers_exprs, 
                                        feature_label) %>% dplyr::mutate(max_val_for_feature = max(value), 
                                                                         min_val_for_feature = min(value)) %>% dplyr::mutate(value = 100 * 
                                                                                                                               (value - min_val_for_feature)/(max_val_for_feature - 
                                                                                                                                                                min_val_for_feature))
        expression_legend_label = "% Max"
      }
    }
  }
  if (label_cell_groups && is.null(color_cells_by) == FALSE) {
    if (is.null(data_df$cell_color)) {
      if (is.null(genes)) {
        message(paste(color_cells_by, "not found in colData(cds), cells will", 
                      "not be colored"))
      }
      text_df = NULL
      label_cell_groups = FALSE
    }
    else {
      if (is.character(data_df$cell_color) || is.factor(data_df$cell_color)) {
        if (label_groups_by_cluster && is.null(data_df$cell_group) == 
            FALSE) {
          text_df = data_df %>% dplyr::group_by(cell_group) %>% 
            dplyr::mutate(cells_in_cluster = dplyr::n()) %>% 
            dplyr::group_by(cell_color, add = TRUE) %>% 
            dplyr::mutate(per = dplyr::n()/cells_in_cluster)
          median_coord_df = text_df %>% dplyr::summarize(fraction_of_group = dplyr::n(), 
                                                         text_x = stats::median(x = data_dim_1), text_y = stats::median(x = data_dim_2))
          text_df = suppressMessages(text_df %>% dplyr::select(per) %>% 
                                       dplyr::distinct())
          text_df = suppressMessages(dplyr::inner_join(text_df, 
                                                       median_coord_df))
          text_df = text_df %>% dplyr::group_by(cell_group) %>% 
            dplyr::top_n(labels_per_group, per)
        }
        else {
          text_df = data_df %>% dplyr::group_by(cell_color) %>% 
            dplyr::mutate(per = 1)
          median_coord_df = text_df %>% dplyr::summarize(fraction_of_group = dplyr::n(), 
                                                         text_x = stats::median(x = data_dim_1), text_y = stats::median(x = data_dim_2))
          text_df = suppressMessages(text_df %>% dplyr::select(per) %>% 
                                       dplyr::distinct())
          text_df = suppressMessages(dplyr::inner_join(text_df, 
                                                       median_coord_df))
          text_df = text_df %>% dplyr::group_by(cell_color) %>% 
            dplyr::top_n(labels_per_group, per)
        }
        text_df$label = as.character(text_df %>% dplyr::pull(cell_color))
      }
      else {
        message(paste("Cells aren't colored in a way that allows them to", 
                      "be grouped."))
        text_df = NULL
        label_cell_groups = FALSE
      }
    }
  }
  if (!is.null(markers_exprs) && nrow(markers_exprs) > 0) {
    data_df <- merge(data_df, markers_exprs, by.x = "sample_name", 
                     by.y = "cell_id")
    data_df$value <- with(data_df, ifelse(value > min_expr, 
                                          value, NA))
    ya_sub <- data_df[!is.na(data_df$value), ]
    na_sub <- data_df[is.na(data_df$value), ]
    if (norm_method == "size_only") {
      g <- ggplot(data = data_df, aes(x = data_dim_1, y = data_dim_2)) + 
        plotting_func(aes(data_dim_1, data_dim_2), size = I(cell_size), 
                      stroke = I(cell_stroke), color = na.value, 
                      alpha = alpha, data = na_sub) + plotting_func(aes(color = value), 
                                                                    size = I(cell_size), stroke = I(cell_stroke), 
                                                                    data = ya_sub[order(ya_sub$value), ]) + viridis::scale_color_viridis(option = viridis_option, 
                                                                                                                                         name = expression_legend_label, na.value = na.value, 
                                                                                                                                         end = 1, direction = direction, alpha = alpha) + guides(alpha = FALSE) + 
        facet_wrap(~feature_label)
    }
    # Added in this option to emulate liger gene plotting
    else if (!is.null(normalized_values)){
      ya_sub = ya_sub[order(ya_sub$value),]
      data_df = rbind(na_sub, ya_sub)
      data_df$value = log2(scale.factor * data_df$value + 1)
      g <- ggplot(data = data_df, aes(x = data_dim_1, y = data_dim_2, color = value)) + 
        plotting_func(size = I(cell_size), 
                      stroke = I(cell_stroke),
                      alpha = alpha) + 
        # plotting_func(aes(color = value), size = I(cell_size), stroke = I(cell_stroke), 
        #                                                             data = ya_sub, alpha = alpha) +  # [order(ya_sub$value), ]
        viridis::scale_color_viridis(option = viridis_option, 
                                     name = expression_legend_label, na.value = na.value, 
                                     end = 1, direction = direction, alpha = alpha) + guides(alpha = FALSE) + 
        facet_wrap(~feature_label)
    }
    else {
      g <- ggplot(data = data_df, aes(x = data_dim_1, y = data_dim_2)) + 
        plotting_func(aes(data_dim_1, data_dim_2), size = I(cell_size), 
                      stroke = I(cell_stroke), color = na.value, 
                      data = na_sub, alpha = alpha) + plotting_func(aes(color = log10(value + 
                                                                                        min_expr)), size = I(cell_size), stroke = I(cell_stroke), 
                                                                    data = ya_sub[order(ya_sub$value), ], alpha = alpha) + 
        viridis::scale_color_viridis(option = viridis_option, 
                                     name = expression_legend_label, na.value = na.value, 
                                     end = 1, direction = direction, alpha = alpha) + guides(alpha = FALSE) + 
        facet_wrap(~feature_label)
    }
  }
  else {
    g <- ggplot(data = data_df, aes(x = data_dim_1, y = data_dim_2))
    if (color_cells_by %in% c("cluster", "partition")) {
      if (is.null(data_df$cell_color)) {
        g <- g + geom_point(color = I("gray"), size = I(cell_size), 
                            stroke = I(cell_stroke), na.rm = TRUE, alpha = I(alpha))
        message(paste("cluster_cells() has not been called yet, can't", 
                      "color cells by cluster"))
      }
      else {
        g <- g + geom_point(aes(color = cell_color), 
                            size = I(cell_size), stroke = I(cell_stroke), 
                            na.rm = TRUE, alpha = alpha)
      }
      g <- g + guides(color = guide_legend(title = color_cells_by, 
                                           override.aes = list(size = 4)))
    }
    else if (class(data_df$cell_color) == "numeric") {
      g <- g + geom_point(aes(color = cell_color), size = I(cell_size), 
                          stroke = I(cell_stroke), na.rm = TRUE, alpha = alpha)
      g <- g + viridis::scale_color_viridis(name = color_cells_by, 
                                            option = viridis_option,
                                            direction = direction)
      if (use_gradient) {
        g <- g + scale_color_gradient(name = color_cells_by, 
                                      low = low_col,
                                      high = high_col)
      }
      
    }
    else {
      g <- g + geom_point(aes(color = cell_color), size = I(cell_size), 
                          stroke = I(cell_stroke), na.rm = TRUE, alpha = alpha)
      g <- g + guides(color = guide_legend(title = color_cells_by, 
                                           override.aes = list(size = 4)))
    }
  }
  if (show_trajectory_graph) {
    g <- g + geom_segment(aes_string(x = "source_prin_graph_dim_1", 
                                     y = "source_prin_graph_dim_2", xend = "target_prin_graph_dim_1", 
                                     yend = "target_prin_graph_dim_2"), size = trajectory_graph_segment_size, 
                          color = I(trajectory_graph_color), linetype = "solid", 
                          na.rm = TRUE, data = edge_df)
    if (label_branch_points) {
      mst_branch_nodes <- branch_nodes(cds, reduction_method)
      branch_point_df <- ica_space_df %>% dplyr::slice(match(names(mst_branch_nodes), 
                                                             sample_name)) %>% dplyr::mutate(branch_point_idx = seq_len(dplyr::n()))
      g <- g + geom_point(aes_string(x = "prin_graph_dim_1", 
                                     y = "prin_graph_dim_2"), shape = 21, stroke = I(trajectory_graph_segment_size), 
                          color = "white", fill = "black", size = I(graph_label_size * 
                                                                      1.5), na.rm = TRUE, branch_point_df) + geom_text(aes_string(x = "prin_graph_dim_1", 
                                                                                                                                  y = "prin_graph_dim_2", label = "branch_point_idx"), 
                                                                                                                       size = I(graph_label_size), color = "white", 
                                                                                                                       na.rm = TRUE, branch_point_df)
    }
    if (label_leaves) {
      mst_leaf_nodes <- leaf_nodes(cds, reduction_method)
      leaf_df <- ica_space_df %>% dplyr::slice(match(names(mst_leaf_nodes), 
                                                     sample_name)) %>% dplyr::mutate(leaf_idx = seq_len(dplyr::n()))
      g <- g + geom_point(aes_string(x = "prin_graph_dim_1", 
                                     y = "prin_graph_dim_2"), shape = 21, stroke = I(trajectory_graph_segment_size), 
                          color = "black", fill = "lightgray", size = I(graph_label_size * 
                                                                          1.5), na.rm = TRUE, leaf_df) + geom_text(aes_string(x = "prin_graph_dim_1", 
                                                                                                                              y = "prin_graph_dim_2", label = "leaf_idx"), 
                                                                                                                   size = I(graph_label_size), color = "black", 
                                                                                                                   na.rm = TRUE, leaf_df)
    }
    if (label_roots) {
      mst_root_nodes <- monocle3:::root_nodes(cds, reduction_method)
      root_df <- ica_space_df %>% dplyr::slice(match(names(mst_root_nodes), 
                                                     sample_name)) %>% dplyr::mutate(root_idx = seq_len(dplyr::n()))
      g <- g + geom_point(aes_string(x = "prin_graph_dim_1", 
                                     y = "prin_graph_dim_2"), shape = 21, stroke = I(trajectory_graph_segment_size), 
                          color = "black", fill = "white", size = I(graph_label_size * 
                                                                      1.5), na.rm = TRUE, root_df) + geom_text(aes_string(x = "prin_graph_dim_1", 
                                                                                                                          y = "prin_graph_dim_2", label = "root_idx"), 
                                                                                                               size = I(graph_label_size), color = "black", 
                                                                                                               na.rm = TRUE, root_df)
    }
  }
  if (label_cell_groups) {
    g <- g + ggrepel::geom_text_repel(data = text_df, mapping = aes_string(x = "text_x", 
                                                                           y = "text_y", label = "label"), size = I(group_label_size))
    if (is.null(markers_exprs)) 
      g <- g + theme(legend.position = "none")
  }
  g <- g + monocle3:::monocle_theme_opts() + xlab(paste(reduction_method, 
                                                        x)) + ylab(paste(reduction_method, y)) + theme(legend.key = element_blank()) + 
    theme(panel.background = element_rect(fill = "white"))
  g
}