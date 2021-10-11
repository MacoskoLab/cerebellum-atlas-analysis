# Reproducing Extended Data Fig 1 panels
library(ggplot2)
library(Seurat)
# v2.3.4
library(liger)
# v0.4.1
source('../cb_custom_functions.R')

# Define utility function (based on DimPlot in Seurat)
# new version of DimPlot for shuffled points
DimPlot2 = function (object, reduction.use = "pca", dim.1 = 1, dim.2 = 2, 
                     cells.use = NULL, pt.size = 1, do.return = FALSE, do.bare = FALSE, 
                     cols.use = NULL, group.by = "ident", pt.shape = NULL, do.hover = FALSE, 
                     data.hover = "ident", do.identify = FALSE, do.label = FALSE, 
                     label.size = 4, no.legend = FALSE, coord.fixed = FALSE, no.axes = FALSE, 
                     dark.theme = FALSE, plot.order = NULL, cells.highlight = NULL, 
                     cols.highlight = "red", sizes.highlight = 1, plot.title = NULL, 
                     vector.friendly = FALSE, shuffle_points = F, png.file = NULL, png.arguments = c(10, 
                                                                                                     10, 100), na.value = "grey50", ...) 
{
  if (vector.friendly) {
    previous_call <- blank_call <- png_call <- match.call()
    blank_call$pt.size <- -1
    blank_call$do.return <- TRUE
    blank_call$vector.friendly <- FALSE
    png_call$no.axes <- TRUE
    png_call$no.legend <- TRUE
    png_call$do.return <- TRUE
    png_call$vector.friendly <- FALSE
    png_call$plot.title <- NULL
    blank_plot <- eval(blank_call, sys.frame(sys.parent()))
    png_plot <- eval(png_call, sys.frame(sys.parent()))
    png.file <- SetIfNull(x = png.file, default = paste0(tempfile(), 
                                                         ".png"))
    ggsave(filename = png.file, plot = png_plot, width = png.arguments[1], 
           height = png.arguments[2], dpi = png.arguments[3])
    to_return <- AugmentPlot(plot1 = blank_plot, imgFile = png.file)
    file.remove(png.file)
    if (do.return) {
      return(to_return)
    }
    else {
      print(to_return)
    }
  }
  embeddings.use <- GetDimReduction(object = object, reduction.type = reduction.use, 
                                    slot = "cell.embeddings")
  if (length(x = embeddings.use) == 0) {
    stop(paste(reduction.use, "has not been run for this object yet."))
  }
  cells.use <- SetIfNull(x = cells.use, default = colnames(x = object@data))
  dim.code <- GetDimReduction(object = object, reduction.type = reduction.use, 
                              slot = "key")
  dim.codes <- paste0(dim.code, c(dim.1, dim.2))
  data.plot <- as.data.frame(x = embeddings.use)
  cells.use <- intersect(x = cells.use, y = rownames(x = data.plot))
  data.plot <- data.plot[cells.use, dim.codes]
  ident.use <- as.factor(x = object@ident[cells.use])
  if (group.by != "ident") {
    ident.use <- as.factor(x = FetchData(object = object, 
                                         vars.all = group.by)[cells.use, 1])
  }
  data.plot$ident <- ident.use
  data.plot$x <- data.plot[, dim.codes[1]]
  data.plot$y <- data.plot[, dim.codes[2]]
  data.plot$pt.size <- pt.size
  if (!is.null(x = cells.highlight)) {
    if (is.character(x = cells.highlight)) {
      cells.highlight <- list(cells.highlight)
    }
    else if (is.data.frame(x = cells.highlight) || !is.list(x = cells.highlight)) {
      cells.highlight <- as.list(x = cells.highlight)
    }
    cells.highlight <- lapply(X = cells.highlight, FUN = function(cells) {
      cells.return <- if (is.character(x = cells)) {
        cells[cells %in% rownames(x = data.plot)]
      }
      else {
        cells <- as.numeric(x = cells)
        cells <- cells[cells <= nrow(x = data.plot)]
        rownames(x = data.plot)[cells]
      }
      return(cells.return)
    })
    cells.highlight <- Filter(f = length, x = cells.highlight)
    if (length(x = cells.highlight) > 0) {
      if (!no.legend) {
        no.legend <- is.null(x = names(x = cells.highlight))
      }
      names.highlight <- if (is.null(x = names(x = cells.highlight))) {
        paste0("Group_", 1L:length(x = cells.highlight))
      }
      else {
        names(x = cells.highlight)
      }
      sizes.highlight <- rep_len(x = sizes.highlight, length.out = length(x = cells.highlight))
      cols.highlight <- rep_len(x = cols.highlight, length.out = length(x = cells.highlight))
      highlight <- rep_len(x = NA_character_, length.out = nrow(x = data.plot))
      if (is.null(x = cols.use)) {
        cols.use <- "black"
      }
      cols.use <- c(cols.use[1], cols.highlight)
      size <- rep_len(x = pt.size, length.out = nrow(x = data.plot))
      for (i in 1:length(x = cells.highlight)) {
        cells.check <- cells.highlight[[i]]
        index.check <- match(x = cells.check, rownames(x = data.plot))
        highlight[index.check] <- names.highlight[i]
        size[index.check] <- sizes.highlight[i]
      }
      plot.order <- sort(x = unique(x = highlight), na.last = TRUE)
      plot.order[is.na(x = plot.order)] <- "Unselected"
      highlight[is.na(x = highlight)] <- "Unselected"
      highlight <- as.factor(x = highlight)
      data.plot$ident <- highlight
      data.plot$pt.size <- size
      if (dark.theme) {
        cols.use[1] <- "white"
      }
    }
  }
  # if (!is.null(x = plot.order)) {
  #   if (any(!plot.order %in% data.plot$ident)) {
  #     stop("invalid ident in plot.order")
  #   }
  #   plot.order <- rev(x = c(plot.order, setdiff(x = unique(x = data.plot$ident), 
  #                                               y = plot.order)))
  #   data.plot$ident <- factor(x = data.plot$ident, levels = plot.order)
  #   data.plot <- data.plot[order(data.plot$ident), ]
  # }
  if (shuffle_points) {
    data.plot = data.plot[sample(row.names(data.plot)),]
  }
  p <- ggplot(data = data.plot, mapping = aes(x = x, y = y)) + 
    geom_point(mapping = aes(colour = factor(x = ident), 
                             size = pt.size))
  if (!is.null(x = pt.shape)) {
    shape.val <- FetchData(object = object, vars.all = pt.shape)[cells.use, 
                                                                 1]
    if (is.numeric(shape.val)) {
      shape.val <- cut(x = shape.val, breaks = 5)
    }
    data.plot[, "pt.shape"] <- shape.val
    p <- ggplot(data = data.plot, mapping = aes(x = x, y = y)) + 
      geom_point(mapping = aes(colour = factor(x = ident), 
                               shape = factor(x = pt.shape), size = pt.size))
  }
  if (!is.null(x = cols.use)) {
    p <- p + scale_colour_manual(values = cols.use, na.value = na.value)
  }
  if (coord.fixed) {
    p <- p + coord_fixed()
  }
  p <- p + guides(size = FALSE)
  p2 <- p + xlab(label = dim.codes[[1]]) + ylab(label = dim.codes[[2]]) + 
    scale_size(range = c(min(data.plot$pt.size), max(data.plot$pt.size)))
  p3 <- p2 + Seurat:::SetXAxisGG() + Seurat:::SetYAxisGG() + Seurat:::SetLegendPointsGG(x = 6) + 
    Seurat:::SetLegendTextGG(x = 12) +Seurat:::no.legend.title + theme_bw() + 
    Seurat:::NoGrid()
  if (dark.theme) {
    p <- p + Seurat:::DarkTheme()
    p3 <- p3 + Seurat:::DarkTheme()
  }
  p3 <- p3 + theme(legend.title = element_blank())
  if (!is.null(plot.title)) {
    p3 <- p3 + ggtitle(plot.title) + theme(plot.title = element_text(hjust = 0.5))
  }
  if (do.label) {
    centers <- data.plot %>% dplyr::group_by(ident) %>% summarize(x = median(x = x), 
                                                                  y = median(x = y))
    p3 <- p3 + geom_point(data = centers, mapping = aes(x = x, 
                                                        y = y), size = 0, alpha = 0) + geom_text(data = centers, 
                                                                                                 mapping = aes(label = ident), size = label.size)
  }
  if (no.legend) {
    p3 <- p3 + theme(legend.position = "none")
  }
  if (no.axes) {
    p3 <- p3 + theme(axis.line = element_blank(), axis.text.x = element_blank(), 
                     axis.text.y = element_blank(), axis.ticks = element_blank(), 
                     axis.title.x = element_blank(), axis.title.y = element_blank(), 
                     panel.background = element_blank(), panel.border = element_blank(), 
                     panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                     plot.background = element_blank())
  }
  if (do.identify || do.hover) {
    if (do.bare) {
      plot.use <- p
    }
    else {
      plot.use <- p3
    }
    if (do.hover) {
      if (is.null(x = data.hover)) {
        features.info <- NULL
      }
      else {
        features.info <- FetchData(object = object, vars.all = data.hover)
      }
      return(HoverLocator(plot = plot.use, data.plot = data.plot, 
                          features.info = features.info, dark.theme = dark.theme))
    }
    else if (do.identify) {
      return(FeatureLocator(plot = plot.use, data.plot = data.plot, 
                            dark.theme = dark.theme, ...))
    }
  }
  if (do.return) {
    if (do.bare) {
      return(p)
    }
    else {
      return(p3)
    }
  }
  if (do.bare) {
    print(p)
  }
  else {
    print(p3)
  }
}

#########################################################################################
# function to generate pseudo replicates across regions 
split_rep1_rep2_meta = function(seurat_meta) {
  seurat_meta[['individual']] = sapply(row.names(seurat_meta), function(x) {
    strsplit(x, '_')[[1]][2]
  })
  cross = table(seurat_meta$individual, seurat_meta$region)
  
  # rep1 is the max
  rep1_inds = c()
  # rep 2 is the second highest
  rep2_inds = c()
  
  for (reg in colnames(cross)) {
    sorted = sort(cross[,reg], decreasing = T)
    if (length(which(sorted == 0)) > 4) {
      print(paste0('Region ', reg, ' has representation from fewer than 2 replicates.'))
    }
    rep1_inds = c(rep1_inds, names(sorted)[1])
    rep2_inds = c(rep2_inds, names(sorted)[2])
  }
  names(rep1_inds) = names(rep2_inds) = colnames(cross)
  
  seurat_meta[['rep_desig']] = 'NA'
  for (lobe in names(rep1_inds)) {
    rep1_ind_val = rep1_inds[lobe]
    seurat_meta[seurat_meta$individual == rep1_ind_val & seurat_meta$region == lobe, 'rep_desig'] = 'rep1'
  } 
  for (lobe in names(rep2_inds)) {
    rep2_ind_val = rep2_inds[lobe]
    seurat_meta[seurat_meta$individual == rep2_ind_val & seurat_meta$region == lobe, 'rep_desig'] = 'rep2'
  } 
  
  return(seurat_meta)
}


#########################################################################################
# Read in metadata corresponding to annotated cb object
cluster_metadata = read.csv('../data/cluster_metadata.csv', header = T, row.names = 1,
                            stringsAsFactors = F)

cluster_metadata$region = factor(cluster_metadata$region, levels = c("I","II", "III", "CUL", "VI", "VII", "VIII",  "IX","X","AN1", "AN2",
                                                                       "PRM", "SIM", "COP", "F",  "PF"))

cluster_metadata$individual = sapply(row.names(cluster_metadata), function(x) {
  strsplit(x, '_')[[1]][2]
})

# first two are "female" colors, last 4 are "male" colors 
colors.use = c('#CC0066', '#FF66CC', '#3366FF', '#000099', '#3399CC', '#3333FF')

# pdf('ExtendedDataFig1_panelA.pdf', useDingbats = F, width = 5, height = 3)
ggplot(cluster_metadata, aes(x = region, fill = individual)) + geom_bar() + 
  scale_fill_manual(values=colors.use) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.95))+
  xlab('Region') + ylab('Number of nuclei') + theme(legend.title = element_blank())
# dev.off()

# Reproducing panel B

cluster_metadata$log_nUMI = log10(cluster_metadata$nUMI)

# pdf('ExtendedDataFig1_panelB.pdf', useDingbats = F,
#     width = 4, height = 3)
ggplot(cluster_metadata, aes(x = region, y = log_nUMI))  + geom_violin() +  
  geom_boxplot(width = 0.25,  outlier.shape = NA) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.95)) +
  xlab('Region') + ylab('log10(nUM)')
# dev.off()

# Reproducing panel C

cluster_metadata$ordered_cluster = factor(cluster_metadata$cluster, 
                                          levels = c('Astrocyte', 'Bergmann', 'Choroid', 'Fibroblast',
                                                     'Endothelial_mural', 'Endothelial_stalk', 'Ependymal',
                                                     'Golgi', 'Granule', 'Macrophage', 'Microglia', 'MLI1',
                                                     'MLI2', 'ODC', 'PLI', 'OPC', 'Purkinje', 'UBC'))

# pdf('ExtendedDataFig1_panelC.pdf', useDingbats = F,
#     width = 4, height = 3)
ggplot(cluster_metadata, aes(x = ordered_cluster, y = log_nUMI))  + geom_violin() +  
  geom_boxplot(width = 0.25, outlier.shape = NA) +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.95)) +
  ylab('log10(nUMI)')
# dev.off()

# Reproducing panel D
# Requires annotated Seurat object with downsampled granules
full_cb = readRDS('../data/cb_annotated_object.RDS')
cluster_metadata_ds = read.csv('../data/cluster_metadata_ds.csv', header = T, row.names = 1,
                               stringsAsFactors = F)
cb_ds = SubsetData(full_cb, cells.use = row.names(cluster_metadata_ds), subset.raw = T)
cb_ds@meta.data$region = factor(cb_ds@meta.data$region)

cb_ds@meta.data$individual = sapply(row.names(cb_ds@meta.data), function(x) {
  strsplit(x, '_')[[1]][2]
})

# Computing alignments for each regional subset after applying standard liger pipeline
# this takes several hours to run
region_alignments = c()
for (region in levels(cb_ds@meta.data$region)) {
  reg_s = SubsetData(cb_ds, subset.name = 'region', accept.value = region, subset.raw = T)

  print(dim(reg_s@meta.data))
  print(table(reg_s@meta.data$individual))

  if (length(unique(reg_s@meta.data$sex)) == 1) {
    print(paste0('Skipping region ', region, ' as all same sex.'))
    region_alignments = c(region_alignments, NA)
  } else {
    reg_l = seuratToLiger(reg_s, combined.seurat = T, names = 'use-meta', meta.var = 'sex')
    reg_l = normalize(reg_l)
    reg_l = liger::selectGenes(reg_l, var.thresh = 0.2)
    print(length(reg_l@var.genes))
    reg_l = liger::scaleNotCenter(reg_l)
    reg_l = optimizeALS(reg_l, k = 30)
    reg_l@H.norm = do.call(rbind, reg_l@H)
    reg_l = runTSNE(reg_l)
    plotByDatasetAndCluster(reg_l)
    reg_s_back = ligerToSeurat2(reg_l)
    row.names(reg_s_back@meta.data) = row.names(reg_l@cell.data)
    reg_s_back@meta.data[['individual']] = reg_s@meta.data[row.names(reg_s_back@meta.data),]$individual
    
    # use at least 30 neighbors
    orig.min = ceiling(x = min(table(reg_s_back@meta.data$individual)) * 0.01 * length(table(reg_s_back@meta.data$individual)))
    align = CalcAlignmentMetric(reg_s_back, reduction.use = 'inmf', grouping.var = 'individual', dims.use = 1:30,
                                nn = max(orig.min, 30))
    region_alignments = c(region_alignments, align)
  }
}

names(region_alignments) = levels(cb_ds@meta.data$region)
region_align_df = data.frame(align = region_alignments,
                             region = factor(names(region_alignments), levels = c("I","II", "III", "CUL", "VI", "VII", "VIII",  
                                                                                  "IX","X","AN1", "AN2",
                                                                                  "PRM", "SIM", "COP", "F",  "PF")))
# pdf('ExtendedDataFig1_panelD.pdf', useDingbats = F, width = 4, height = 3)
ggplot(na.omit(region_align_df),  aes(region, align)) + geom_bar(stat = 'identity', na.rm = T) +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.95)) +
  xlab('Region') + ylab('Alignment')
# dev.off()

# Reproducing panel E
# Requires use of Seurat subanalysis objects or t-SNE/UMAP coordinates for these objects

purk_s = readRDS('purkinje_seurat.RDS')
purk_s@meta.data = split_rep1_rep2_meta(purk_s@meta.data)
purk_s = SetIdent(purk_s, ident.use = purk_s@meta.data$rep_desig)

purk_reps = SubsetData(purk_s, ident.use = c('rep1', 'rep2'), subset.raw = T)
purk_align = CalcAlignmentMetric(purk_reps, reduction.use = 'inmf', dims.use = c(1:20), grouping.var = 'ident')
purk_reps_plot = DimPlot2(purk_reps, reduction.use = 'tsne', shuffle_points = T, no.axes = T, 
                          no.legend = T, pt.size = 0.05, do.return = T)
purk_reps_plot = purk_reps_plot + ggtitle(paste0('Purkinje, a=', round(purk_align, 2)))

bergmann_s = readRDS('bergmann_seurat.RDS')
bergmann_s@meta.data = split_rep1_rep2_meta(bergmann_s@meta.data)
bergmann_s = SetIdent(bergmann_s, ident.use = bergmann_s@meta.data$rep_desig)
bergmann_reps = SubsetData(bergmann_s, ident.use = c('rep1', 'rep2'), subset.raw = T)
berg_align = CalcAlignmentMetric(bergmann_reps, reduction.use = 'inmf', dims.use = c(1:20), grouping.var = 'ident')

berg_reps_plot = DimPlot2(bergmann_reps, reduction.use = 'tsne', shuffle_points = T, no.axes = T, 
                          no.legend = T, pt.size = 0.05, do.return = T)
berg_reps_plot = berg_reps_plot + ggtitle(paste0('Bergmann, a=', round(berg_align, 2)))

ints_s = readRDS('interneurons_seurat.RDS')
ints_s@meta.data = split_rep1_rep2_meta(ints_s@meta.data)
ints_s = SetIdent(ints_s, ident.use = ints_s@meta.data$rep_desig)
ints_reps = SubsetData(ints_s, ident.use = c('rep1', 'rep2'), subset.raw = T)
ints_align = CalcAlignmentMetric(ints_reps, reduction.use = 'inmf', dims.use = c(1:35), grouping.var = 'ident')

ints_reps_plot = DimPlot2(ints_reps, reduction.use = 'tsne', shuffle_points = T, no.axes = T, 
                          no.legend = T, pt.size = 0.05, do.return = T)
ints_reps_plot = ints_reps_plot + ggtitle(paste0('MLI/PLI, a=', round(ints_align, 2)))

granule_s = readRDS('granule_seurat.RDS')
granule_s@meta.data = split_rep1_rep2_meta(granule_s@meta.data)
granule_s = SetIdent(granule_s, ident.use = granule_s@meta.data$rep_desig)
gran_reps = SubsetData(granule_s, ident.use = c('rep1', 'rep2'), subset.raw = T)
gran_align = CalcAlignmentMetric(gran_reps, reduction.use = 'inmf', dims.use = c(1:18), grouping.var = 'ident')

# note that this is actually UMAP but is labeled as tSNE
gran_reps_plot = DimPlot2(gran_reps, reduction.use = 'tsne', shuffle_points = T, no.axes = T, 
                          no.legend = T, pt.size = 0.05, do.return = T)
gran_reps_plot = gran_reps_plot + ggtitle(paste0('Granule, a=', round(gran_align, 2)))

# png("ExtendedDataFig1_panelE.pdf", units="in", width=16, height=4, res=300)
plot_grid(purk_reps_plot, gran_reps_plot, ints_reps_plot, berg_reps_plot, nrow = 1)
# dev.off()

# Reproducing panel F
# using binomial model suggested by the howmanycells tools of Satija lab

# Number of cells considered
x = seq(4.5e5, 6.7e5, length.out = 1000)
rare_type_prop = 0.00015
min_cells_detection = 70
num_rare_pops = 10
atlas_count = nrow(cluster_metadata)
y = (pnbinom(x - min_cells_detection, min_cells_detection, rare_type_prop))^num_rare_pops

curr_prob = (pnbinom(atlas_count - min_cells_detection, min_cells_detection, rare_type_prop))^num_rare_pops
annotation = data.frame(
  x = c(600000),
  y = c(curr_prob + 0.04),
  label = c(round(curr_prob, 3))
)

estimated = data.frame(num_cells = x, 
                       prob_detection = y)

# pdf('ExtendedDataFig1_panelF.pdf', useDingbats = F,
#     width = 5, height = 3)
ggplot(estimated, aes(num_cells, prob_detection)) + geom_line(col = 'blue') + 
  geom_vline(xintercept = atlas_count, linetype = 'dashed') + geom_text(data=annotation,
                                                                        aes(x=x, y=y, label=label)) +
  xlab('Number of cells sampled') + ylab('Probability')
# dev.off()





