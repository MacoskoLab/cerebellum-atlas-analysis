# Figure 4

# read in full_cb object as in previous figures
ints_s = SubsetData(full_cb, subset.name = 'cluster', accept.value = c('MLI1', 'MLI2', 'PLI'))
ints_s = NormalizeData(ints_s)

# panel a
genes.plot = c("Htr2a","Slc6a5", "Tns1","Slc16a12","Acvr1c", "Npas3",
               "Grm8",'Cdh22', 'Nxph1', 'Prkcd', 'Sorcs3',  'Ptprk')


pdf("Figure4_panelA.pdf",useDingbats = F,height = 3.5, width = 6)
DotPlot2(ints_s, genes.plot = genes.plot,dot.scale = 10, x.lab.rot = 45, plot.legend = T)
dev.off()

# panel c
library(dplyr)
# reshape2 also required

# read in HCR data for SORCS3, NXPH1, GRM8
grm_data = read.csv('sorcs_nxph_grm_data.csv', stringsAsFactors = F)
grm_melted = reshape2::melt(grm_data, variable.name = 'slide', value.name = 'count')

summed2 =  grm_melted %>% group_by(slide, position, mli_ident) %>% summarise(sum(count)) %>% ungroup()
colnames(summed2) = c('slide', 'position', 'mli_ident', 'sum')

summed2$position = factor(summed2$position, levels = c('distal', 'middle', 'inner'))
merged2 = summed2 %>% group_by(position, slide) %>% summarise(sum(sum)) %>%
  rename('total' = 'sum(sum)') %>% right_join(summed2)
merged2[['percent_pos']] = merged2$sum / merged2$total

merged2$mli_ident = factor(merged2$mli_ident, levels = rev(c('sorcs', 'nxph', 'nxph_sorcs')))

pdf('Figure4_panelC.pdf', useDingbats = F, width = 6, height = 4)
# set seed for jittering reproducibility
set.seed(1)
merged2 %>% group_by(position, mli_ident) %>% summarise(mean(percent_pos), sd(percent_pos)) %>%
  ggplot(aes(x = position, y = `mean(percent_pos)`, fill = mli_ident)) +
  geom_bar(stat = 'identity',
           position=position_dodge(), width = 0.98) + 
  geom_errorbar(aes(x=position, ymin=`mean(percent_pos)` + `sd(percent_pos)`, 
                    ymax= `mean(percent_pos)` - `sd(percent_pos)`),
                position=position_dodge(0.95), width = 0.3, size = 0.3) + 
  geom_point(data = merged2, aes(x = position, y = percent_pos,
                                 group = interaction(factor(position),factor(mli_ident))),
             position = position_jitterdodge(dodge.width=0.95), alpha = 0.6, shape = 1) + 
  scale_shape(solid = FALSE) + scale_y_continuous(limits = c(0, 0.84), expand = c(0,0)) + 
  theme(axis.text.y = element_text(angle = 90)) + 
  coord_flip()
dev.off()


# for panel d see separate file 


# panel i
# convert to liger object

cluster_plots = plotByDatasetAndCluster(ints_l, return.plots = T, do.legend = F, pt.size = 0.005, text.size = 0)
cluster_plots[[2]] = cluster_plots[[2]] + theme(
  axis.line = element_blank(), axis.text.x = element_blank(),
  axis.text.y = element_blank(), axis.ticks = element_blank(),
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  panel.background = element_blank(), panel.border = element_blank(),
  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  plot.background = element_blank(), plot.title = element_blank()) +
  scale_color_manual(values = c('#C13EB2', '#3EC14D', 'blue', 'blue', 'blue'))

# ordered so positive expression in front
gjd2 = plotGene_2020(ints_l, 'Gjd2', plot.by = 'none', points.only = T, do.legend = F, 
                     max.clip = 5, clip.absolute = T, pt.size = 0.001, return.plots = T)


# plot as png
png("Figure4_panelI.png", units="in", width=6, height=8, res=300)
plot_grid(cluster_plots[[2]], gjd2, nrow = 2)
dev.off()

# panel k

# read in HCR data for SORCS3, NXPH1, GJD2
gjd_data = read.csv('sorcs_nxph_gjd_data.csv', stringsAsFactors = F)
gjd_melted = reshape2::melt(gjd_data, variable.name = 'slide', value.name = 'count')

# summarise on two levels
summed =  gjd_melted %>% group_by(slide, mli_ident) %>% summarise(sum(count)) %>% ungroup()
summed_gjd = gjd_melted %>% group_by(slide, mli_ident, gjd_pos) %>% summarise(sum(count)) %>% 
  rename('sum_gjd' = 'sum(count)') %>% ungroup()

merged = inner_join(summed, summed_gjd[summed_gjd$gjd_pos == 'TRUE',])
merged[['percent']] = merged$sum_gjd / merged[['sum(count)']]
colnames(merged) =  c('slide', 'mli_ident', 'total', 'gjd_pos', 'sum_gjd', 'percent')

pdf('Figure4_panelK.pdf', useDingbats = F, width = 5, height = 2.5)
set.seed(1)
merged %>% group_by(mli_ident) %>% summarise(mean(percent), sd(percent)) %>%
  ggplot(aes(x = mli_ident, y = `mean(percent)`)) +
  geom_bar(stat = 'identity',
           position=position_dodge(width = 0.98)) + 
  geom_errorbar(aes(x=mli_ident, ymin=`mean(percent)` + `sd(percent)`, ymax= `mean(percent)` - `sd(percent)`),
                position=position_dodge(0.95), width = 0.3, size = 0.3) + 
  coord_flip() +
  geom_point(data = merged, aes(x = mli_ident, y = percent),
             position = position_jitter(width = 0.3), alpha = 0.6, shape = 1, size = 2) + 
  scale_shape(solid = FALSE) + scale_y_continuous(limits = c(0, 1.03), expand = c(0,0)) + 
  theme(axis.text.y = element_text(angle = 90)) 
dev.off()

