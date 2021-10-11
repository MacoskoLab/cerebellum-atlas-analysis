# Reproducing Extended Data Figure 9 panels
library(liger)
# v0.4.2
library(ggplot2)

# Panel A

# Read in golgi as liger object 
golgi_l = readRDS('../data/golgi_liger.RDS')
plot_list = list()
for (gene in c('Lgi2', 'Slc6a5', 'Gad2', 'Grm2', 'Sorcs3', 'Gjd2', 'Ptprt', 'Sst', 'Nxph1', 'Ptprk')) {
  plot_list[[gene]] = plotGene_2020(golgi_l, gene, plot.by = 'none', return.plots = T, points.only = T)
}

# png('ExtendedDataFig9_panelA.png', units="in", width=10, height=4, res=300)
plot_grid(plotlist = plot_list, nrow = 2)
# dev.off()

# Panel C
# Read in intensity data for golgi cells 

golgi_sst_gjd = read.csv('../data/SST_GJD_LGI_ROI_intensities.csv')

corr_val = cor.test(golgi_sst_gjd$Mean_56_..GJD2., golgi_sst_gjd$Mean_647_.SST., 
                    method = 'spearman')
corr_p_val = corr_val$p.value
corr_stat = corr_val$estimate

# pdf('ExtendedDataFig9_panelC.pdf', useDingbats = F, width = 5, height = 3)
ggplot(golgi_sst_gjd, aes(x = Mean_647_.SST., y = Mean_56_..GJD2.)) + geom_point() + 
  xlab('Sst') + ylab('Gjd2') + 
  geom_text(x=5000, y=3000, label=paste0("rho = ", round(corr_stat,2)), size = 5) + theme_classic() +
  xlab('Mean pixel intensity, Sst') + ylab('Mean pixel intensity, Gjd2')
# dev.off()


