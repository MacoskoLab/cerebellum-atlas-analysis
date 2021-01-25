# Figure 3 

# read in full_cb object as in previous panels 
# compute logistic fits for all cell type combinations shown

mli_comb = SubsetData(full_cb, subset.name = 'cluster', accept.value = c('MLI1', 'MLI2'))

# prepare basket_stellate results with new function 
mli_comb = SetIdent(mli_comb, ident.use = mli_comb@meta.data$cluster)
mli_returned = generateObjectContinuityFits(mli_comb, 'MLI1', use.raw.order = F, 
                                            'MLI2', reverse_order = T)


# golgi-mli already downsampled object 
# can produce similar plot (different seed) by using downsampling functionality 
# golgi_mli = readRDS('golgi_mli2_s_new.RDS')
# golgi_mli = SetIdent(golgi_mli2, ident.use = golgi_mli2@meta.data$cell_type)

golgi_mli = SubsetData(full_cb, subset.name = 'cluster', accept.value = c('MLI1', 'Golgi'))

gmli_returned = generateObjectContinuityFits(golgi_mli, 'MLI1',  
                                             'Golgi', use.raw.order = F, reverse_order = T)


# use this to exactly reproduce figure
# mli1_ds = readRDS('mli1_s_ds5k.RDS')
# mli1_ds = SetIdent(mli1_ds, ident.use = mli1_ds@meta.data$named_liger_ident)

mli1_s = SubsetData(full_cb, subset.name = 'cluster', accept.value = c('MLI1'))
mli1_s = SetIdent(mli1_s, ident.use = mli1_s@meta.data$subcluster)

mli1_returned = generateObjectContinuityFits(mli1_ds, 'MLI1_1',  
                                             'MLI1_2', use.raw.order = F, reverse_order = T)

golgi_s = SubsetData(full_cb, subset.name = 'cluster', accept.value = c('Golgi'))
golgi_s = SetIdent(golgi_s, ident.use = golgi_s@meta.data$subcluster)

golgi_returned = generateObjectContinuityFits(golgi_s, 'Golgi_1',  
                                              'Golgi_2', use.raw.order = F, reverse_order = F)

ubc_s = SubsetData(full_cb, subset.name = 'cluster', accept.value = c('UBC'))
ubc_s = SetIdent(ubc_s, ident.use = ubc_s@meta.data$subcluster)

ubc_returned = generateObjectContinuityFits(ubc_s, 'UBC_1',  
                                            'UBC_3', use.raw.order = F, reverse_order = F)


# panel a
grm8 = plotFittedFunc2(mli1_returned$object, 'Grm8', factor_order = mli1_returned[[2]], 
                        fitparams = mli1_returned[[1]], return.plot = T, use_scale = T, do.legend = F)
npas3 = plotFittedFunc2(mli1_returned$object, 'Npas3', factor_order = mli1_returned[[2]], 
                        fitparams = mli1_returned[[1]], return.plot = T, use_scale = T, do.legend = F)

ptprk = plotFittedFunc2(mli_returned$object, 'Ptprk', factor_order = mli_returned[[2]], 
                        fitparams = mli_returned[[1]], return.plot = T, use_scale = T, do.legend = F)
nxph1 = plotFittedFunc2(mli_returned$object, 'Nxph1', factor_order = mli_returned[[2]], 
                        fitparams = mli_returned[[1]], return.plot = T, use_scale = T, do.legend = F)


pdf('Figure3_panelA.pdf', useDingbats = F, width = 4, height = 9)
plot_grid(grm8, npas3, ptprk, nxph1, ncol = 1, align = 'v', axis = 'l')
dev.off()


# panel b
new_comparison_200 = makeCompareDF2(list(gmli_returned$fits,
                                         mli1_returned$fits, 
                                         mli_returned$fits,
                                         ubc_returned$fits,
                                         golgi_returned$fits),
                                    n_cell_counts = c(5000, 5000, 5000, 1613, 3989),
                                    names = c('MLI1_Golgi', 'MLI1_1_MLI1_2', 'MLI1_MLI2', 'UBC', 'Golgi'),
                                    top_ranked_spearman = 200, converged_only = T, 
                                    positive_v = T)

pdf('Figure3_panelB.pdf', useDingbats = F,  width = 7, height = 4)
ggplot(top_200_fits, aes(x = log_midslope, col = dataset)) + stat_ecdf()
dev.off()



# panel d
grm1 = plotFittedFunc2(ubc_returned$object, 'Grm1', factor_order = ubc_returned$pseudo_order, color_by_cluster = F,
                       continuous_color = 'pseudotime', use_gradient = T,
                       fitparams = ubc_returned$fits, return.plot = T, use_scale = T, do.legend = F)
calb2 = plotFittedFunc2(ubc_returned$object, 'Calb2', factor_order = ubc_returned$pseudo_order, color_by_cluster = F,
                        continuous_color = 'pseudotime', use_gradient = T,
                        fitparams = ubc_returned$fits, return.plot = T, use_scale = T, do.legend = F)
plcb4 = plotFittedFunc2(ubc_returned$object, 'Plcb4', factor_order = ubc_returned$pseudo_order, color_by_cluster = F,
                        continuous_color = 'pseudotime', use_gradient = T,
                        fitparams = ubc_returned$fits, return.plot = T, use_scale = T, do.legend = F)
plcb1 = plotFittedFunc2(ubc_returned$object, 'Plcb1', factor_order = ubc_returned$pseudo_order, color_by_cluster = F,
                        continuous_color = 'pseudotime', use_gradient = T,
                        fitparams = ubc_returned$fits, return.plot = T, use_scale = T, do.legend = F)


pdf('Figure3_panelD.pdf', useDingbats = F, width = 4, height = 9)
plot_grid(grm1, calb2, plcb4, plcb1, ncol = 1, align = 'v', axis = 'l')
dev.off()

