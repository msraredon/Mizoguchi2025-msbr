load("/gpfs/gibbs/project/kaminski/hk738/Storage/Satoshi_Data/BASC_Analysis/BASC_Combined/2025.5.9/Spearman_Correlation_Module_Score_data.2025.5.12.RData")

library(seriation)
list_seriation_methods('matrix')

module_correl_df

list_seriation_methods('dist')

data(Zoo)


dist2order = function(corr, method, ...) {
  d_corr = as.dist(1 - corr)
  s = seriate(d_corr, method = method, ...)
  i = get_order(s)
  return(i)
}

i = dist2order(module.correl.matrix, 'OLO')
corrplot(module.correl.matrix[i, i])

module.correl.matrix
i = dist2order(module.correl.matrix, 'OLO_average')
corrplot(module.correl.matrix[i, i])


png(file = 'Module_Score_Correlation_Matrix.OLO.png', width=8.5, height=8.5, units = 'in',res=300, type = "cairo")
i = dist2order(module.correl.matrix, 'OLO')
corrplot(module.correl.matrix[i, i], method = "circle")
dev.off()

png(file = 'Module_Score_Correlation_Matrix.OLO.avg.png', width=8.5, height=8.5, units = 'in',res=300, type = "cairo")
i = dist2order(module.correl.matrix, 'OLO_average')
corrplot(module.correl.matrix[i, i], method = "circle")
dev.off()


