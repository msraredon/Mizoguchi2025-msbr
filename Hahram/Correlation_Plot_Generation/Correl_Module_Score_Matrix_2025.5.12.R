ATI <- ATI_Module_Values
ATII <- ATII_Module_Values
ATII_ATI <- ATII_ATI_Module_Values
Basal <- Basal_Module_Values
BASC <- BASC_Module_Values
Ciliated <- Ciliated_Module_Values
Hillock <- Hillock_Module_Values
Secretory <- Secretory_Module_Values
Tuft <-Tuft_Module_Values
AbBasal <-Aberrant_Basaloid_Module_Values
 ref.numbers <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
 raw_module_values <- list(ATI, ATII, ATII_ATI, Basal, BASC, Ciliated, Hillock, Secretory, Tuft, AbBasal)
raw_module_values_2 <- list(ATI, ATII, ATII_ATI, Basal, BASC, Ciliated, Hillock, Secretory, Tuft, AbBasal)
ref.names <- c("ATI", "ATII", "ATII_ATI", "Basal", "BASC", "Ciliated", "Hillock", "Secretory", "Tuft", "AbBasal")                   
 for (i in ref.numbers) {
for (j in ref.numbers){
# if (i == j){ 
# j = j+1
# return
#  }
# else{
x = unlist(raw_module_values[i])
y =  unlist(raw_module_values_2[j])
 result <- cor(x, y, method = "spearman")
assign(paste0(ref.names[i], "plus", ref.names[j], "SpearmanValue"), result)
 }
 }
# }
save.image("/gpfs/gibbs/project/kaminski/hk738/Storage/Satoshi_Data/BASC_Analysis/BASC_Combined/2025.5.9/Spearman_Correlation_Module_Score_Data.2025.5.9.RData")

####matrix shape will be this:
#          AbBasal ATI ATII ATII_ATI Basal BASC Ciliated Hillock Secretory Tuft
# AbBasal
# 
# ATI
# 
# ATII
# 
# ATII_ATI
# 
# Basal
# 
# BASC
# 
# Ciliated
# 
# Hillock
# 
# Secretory
# 
# Tuft
matrix.row.1 <- c(AbBasalplusAbBasalSpearmanValue, AbBasalplusATISpearmanValue, AbBasalplusATIISpearmanValue, AbBasalplusATII_ATISpearmanValue,AbBasalplusBasalSpearmanValue, AbBasalplusBASCSpearmanValue, AbBasalplusCiliatedSpearmanValue, AbBasalplusHillockSpearmanValue, AbBasalplusSecretorySpearmanValue, AbBasalplusTuftSpearmanValue)

matrix.row.2 <- c(ATIplusAbBasalSpearmanValue, ATIplusATISpearmanValue, ATIplusATIISpearmanValue, ATIplusATII_ATISpearmanValue, ATIplusBasalSpearmanValue, ATIplusBASCSpearmanValue, ATIplusCiliatedSpearmanValue, ATIplusHillockSpearmanValue, ATIplusSecretorySpearmanValue, ATIplusTuftSpearmanValue)

matrix.row.3 <- c(ATIIplusAbBasalSpearmanValue, ATIIplusATISpearmanValue, ATIIplusATIISpearmanValue, ATIIplusATII_ATISpearmanValue, ATIIplusBasalSpearmanValue, ATIIplusBASCSpearmanValue, ATIIplusCiliatedSpearmanValue, ATIIplusHillockSpearmanValue, ATIIplusSecretorySpearmanValue, ATIIplusTuftSpearmanValue)

matrix.row.4 <- c(ATII_ATIplusAbBasalSpearmanValue, ATII_ATIplusATISpearmanValue, ATII_ATIplusATIISpearmanValue, ATII_ATIplusATII_ATISpearmanValue, ATII_ATIplusBasalSpearmanValue, ATII_ATIplusBASCSpearmanValue, ATII_ATIplusCiliatedSpearmanValue, ATII_ATIplusHillockSpearmanValue, ATII_ATIplusSecretorySpearmanValue, ATII_ATIplusTuftSpearmanValue)

matrix.row.5 <- c(BasalplusAbBasalSpearmanValue, BasalplusATISpearmanValue, BasalplusATIISpearmanValue, BasalplusATII_ATISpearmanValue, BasalplusBasalSpearmanValue, BasalplusBASCSpearmanValue, BasalplusCiliatedSpearmanValue, BasalplusHillockSpearmanValue, BasalplusSecretorySpearmanValue, BasalplusTuftSpearmanValue)

matrix.row.6 <- c(BASCplusAbBasalSpearmanValue, BASCplusATISpearmanValue, BASCplusATIISpearmanValue, BASCplusATII_ATISpearmanValue, BASCplusBasalSpearmanValue, BASCplusBASCSpearmanValue, BASCplusCiliatedSpearmanValue, BASCplusHillockSpearmanValue, BASCplusSecretorySpearmanValue, BASCplusTuftSpearmanValue)

matrix.row.7 <- c(CiliatedplusAbBasalSpearmanValue, CiliatedplusATISpearmanValue, CiliatedplusATIISpearmanValue, CiliatedplusATII_ATISpearmanValue, CiliatedplusBasalSpearmanValue, CiliatedplusBASCSpearmanValue, CiliatedplusCiliatedSpearmanValue, CiliatedplusHillockSpearmanValue, CiliatedplusSecretorySpearmanValue, CiliatedplusTuftSpearmanValue)

matrix.row.8 <- c(HillockplusAbBasalSpearmanValue, HillockplusATISpearmanValue, HillockplusATIISpearmanValue, HillockplusATII_ATISpearmanValue, HillockplusBasalSpearmanValue, HillockplusBASCSpearmanValue, HillockplusCiliatedSpearmanValue, HillockplusHillockSpearmanValue, HillockplusSecretorySpearmanValue, HillockplusTuftSpearmanValue)

matrix.row.9 <- c(SecretoryplusAbBasalSpearmanValue, SecretoryplusATISpearmanValue, SecretoryplusATIISpearmanValue, SecretoryplusATII_ATISpearmanValue, SecretoryplusBasalSpearmanValue, SecretoryplusBASCSpearmanValue, SecretoryplusCiliatedSpearmanValue, SecretoryplusHillockSpearmanValue, SecretoryplusSecretorySpearmanValue, SecretoryplusTuftSpearmanValue)

matrix.row.10 <- c(TuftplusAbBasalSpearmanValue, TuftplusATISpearmanValue, TuftplusATIISpearmanValue, TuftplusATII_ATISpearmanValue, TuftplusBasalSpearmanValue, TuftplusBASCSpearmanValue, TuftplusCiliatedSpearmanValue, TuftplusHillockSpearmanValue, TuftplusSecretorySpearmanValue, TuftplusTuftSpearmanValue)


module.correl.matrix = matrix(
  
  # Taking sequence of elements  
  c(matrix.row.1, matrix.row.2, matrix.row.3, matrix.row.4, matrix.row.5, matrix.row.6, matrix.row.7, matrix.row.8, matrix.row.9, matrix.row.10), 
  
  # No of rows
  nrow = 10,   
  
  # No of columns
  ncol = 10,         
  
  # By default matrices are in column-wise order
  # So this parameter decides how to arrange the matrix
  byrow = TRUE          
)

rownames(module.correl.matrix) = c("AbBasal", "ATI", "ATII", "ATII_ATI", "Basal", "BASC", "Ciliated", "Hillock", "Secretory", "Tuft")

# Naming columns
colnames(module.correl.matrix) = c("AbBasal", "ATI", "ATII", "ATII_ATI", "Basal", "BASC", "Ciliated", "Hillock", "Secretory", "Tuft")

cat("The 3x3 matrix:\n")
print(module.correl.matrix)


install.packages("corrplot")
library(corrplot)

M<-cor(mtcars)

setwd("/gpfs/gibbs/project/kaminski/hk738/Storage/Satoshi_Data/BASC_Analysis/BASC_Combined/2025.5.9")
png(file = 'Module_Score_Correlation_Matrix.png', width=8.5, height=8.5, units = 'in',res=300, type = "cairo")
corrplot(module.correl.matrix, method="circle")
dev.off()


module_correl_df <- as.data.frame(module.correl.matrix)
module_correl_df<- signif(module_correl_df, 3)

gridExtra::grid.table(module_correl_df)
myTable <- tableGrob(module_correl_df, rows = NULL)

myTable <- tableGrob(
  module_correl_df, 
  theme = ttheme_default(core = list(bg_params = list(fill = "grey99")), base_size = 10)
  
)

png(file = "Correl_Validation.png", width=9, height=5, units = 'in',res=300, type = "cairo")
grid.draw(myTable)
dev.off()