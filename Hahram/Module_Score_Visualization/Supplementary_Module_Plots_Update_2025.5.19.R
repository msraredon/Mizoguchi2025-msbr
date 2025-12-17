load("/gpfs/gibbs/project/kaminski/hk738/Storage/Satoshi_Data/BASC_Analysis/BASC_Combined/2025.5.5/HK_tissue_engineering_epithelial_merged.dataset.2025.5.5.Robj")

####Redoing the Text Sizes

ATI_Module_Values <- HK_tissue_engineering_epithelial_merged.dataset.2025.5.5@meta.data$ATI_Module1
ATII_Module_Values <- HK_tissue_engineering_epithelial_merged.dataset.2025.5.5@meta.data$ATII_Module1
ATII_ATI_Module_Values <- HK_tissue_engineering_epithelial_merged.dataset.2025.5.5@meta.data$ATII_ATI_Module1
Basal_Module_Values <- HK_tissue_engineering_epithelial_merged.dataset.2025.5.5@meta.data$Basal_Module1
BASC_Module_Values <- HK_tissue_engineering_epithelial_merged.dataset.2025.5.5@meta.data$BASC_Module1
Ciliated_Module_Values <- HK_tissue_engineering_epithelial_merged.dataset.2025.5.5@meta.data$Ciliated_Module1
Hillock_Module_Values <- HK_tissue_engineering_epithelial_merged.dataset.2025.5.5@meta.data$Hillock_Module1
Secretory_Module_Values <- HK_tissue_engineering_epithelial_merged.dataset.2025.5.5@meta.data$Secretory_Module1
Tuft_Module_Values <- HK_tissue_engineering_epithelial_merged.dataset.2025.5.5@meta.data$Tuft_Module1
Aberrant_Basaloid_Module_Values <- HK_tissue_engineering_epithelial_merged.dataset.2025.5.5@meta.data$Aberrant_Basaloid_Module1

Condition_Names <- HK_tissue_engineering_epithelial_merged.dataset.2025.5.5@meta.data$ConditionRearrangedNewName

module_trend_plot <- tibble(ATI_Module_Values, ATII_Module_Values, 
                            ATII_ATI_Module_Values, Basal_Module_Values, 
                            BASC_Module_Values, Ciliated_Module_Values,
                            Hillock_Module_Values, Secretory_Module_Values, 
                            Tuft_Module_Values, Aberrant_Basaloid_Module_Values, Condition_Names)
module_trend_plot.2 <- pivot_longer(module_trend_plot, 
                                    cols = c("ATI_Module_Values", 
                                             "ATII_ATI_Module_Values", 
                                             "ATII_Module_Values",
                                             "Basal_Module_Values", 
                                             "BASC_Module_Values",
                                             "Ciliated_Module_Values",
                                             "Tuft_Module_Values",
                                             "Secretory_Module_Values", 
                                             "Aberrant_Basaloid_Module_Values",
                                             "Hillock_Module_Values"), 
                                    
                                    names_to = "Module_Type", 
                                    
                                    values_to = "Module_Score_Value"
)

module_trend_df <- module_trend_plot.2

module_trend_point_df.ATI <- filter(module_trend_df, Module_Type == "ATI_Module_Values")
module_trend_point_df.ATII <- filter(module_trend_df, Module_Type == "ATII_Module_Values")
module_trend_point_df.ATII_ATI <- filter(module_trend_df, Module_Type == "ATII_ATI_Module_Values")
module_trend_point_df.Basal <- filter(module_trend_df, Module_Type == "Basal_Module_Values")
module_trend_point_df.BASC <- filter(module_trend_df, Module_Type == "BASC_Module_Values")
module_trend_point_df.Ciliated <- filter(module_trend_df, Module_Type == "Ciliated_Module_Values")
module_trend_point_df.Hillock <- filter(module_trend_df, Module_Type == "Hillock_Module_Values")
module_trend_point_df.Secretory <- filter(module_trend_df, Module_Type == "Secretory_Module_Values")
module_trend_point_df.Tuft <- filter(module_trend_df, Module_Type == "Tuft_Module_Values")
module_trend_point_df.AbBasal <- filter(module_trend_df, Module_Type == "Aberrant_Basaloid_Module_Values")



module_trend_plot


######

setwd("/gpfs/gibbs/project/kaminski/hk738/Storage/Satoshi_Data/BASC_Analysis/BASC_Combined/2025.5.19/Box_Plot")
#####Box Plot
######
png(file = "AbBasal_Module_Score_Box_Plot.png", width=18, height=12, units = 'in',res=300, type = "cairo")
p <- ggplot(module_trend_point_df.AbBasal, aes(Condition_Names, Module_Score_Value)) +
  geom_boxplot(fill="#FF6666", alpha=1)
p + theme(axis.text=element_text(size=15, face = "bold"), axis.title=element_text(size=30,face="bold")) + theme(strip.text.x = element_text(size = 20, face="bold")) + labs(y = "Module Score Value", x = "Condition")
dev.off()

png(file = "ATI_Module_Score_Box_Plot.png", width=18, height=12, units = 'in',res=300, type = "cairo")
p <- ggplot(module_trend_point_df.ATI, aes(Condition_Names, Module_Score_Value)) +
  geom_boxplot(fill="#987284", alpha= 1)
p + theme(axis.text=element_text(size=15, face = "bold"), axis.title=element_text(size=30,face="bold")) + theme(strip.text.x = element_text(size = 20, face="bold")) + labs(y = "Module Score Value", x = "Condition")
dev.off()

png(file = "ATII_Module_Score_Box_Plot.png", width=18, height=12, units = 'in',res=300, type = "cairo")
p <- ggplot(module_trend_point_df.ATII, aes(Condition_Names, Module_Score_Value)) +
  geom_boxplot(fill="#156064", alpha=1)
p + theme(axis.text=element_text(size=15, face = "bold"), axis.title=element_text(size=30,face="bold")) + theme(strip.text.x = element_text(size = 20, face="bold")) + labs(y = "Module Score Value", x = "Condition")
dev.off()

png(file = "ATII_ATI_Module_Score_Box_Plot.png", width=18, height=12, units = 'in',res=300, type = "cairo")
p <- ggplot(module_trend_point_df.ATII_ATI, aes(Condition_Names, Module_Score_Value)) +
  geom_boxplot(fill="pink", alpha=1)
p + theme(axis.text=element_text(size=15, face = "bold"), axis.title=element_text(size=30,face="bold")) + theme(strip.text.x = element_text(size = 20, face="bold")) + labs(y = "Module Score Value", x = "Condition")
dev.off()

png(file = "Basal_Module_Score_Box_Plot.png", width=18, height=12, units = 'in',res=300, type = "cairo")
p <- ggplot(module_trend_point_df.Basal, aes(Condition_Names, Module_Score_Value)) +
  geom_boxplot(fill="#38AECC", alpha=1)
p  + theme(axis.text=element_text(size=15, face = "bold"), axis.title=element_text(size=30,face="bold")) + theme(strip.text.x = element_text(size = 20, face="bold")) + labs(y = "Module Score Value", x = "Condition")
dev.off()

png(file = "BASC_Module_Score_Box_Plot.png", width=18, height=12, units = 'in',res=300, type = "cairo")
p <- ggplot(module_trend_point_df.BASC, aes(Condition_Names, Module_Score_Value)) +
  geom_boxplot(fill="#D741A7", alpha=1)
p + theme(axis.text=element_text(size=15, face = "bold"), axis.title=element_text(size=30,face="bold")) + theme(strip.text.x = element_text(size = 20, face="bold")) + labs(y = "Module Score Value", x = "Condition")
dev.off()


png(file = "Ciliated_Module_Score_Box_Plot.png", width=18, height=12, units = 'in',res=300, type = "cairo")
p <- ggplot(module_trend_point_df.Ciliated, aes(Condition_Names, Module_Score_Value)) +
  geom_boxplot(fill="#A3E51F", alpha=1)
p + theme(axis.text=element_text(size=15, face = "bold"), axis.title=element_text(size=30,face="bold")) + theme(strip.text.x = element_text(size = 20, face="bold")) + labs(y = "Module Score Value", x = "Condition")
dev.off()

png(file = "Hillock_Module_Score_Box_Plot.png", width=18, height=12, units = 'in',res=300, type = "cairo")
p <- ggplot(module_trend_point_df.Hillock, aes(Condition_Names, Module_Score_Value)) +
  geom_boxplot(fill="#F5B700", alpha=1)
p + theme(axis.text=element_text(size=15, face = "bold"), axis.title=element_text(size=30,face="bold")) + theme(strip.text.x = element_text(size = 20, face="bold")) + labs(y = "Module Score Value", x = "Condition")
dev.off()

png(file = "Tuft_Module_Score_Box_Plot.png", width=18, height=12, units = 'in',res=300, type = "cairo")
p <- ggplot(module_trend_point_df.Tuft, aes(Condition_Names, Module_Score_Value)) +
  geom_boxplot(fill="#6DECAF", alpha=1)
p + theme(axis.text=element_text(size=15, face = "bold"), axis.title=element_text(size=30,face="bold")) + theme(strip.text.x = element_text(size = 20, face="bold")) + labs(y = "Module Score Value", x = "Condition")
dev.off()

png(file = "Secretory_Module_Score_Box_Plot.png", width=18, height=12, units = 'in',res=300, type = "cairo")
p <- ggplot(module_trend_point_df.Secretory, aes(Condition_Names, Module_Score_Value)) +
  geom_boxplot(fill="yellow2", alpha=1)
p + theme(axis.text=element_text(size=15, face = "bold"), axis.title=element_text(size=30,face="bold")) + theme(strip.text.x = element_text(size = 20, face="bold")) + labs(y = "Module Score Value", x = "Condition")
dev.off()

######
setwd("/gpfs/gibbs/project/kaminski/hk738/Storage/Satoshi_Data/BASC_Analysis/BASC_Combined/2025.5.19/Density_Plot")
#######Density Plot
######
condition_cols <- c("red", "blue", "yellow", "green", "purple", "orange")
png(file = "AbBasal_Module_Score_Density_Plot.png", width=18, height=12, units = 'in',res=300, type = "cairo")
p <- ggplot(module_trend_point_df.AbBasal, aes(Module_Score_Value, fill = Condition_Names, colors = Condition_Names)) +
  geom_density(alpha = 0.5) + scale_fill_manual(values = condition_cols, name = "Condition") + labs(y = "Density", x = "Module Score Value")  + theme(axis.text=element_text(size=25, face = "bold"), axis.title=element_text(size=30,face="bold")) + theme(strip.text.x = element_text(size = 20, face="bold")) 

p +theme(legend.text = element_text(size=25), legend.title = element_text(size = 30))
dev.off()

png(file = "ATI_Module_Score_Density_Plot.png", width=18, height=12, units = 'in',res=300, type = "cairo")
p <- ggplot(module_trend_point_df.ATI, aes(Module_Score_Value, fill = Condition_Names, colors = Condition_Names)) +
  geom_density(alpha = 0.5) + scale_fill_manual(values = condition_cols, name = "Condition") + labs(y = "Density", x = "Module Score Value")  + theme(axis.text=element_text(size=25, face = "bold"), axis.title=element_text(size=30,face="bold")) + theme(strip.text.x = element_text(size = 20, face="bold")) 

p +theme(legend.text = element_text(size=25), legend.title = element_text(size = 30))
dev.off()

png(file = "ATII_Module_Score_Density_Plot.png", width=18, height=12, units = 'in',res=300, type = "cairo")
p <- ggplot(module_trend_point_df.ATII, aes(Module_Score_Value, fill = Condition_Names, colors = Condition_Names)) +
  geom_density(alpha = 0.5) + scale_fill_manual(values = condition_cols, name = "Condition") + labs(y = "Density", x = "Module Score Value")  + theme(axis.text=element_text(size=25, face = "bold"), axis.title=element_text(size=30,face="bold")) + theme(strip.text.x = element_text(size = 20, face="bold")) 

p +theme(legend.text = element_text(size=25), legend.title = element_text(size = 30))
dev.off()

png(file = "ATII_ATI_Module_Score_Density_Plot.png", width=18, height=12, units = 'in',res=300, type = "cairo")
p <- ggplot(module_trend_point_df.ATII_ATI, aes(Module_Score_Value, fill = Condition_Names, colors = Condition_Names)) +
  geom_density(alpha = 0.5) + scale_fill_manual(values = condition_cols, name = "Condition") + labs(y = "Density", x = "Module Score Value")  + theme(axis.text=element_text(size=25, face = "bold"), axis.title=element_text(size=30,face="bold")) + theme(strip.text.x = element_text(size = 20, face="bold")) 

p +theme(legend.text = element_text(size=25), legend.title = element_text(size = 30))
dev.off()

png(file = "Basal_Module_Score_Density_Plot.png", width=18, height=12, units = 'in',res=300, type = "cairo")
p <- ggplot(module_trend_point_df.Basal, aes(Module_Score_Value, fill = Condition_Names, colors = Condition_Names)) +
  geom_density(alpha = 0.5) + scale_fill_manual(values = condition_cols, name = "Condition") + labs(y = "Density", x = "Module Score Value")  + theme(axis.text=element_text(size=25, face = "bold"), axis.title=element_text(size=30,face="bold")) + theme(strip.text.x = element_text(size = 20, face="bold")) 

p +theme(legend.text = element_text(size=25), legend.title = element_text(size = 30))
dev.off()

png(file = "BASC_Module_Score_Density_Plot.png", width=18, height=12, units = 'in',res=300, type = "cairo")
p <- ggplot(module_trend_point_df.BASC, aes(Module_Score_Value, fill = Condition_Names, colors = Condition_Names)) +
  geom_density(alpha = 0.5) + scale_fill_manual(values = condition_cols, name = "Condition") + labs(y = "Density", x = "Module Score Value")  + theme(axis.text=element_text(size=25, face = "bold"), axis.title=element_text(size=30,face="bold")) + theme(strip.text.x = element_text(size = 20, face="bold")) 

p +theme(legend.text = element_text(size=25), legend.title = element_text(size = 30))
dev.off()


png(file = "Ciliated_Module_Score_Density_Plot.png", width=18, height=12, units = 'in',res=300, type = "cairo")
p <- ggplot(module_trend_point_df.Ciliated, aes(Module_Score_Value, fill = Condition_Names, colors = Condition_Names)) +
  geom_density(alpha = 0.5) + scale_fill_manual(values = condition_cols, name = "Condition") + labs(y = "Density", x = "Module Score Value")  + theme(axis.text=element_text(size=25, face = "bold"), axis.title=element_text(size=30,face="bold")) + theme(strip.text.x = element_text(size = 20, face="bold")) 

p +theme(legend.text = element_text(size=25), legend.title = element_text(size = 30))
dev.off()

png(file = "Hillock_Module_Score_Density_Plot.png", width=18, height=12, units = 'in',res=300, type = "cairo")
p <- ggplot(module_trend_point_df.Hillock, aes(Module_Score_Value, fill = Condition_Names, colors = Condition_Names)) +
  geom_density(alpha = 0.5) + scale_fill_manual(values = condition_cols, name = "Condition") + labs(y = "Density", x = "Module Score Value")  + theme(axis.text=element_text(size=25, face = "bold"), axis.title=element_text(size=30,face="bold")) + theme(strip.text.x = element_text(size = 20, face="bold")) 

p +theme(legend.text = element_text(size=25), legend.title = element_text(size = 30))
dev.off()

png(file = "Tuft_Module_Score_Density_Plot.png", width=18, height=12, units = 'in',res=300, type = "cairo")
p <- ggplot(module_trend_point_df.Tuft, aes(Module_Score_Value, fill = Condition_Names, colors = Condition_Names)) +
  geom_density(alpha = 0.5) + scale_fill_manual(values = condition_cols, name = "Condition") + labs(y = "Density", x = "Module Score Value")  + theme(axis.text=element_text(size=25, face = "bold"), axis.title=element_text(size=30,face="bold")) + theme(strip.text.x = element_text(size = 20, face="bold")) 

p +theme(legend.text = element_text(size=25), legend.title = element_text(size = 30))
dev.off()


png(file = "Secretory_Module_Score_Density_Plot.png", width=18, height=12, units = 'in',res=300, type = "cairo")
p <- ggplot(module_trend_point_df.Secretory, aes(Module_Score_Value, fill = Condition_Names, colors = Condition_Names)) +
  geom_density(alpha = 0.5) + scale_fill_manual(values = condition_cols, name = "Condition") + labs(y = "Density", x = "Module Score Value")  + theme(axis.text=element_text(size=25, face = "bold"), axis.title=element_text(size=30,face="bold")) + theme(strip.text.x = element_text(size = 20, face="bold")) 

p +theme(legend.text = element_text(size=25), legend.title = element_text(size = 30))
dev.off()
######


module_trend_df.avg <- module_trend_df %>%
  group_by(Condition_Names, Module_Type) %>%
  summarize(mean = mean(Module_Score_Value),
            sd   = sd(Module_Score_Value),
            mean_p2sd = mean + 2 * sd,
            mean_m2sd = mean - 2 * sd) %>%
  ungroup()



module_trend_df.ATI <- filter(module_trend_df.avg, Module_Type == "ATI_Module_Values")
module_trend_df.ATII <- filter(module_trend_df.avg, Module_Type == "ATII_Module_Values")
module_trend_df.ATII_ATI <- filter(module_trend_df.avg, Module_Type == "ATII_ATI_Module_Values")
module_trend_df.Basal <- filter(module_trend_df.avg, Module_Type == "Basal_Module_Values")
module_trend_df.BASC <- filter(module_trend_df.avg, Module_Type == "BASC_Module_Values")
module_trend_df.Ciliated <- filter(module_trend_df.avg, Module_Type == "Ciliated_Module_Values")
module_trend_df.Hillock <- filter(module_trend_df.avg, Module_Type == "Hillock_Module_Values")
module_trend_df.Secretory <- filter(module_trend_df.avg, Module_Type == "Secretory_Module_Values")
module_trend_df.Tuft <- filter(module_trend_df.avg, Module_Type == "Tuft_Module_Values")
module_trend_df.AbBasal <- filter(module_trend_df.avg, Module_Type == "Aberrant_Basaloid_Module_Values")
setwd("/gpfs/gibbs/project/kaminski/hk738/Storage/Satoshi_Data/BASC_Analysis/BASC_Combined/2025.5.19/Trend_Plot")

plot.2 <-ggplot(data = module_trend_df.AbBasal, 
                aes(x = Condition_Names, 
                    y= mean, 
                    group = Module_Type)) +
  geom_line(aes(color=Module_Type), size = 1.5)+
  geom_point(aes(color=Module_Type)) +
  scale_color_manual(values = c("#FF6666"), labels = "Aberrant Basaloid Module Value") + 
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd, color = Module_Type), width=.2, size = 1.5) 


png(file = "AbBasal_Module_Score_Trend_ErrorBars.png", width=16, height=10, units = 'in',res=300, type = "cairo")
plot.2 + theme_bw() + labs(y = "Average Module Score", x = "Condition", color = "Module Type") + theme(axis.text=element_text(size=15, face = "bold"), axis.title=element_text(size=30,face="bold")) + theme(strip.text.x = element_text(size = 20, face="bold")) + theme(legend.position = "top") +theme(legend.text = element_text(size=15), legend.title = element_text(size = 20))

dev.off()

plot.2 <-ggplot(data = module_trend_df.ATI, 
                aes(x = Condition_Names, 
                    y= mean, 
                    group = Module_Type)) +
  geom_line(aes(color=Module_Type), size = 1.5)+
  geom_point(aes(color=Module_Type)) +
  scale_color_manual(values = c("#987284"), labels = "ATI Module Value") + 
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd, color = Module_Type), width=.2, size = 1.5) 

png(file = "ATI_Module_Score_Trend_ErrorBars.png", width=16, height=10, units = 'in',res=300, type = "cairo")
plot.2 + theme_bw() + labs(y = "Average Module Score", x = "Condition", color = "Module Type")+ theme(axis.text=element_text(size=15, face = "bold"), axis.title=element_text(size=30,face="bold")) + theme(strip.text.x = element_text(size = 20, face="bold")) + theme(legend.position = "top") +theme(legend.text = element_text(size=15), legend.title = element_text(size = 20))
dev.off()

plot.2 <-ggplot(data = module_trend_df.ATII, 
                aes(x = Condition_Names, 
                    y= mean, 
                    group = Module_Type)) +
  geom_line(aes(color=Module_Type), size = 1.5)+
  geom_point(aes(color=Module_Type)) +
  scale_color_manual(values = c("#156064"), labels = "ATII Module Value") + 
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd, color = Module_Type), width=.2, size = 1.5) 

png(file = "ATII_Module_Score_Trend_ErrorBars.png", width=16, height=10, units = 'in',res=300, type = "cairo")
plot.2 + theme_bw() + labs(y = "Average Module Score", x = "Condition", color = "Module Type")+ theme(axis.text=element_text(size=15, face = "bold"), axis.title=element_text(size=30,face="bold")) + theme(strip.text.x = element_text(size = 20, face="bold")) + theme(legend.position = "top") +theme(legend.text = element_text(size=15), legend.title = element_text(size = 20))
dev.off()


plot.2 <-ggplot(data = module_trend_df.ATII_ATI, 
                aes(x = Condition_Names, 
                    y= mean, 
                    group = Module_Type)) +
  geom_line(aes(color=Module_Type), size = 1.5)+
  geom_point(aes(color=Module_Type)) +
  scale_color_manual(values = c("pink"), labels = "ATII-ATI Module Value") + 
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd, color = Module_Type), width=.2, size = 1.5) 

png(file = "ATII_ATI_Module_Score_Trend_ErrorBars.png", width=16, height=10, units = 'in',res=300, type = "cairo")
plot.2 + theme_bw() + labs(y = "Average Module Score", x = "Condition", color = "Module Type")+ theme(axis.text=element_text(size=15, face = "bold"), axis.title=element_text(size=30,face="bold")) + theme(strip.text.x = element_text(size = 20, face="bold")) + theme(legend.position = "top") +theme(legend.text = element_text(size=15), legend.title = element_text(size = 20))
dev.off()

plot.2 <-ggplot(data = module_trend_df.Basal, 
                aes(x = Condition_Names, 
                    y= mean, 
                    group = Module_Type)) +
  geom_line(aes(color=Module_Type), size = 1.5)+
  geom_point(aes(color=Module_Type)) +
  scale_color_manual(values = c("#38AECC"), labels = "Basal Module Value") + 
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd, color = Module_Type), width=.2, size = 1.5) 

png(file = "Basal_Module_Score_Trend_ErrorBars.png", width=16, height=10, units = 'in',res=300, type = "cairo")
plot.2 + theme_bw() + labs(y = "Average Module Score", x = "Condition", color = "Module Type")+ theme(axis.text=element_text(size=15, face = "bold"), axis.title=element_text(size=30,face="bold")) + theme(strip.text.x = element_text(size = 20, face="bold")) + theme(legend.position = "top") +theme(legend.text = element_text(size=15), legend.title = element_text(size = 20))
dev.off()

plot.2 <-ggplot(data = module_trend_df.BASC, 
                aes(x = Condition_Names, 
                    y= mean, 
                    group = Module_Type)) +
  geom_line(aes(color=Module_Type), size = 1.5)+
  geom_point(aes(color=Module_Type)) +
  scale_color_manual(values = c("#D741A7"), labels = "BASC Module Value") + 
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd, color = Module_Type), width=.2, size = 1.5) 

png(file = "BASC_Module_Score_Trend_ErrorBars.png", width=16, height=10, units = 'in',res=300, type = "cairo")
plot.2 + theme_bw() + labs(y = "Average Module Score", x = "Condition", color = "Module Type")+ theme(axis.text=element_text(size=15, face = "bold"), axis.title=element_text(size=30,face="bold")) + theme(strip.text.x = element_text(size = 20, face="bold")) + theme(legend.position = "top") +theme(legend.text = element_text(size=15), legend.title = element_text(size = 20))
dev.off()

plot.2 <-ggplot(data = module_trend_df.Ciliated, 
                aes(x = Condition_Names, 
                    y= mean, 
                    group = Module_Type)) +
  geom_line(aes(color=Module_Type), size = 1.5)+
  geom_point(aes(color=Module_Type)) +
  scale_color_manual(values = c("#A3E51F"), labels = "Ciliated Module Value") + 
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd, color = Module_Type), width=.2, size = 1.5) 

png(file = "Ciliated_Module_Score_Trend_ErrorBars.png", width=16, height=10, units = 'in',res=300, type = "cairo")
plot.2 + theme_bw() + labs(y = "Average Module Score", x = "Condition", color = "Module Type")+ theme(axis.text=element_text(size=15, face = "bold"), axis.title=element_text(size=30,face="bold")) + theme(strip.text.x = element_text(size = 20, face="bold")) + theme(legend.position = "top") +theme(legend.text = element_text(size=15), legend.title = element_text(size = 20))
dev.off()

plot.2 <-ggplot(data = module_trend_df.Hillock, 
                aes(x = Condition_Names, 
                    y= mean, 
                    group = Module_Type)) +
  geom_line(aes(color=Module_Type), size = 1.5)+
  geom_point(aes(color=Module_Type)) +
  scale_color_manual(values = c("#F5B700"), labels = "Hillock Module Value") + 
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd, color = Module_Type), width=.2, size = 1.5) 

png(file = "Hillock_Module_Score_Trend_ErrorBars.png", width=16, height=10, units = 'in',res=300, type = "cairo")
plot.2 + theme_bw() + labs(y = "Average Module Score", x = "Condition", color = "Module Type")+ theme(axis.text=element_text(size=15, face = "bold"), axis.title=element_text(size=30,face="bold")) + theme(strip.text.x = element_text(size = 20, face="bold")) + theme(legend.position = "top") +theme(legend.text = element_text(size=15), legend.title = element_text(size = 20))
dev.off()

plot.2 <-ggplot(data = module_trend_df.Secretory, 
                aes(x = Condition_Names, 
                    y= mean, 
                    group = Module_Type)) +
  geom_line(aes(color=Module_Type), size = 1.5)+
  geom_point(aes(color=Module_Type)) +
  scale_color_manual(values = c("yellow2"), labels = "Secretory Module Value") + 
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd, color = Module_Type), width=.2, size = 1.5) 

png(file = "Secretory_Module_Score_Trend_ErrorBars.png", width=16, height=10, units = 'in',res=300, type = "cairo")
plot.2 + theme_bw() + labs(y = "Average Module Score", x = "Condition", color = "Module Type")+ theme(axis.text=element_text(size=15, face = "bold"), axis.title=element_text(size=30,face="bold")) + theme(strip.text.x = element_text(size = 20, face="bold")) + theme(legend.position = "top") +theme(legend.text = element_text(size=15), legend.title = element_text(size = 20))
dev.off()

plot.2 <-ggplot(data = module_trend_df.Tuft, 
                aes(x = Condition_Names, 
                    y= mean, 
                    group = Module_Type)) +
  geom_line(aes(color=Module_Type), size = 1.5)+
  geom_point(aes(color=Module_Type)) +
  scale_color_manual(values = c("#6DECAF"), labels = "Tuft Module Value") + 
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd, color = Module_Type), width=.2, size = 1.5) 


png(file = "Tuft_Module_Score_Trend_ErrorBars.png", width=16, height=10, units = 'in',res=300, type = "cairo")
plot.2 + theme_bw() + labs(y = "Average Module Score", x = "Condition", color = "Module Type")+ theme(axis.text=element_text(size=15, face = "bold"), axis.title=element_text(size=30,face="bold")) + theme(strip.text.x = element_text(size = 20, face="bold")) + theme(legend.position = "top") +theme(legend.text = element_text(size=15), legend.title = element_text(size = 20))
dev.off()


setwd("/gpfs/gibbs/project/kaminski/hk738/Storage/Satoshi_Data/BASC_Analysis/BASC_Combined/2025.5.19/FeaturePlots")

fix.sc.basal <- scale_color_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")),  limits = c(-1, 2), oob = scales::squish)

png(file = "Basal.ModuleScore.Merged_Epithelium.UniversalScale.png", width=18, height=12, units = 'in',res=300, type = "cairo")
Basal_Module_Atlas <- FeaturePlot(HK_tissue_engineering_epithelial_merged.dataset.2025.5.5, features = "Basal_Module1", pt.size = 1.5, order = T, label = T, repel = T, combine = FALSE, label.size = 8.5)
Basal_Module_Atlas.2 <- lapply(Basal_Module_Atlas, function (x) x + fix.sc.basal)
CombinePlots(Basal_Module_Atlas.2) + NoAxes() + ggtitle("Basal Module Score Feature Plot") + theme(legend.key.size = unit(1, 'cm')) + theme(legend.text = element_text(size=15), theme = "bold")
dev.off()

fix.sc.hillock <- scale_color_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")),  limits = c(-0.5, 2), oob = scales::squish)

png(file = "Hillock.ModuleScore.Merged_Epithelium.UniversalScale.png", width=18, height=12, units = 'in',res=300, type = "cairo")
Hillock_Module_Atlas <- FeaturePlot(HK_tissue_engineering_epithelial_merged.dataset.2025.5.5, features = "Hillock_Module1", pt.size = 1.5, order = T, label = T, repel = T, combine = FALSE, label.size = 8.5)
Hillock_Module_Atlas.2 <- lapply(Hillock_Module_Atlas, function (x) x + fix.sc.hillock)
CombinePlots(Hillock_Module_Atlas.2) + NoAxes() + ggtitle("Hillock Module Score Feature Plot") + theme(legend.key.size = unit(1, 'cm')) + theme(legend.text = element_text(size=15), theme = "bold")
dev.off()

fix.sc.at2 <- scale_color_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")),  limits = c(-0.5, 2), oob = scales::squish)

png(file = "AT2.ModuleScore.Merged_Epithelium.UniversalScale.png", width=18, height=12, units = 'in',res=300, type = "cairo")
AT2_Module_Atlas <- FeaturePlot(HK_tissue_engineering_epithelial_merged.dataset.2025.5.5, features = "ATII_Module1", pt.size = 1.5, order = T, label = T, repel = T, combine = FALSE, label.size = 8.5)
AT2_Module_Atlas.2 <- lapply(AT2_Module_Atlas, function (x) x + fix.sc.at2)
CombinePlots(AT2_Module_Atlas.2) + NoAxes() + ggtitle("ATII Module Score Feature Plot") + theme(legend.key.size = unit(1, 'cm')) + theme(legend.text = element_text(size=15), theme = "bold")
dev.off()

fix.sc.at1 <- scale_color_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")),  limits = c(-.25, .75), oob = scales::squish)

png(file = "ATI.ModuleScore.Merged_Epithelium.UniversalScale.png", width=18, height=12, units = 'in',res=300, type = "cairo")
ATI_Module_Atlas <- FeaturePlot(HK_tissue_engineering_epithelial_merged.dataset.2025.5.5, features = "ATI_Module1", pt.size = 1.5, order = T, label = T, repel = T, combine = FALSE, label.size = 8.5)
ATI_Module_Atlas.2 <- lapply(ATI_Module_Atlas, function (x) x + fix.sc.at1)
CombinePlots(ATI_Module_Atlas.2) + NoAxes() + ggtitle("ATI Module Score Feature Plot") + theme(legend.key.size = unit(1, 'cm')) + theme(legend.text = element_text(size=15), theme = "bold")
dev.off()


fix.sc.secretory <- scale_color_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")),  limits = c(-0.5, 2), oob = scales::squish)

png(file = "Secretory.ModuleScore.Merged_Epithelium.UniversalScale.png", width=18, height=12, units = 'in',res=300, type = "cairo")
Secretory_Module_Atlas <- FeaturePlot(HK_tissue_engineering_epithelial_merged.dataset.2025.5.5, features = "Secretory_Module1", pt.size = 1.5, order = T, label = T, repel = T, combine = FALSE, label.size = 8.5)
Secretory_Module_Atlas.2 <- lapply(Secretory_Module_Atlas, function (x) x + fix.sc.secretory)
CombinePlots(Secretory_Module_Atlas.2) + NoAxes() + ggtitle("Secretory Module Score Feature Plot") + theme(legend.key.size = unit(1, 'cm')) + theme(legend.text = element_text(size=15), theme = "bold")
dev.off()


fix.sc.ciliated <- scale_color_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")),  limits = c(-.25, 1), oob = scales::squish)

png(file = "Ciliated.ModuleScore.Merged_Epithelium.UniversalScale.png", width=18, height=12, units = 'in',res=300, type = "cairo")
Ciliated_Module_Atlas <- FeaturePlot(HK_tissue_engineering_epithelial_merged.dataset.2025.5.5, features = "Ciliated_Module1", pt.size = 1.5, order = T, label = T, repel = T, combine = FALSE, label.size = 8.5)
Ciliated_Module_Atlas.2 <- lapply(Ciliated_Module_Atlas, function (x) x + fix.sc.ciliated)
CombinePlots(Ciliated_Module_Atlas.2) + NoAxes() + ggtitle("Ciliated Module Score Feature Plot") + theme(legend.key.size = unit(1, 'cm')) + theme(legend.text = element_text(size=15), theme = "bold")
dev.off()


fix.sc.tuft <- scale_color_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")),  limits = c(-.25, .75), oob = scales::squish)

png(file = "Tuft.ModuleScore.Merged_Epithelium.UniversalScale.png", width=18, height=12, units = 'in',res=300, type = "cairo")
Tuft_Module_Atlas <- FeaturePlot(HK_tissue_engineering_epithelial_merged.dataset.2025.5.5, features = "Tuft_Module1", pt.size = 1.5, order = T, label = T, repel = T, combine = FALSE, label.size = 8.5)
Tuft_Module_Atlas.2 <- lapply(Tuft_Module_Atlas, function (x) x + fix.sc.tuft)
CombinePlots(Tuft_Module_Atlas.2) + NoAxes() + ggtitle("Tuft Module Score Feature Plot") + theme(legend.key.size = unit(1, 'cm')) + theme(legend.text = element_text(size=15), theme = "bold")
dev.off()


fix.sc.basc <- scale_color_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")),  limits = c(-0.25, 0.75), oob = scales::squish)

png(file = "BASC.ModuleScore.Merged_Epithelium.UniversalScale.png", width=18, height=12, units = 'in',res=300, type = "cairo")
BASC_Module_Atlas <- FeaturePlot(HK_tissue_engineering_epithelial_merged.dataset.2025.5.5, features = "BASC_Module1", pt.size = 1.5, order = T, label = T, repel = T, combine = FALSE, label.size = 8.5)
BASC_Module_Atlas.2 <- lapply(BASC_Module_Atlas, function (x) x + fix.sc.basc)
CombinePlots(BASC_Module_Atlas.2) + NoAxes() + ggtitle("BASC Module Score Feature Plot") + theme(legend.key.size = unit(1, 'cm')) + theme(legend.text = element_text(size=15), theme = "bold")
dev.off()


fix.sc.at2_at1 <- scale_color_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")),  limits = c(-.25, 1), oob = scales::squish)

png(file = "ATII_ATI_ModuleScore.Merged_Epithelium.UniversalScale.png", width=18, height=12, units = 'in',res=300, type = "cairo")
ATII_ATI_Module_Atlas <- FeaturePlot(HK_tissue_engineering_epithelial_merged.dataset.2025.5.5, features = "ATII_ATI_Module1", pt.size = 1.5, order = T, label = T, repel = T, combine = FALSE, label.size = 8.5)
ATII_ATI_Module_Atlas.2 <- lapply(ATII_ATI_Module_Atlas, function (x) x + fix.sc.at2_at1)
CombinePlots(ATII_ATI_Module_Atlas.2) + NoAxes() + ggtitle("ATII-ATI Module Score Feature Plot") + theme(legend.key.size = unit(1, 'cm')) + theme(legend.text = element_text(size=15), theme = "bold")
dev.off()

fix.sc.ab_basaloid <- scale_color_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")),  limits = c(-0.25, 0.5), oob = scales::squish)


png(file = "Ab_Basal.ModuleScore.Merged_Epithelium.UniversalScale.png", width=18, height=12, units = 'in',res=300, type = "cairo")
AbBasaloid_Module_Atlas <- FeaturePlot(HK_tissue_engineering_epithelial_merged.dataset.2025.5.5, features = "Aberrant_Basaloid_Module1", pt.size = 1.5, order = T, label = T, repel = T, combine = FALSE, label.size = 8.5)
AbBasaloid_Module_Atlas.2 <- lapply(AbBasaloid_Module_Atlas, function (x) x + fix.sc.ab_basaloid)
CombinePlots(AbBasaloid_Module_Atlas.2) + NoAxes() + ggtitle("Aberrant Basal Module Score Feature Plot") + theme(legend.key.size = unit(1, 'cm')) + theme(legend.text = element_text(size=15), theme = "bold")
dev.off()

dat <- HK_tissue_engineering_epithelial_merged.dataset.2025.5.5
rownames(dat@assays$RNA@layers$dat) <- rownames(dat)

Idents(dat) <- dat$ConditionRearrangedNewName

to.plot <- dat@meta.data



####Vln Plot

png(file = 'Aberrant_Basaloid_Mod_VlnPlot.png', width=12, height=8, units = 'in',res=300, type = "cairo")
plot <- ggplot(to.plot, 
               aes(x=ConditionRearrangedNewName, 
                   y=Aberrant_Basaloid_Module1, fill = ConditionRearrangedNewName)) +
  geom_violin(aes(group = ConditionRearrangedNewName),linewidth = 0.1)+
  stat_summary(fun = "mean",
               geom = "point",
               color = "black",size = 0.75)+
  theme_classic()+ 
  theme(axis.text.x = element_text(angle = 35, vjust = 1, hjust=1, face = "bold", size = 10))+
  ylab('Normalized Expression')+ xlab('Condition') +
  ggtitle("Aberrant Basaloid Module Expression Across Conditions") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(axis.title.x = element_text(size = 20, face = "bold")) +
  theme(axis.title.y = element_text(size = 20, face = "bold")) +
  theme(title = element_text(size = 20, face = "bold")) +
  theme(axis.text.y = element_text(size = 10, face = "bold"))

plot+ scale_fill_manual(values=c("red", "blue", "green", "yellow", "purple", "orange")) + NoLegend()
dev.off()

png(file = 'Basal_Mod_VlnPlot.png', width=12, height=8, units = 'in',res=300, type = "cairo")
plot <- ggplot(to.plot, 
               aes(x=ConditionRearrangedNewName, 
                   y=Basal_Module1, fill = ConditionRearrangedNewName)) +
  geom_violin(aes(group = ConditionRearrangedNewName),linewidth = 0.1)+
  stat_summary(fun = "mean",
               geom = "point",
               color = "black",size = 0.75)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 35, vjust = 1, hjust=1, face = "bold", size = 10))+
  ylab('Normalized Expression')+ xlab('Condition') +
  ggtitle("Basal Module Expression Across Conditions") +
  theme(plot.title = element_text(hjust = 0.5))+ 
  theme(axis.title.x = element_text(size = 20, face = "bold")) +
  theme(axis.title.y = element_text(size = 20, face = "bold"))+
  theme(title = element_text(size = 20, face = "bold")) +
  theme(axis.text.y = element_text(size = 10, face = "bold"))

plot + scale_fill_manual(values=c("red", "blue", "green", "yellow", "purple", "orange")) + NoLegend()
dev.off()


png(file = 'ATI_Mod_VlnPlot.png', width=12, height=8, units = 'in',res=300, type = "cairo")
plot <- ggplot(to.plot, 
               aes(x=ConditionRearrangedNewName, 
                   y=ATI_Module1, fill = ConditionRearrangedNewName)) +
  geom_violin(aes(group = ConditionRearrangedNewName),linewidth = 0.1)+
  stat_summary(fun = "mean",
               geom = "point",
               color = "black",size = 0.75)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 35, vjust = 1, hjust=1, face = "bold", size = 10))+
  ylab('Normalized Expression')+ xlab('Condition') +
  ylim(-0.25, 0.85) +
  ggtitle("ATI Module Expression Across Conditions") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(axis.title.x = element_text(size = 20, face = "bold")) +
  theme(axis.title.y = element_text(size = 20, face = "bold"))+
  theme(title = element_text(size = 20, face = "bold")) +
  theme(axis.text.y = element_text(size = 10, face = "bold"))

plot + scale_fill_manual(values=c("red", "blue", "green", "yellow", "purple", "orange")) + NoLegend()
dev.off()


png(file = 'ATII_Mod_VlnPlot.png', width=12, height=8, units = 'in',res=300, type = "cairo")
plot <- ggplot(to.plot, 
               aes(x=ConditionRearrangedNewName, 
                   y=ATII_Module1, fill = ConditionRearrangedNewName)) +
  geom_violin(aes(group = ConditionRearrangedNewName),linewidth = 0.1)+
  stat_summary(fun = "mean",
               geom = "point",
               color = "black",size = 0.75)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 35, vjust = 1, hjust=1, face = "bold", size = 10))+
  ylab('Normalized Expression')+ xlab('Condition') +
  ggtitle("ATII Module Expression Across Conditions") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(axis.title.x = element_text(size = 20, face = "bold")) +
  theme(axis.title.y = element_text(size = 20, face = "bold"))+
  theme(title = element_text(size = 20, face = "bold")) +
  theme(axis.text.y = element_text(size = 10, face = "bold"))

plot + scale_fill_manual(values=c("red", "blue", "green", "yellow", "purple", "orange")) + NoLegend()
dev.off()


png(file = 'ATII_ATI_Mod_VlnPlot.png', width=12, height=8, units = 'in',res=300, type = "cairo")
plot <- ggplot(to.plot, 
               aes(x=ConditionRearrangedNewName, 
                   y=ATII_ATI_Module1, fill = ConditionRearrangedNewName)) +
  geom_violin(aes(group = ConditionRearrangedNewName),linewidth = 0.1)+
  stat_summary(fun = "mean",
               geom = "point",
               color = "black",size = 0.75)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 35, vjust = 1, hjust=1, face = "bold", size = 10))+
  ylab('Normalized Expression')+ xlab('Condition') +
  ggtitle("ATII-ATI Module Expression Across Conditions") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(axis.title.x = element_text(size = 20, face = "bold")) +
  theme(axis.title.y = element_text(size = 20, face = "bold"))+
  theme(title = element_text(size = 20, face = "bold")) +
  theme(axis.text.y = element_text(size = 10, face = "bold"))

plot + scale_fill_manual(values=c("red", "blue", "green", "yellow", "purple", "orange")) + NoLegend()
dev.off()


png(file = 'BASC_Mod_VlnPlot.png', width=12, height=8, units = 'in',res=300, type = "cairo")
plot <- ggplot(to.plot, 
               aes(x=ConditionRearrangedNewName, 
                   y=BASC_Module1, fill = ConditionRearrangedNewName)) +
  geom_violin(aes(group = ConditionRearrangedNewName),linewidth = 0.1)+
  stat_summary(fun = "mean",
               geom = "point",
               color = "black",size = 0.75)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 35, vjust = 1, hjust=1, face = "bold", size = 10))+
  ylab('Normalized Expression')+ xlab('Condition') +
  ggtitle("BASC Module Expression Across Conditions") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(axis.title.x = element_text(size = 20, face = "bold")) +
  theme(axis.title.y = element_text(size = 20, face = "bold"))+
  theme(title = element_text(size = 20, face = "bold")) +
  theme(axis.text.y = element_text(size = 10, face = "bold"))

plot + scale_fill_manual(values=c("red", "blue", "green", "yellow", "purple", "orange")) + NoLegend()
dev.off()


png(file = 'Ciliated_Mod_VlnPlot.png', width=12, height=8, units = 'in',res=300, type = "cairo")
plot <- ggplot(to.plot, 
               aes(x=ConditionRearrangedNewName, 
                   y=Ciliated_Module1, fill = ConditionRearrangedNewName)) +
  geom_violin(aes(group = ConditionRearrangedNewName),linewidth = 0.1)+
  stat_summary(fun = "mean",
               geom = "point",
               color = "black",size = 0.75)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 35, vjust = 1, hjust=1, face = "bold", size = 10))+
  ylab('Normalized Expression')+ xlab('Condition') +
  ggtitle("Ciliated Module Expression Across Conditions") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(axis.title.x = element_text(size = 20, face = "bold")) +
  theme(axis.title.y = element_text(size = 20, face = "bold"))+
  theme(title = element_text(size = 20, face = "bold")) +
  theme(axis.text.y = element_text(size = 10, face = "bold"))

plot + scale_fill_manual(values=c("red", "blue", "green", "yellow", "purple", "orange")) + NoLegend()
dev.off()


png(file = 'Hillock_Mod_VlnPlot.png', width=12, height=8, units = 'in',res=300, type = "cairo")
plot <- ggplot(to.plot, 
               aes(x=ConditionRearrangedNewName, 
                   y=Hillock_Module1, fill = ConditionRearrangedNewName)) +
  geom_violin(aes(group = ConditionRearrangedNewName),linewidth = 0.1)+
  stat_summary(fun = "mean",
               geom = "point",
               color = "black",size = 0.75)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 35, vjust = 1, hjust=1, face = "bold", size = 10))+
  ylab('Normalized Expression')+ xlab('Condition') +
  ggtitle("Hillock Module Expression Across Conditions") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(axis.title.x = element_text(size = 20, face = "bold")) +
  theme(axis.title.y = element_text(size = 20, face = "bold"))+
  theme(title = element_text(size = 20, face = "bold")) +
  theme(axis.text.y = element_text(size = 10, face = "bold"))

plot + scale_fill_manual(values=c("red", "blue", "green", "yellow", "purple", "orange")) + NoLegend()
dev.off()


png(file = 'Secretory_Mod_VlnPlot.png', width=12, height=8, units = 'in',res=300, type = "cairo")
plot <- ggplot(to.plot, 
               aes(x=ConditionRearrangedNewName, 
                   y=Secretory_Module1, fill = ConditionRearrangedNewName)) +
  geom_violin(aes(group = ConditionRearrangedNewName),linewidth = 0.1)+
  stat_summary(fun = "mean",
               geom = "point",
               color = "black",size = 0.75)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 35, vjust = 1, hjust=1, face = "bold", size = 10))+
  ylab('Normalized Expression')+ xlab('Condition') +
  ggtitle("Secretory Module Expression Across Conditions") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(axis.title.x = element_text(size = 20, face = "bold")) +
  theme(axis.title.y = element_text(size = 20, face = "bold"))+
  theme(title = element_text(size = 20, face = "bold")) +
  theme(axis.text.y = element_text(size = 10, face = "bold"))

plot + scale_fill_manual(values=c("red", "blue", "green", "yellow", "purple", "orange")) + NoLegend()
dev.off()


png(file = 'Tuft_Mod_VlnPlot.png', width=12, height=8, units = 'in',res=300, type = "cairo")
plot <- ggplot(to.plot, 
               aes(x=ConditionRearrangedNewName, 
                   y=Tuft_Module1, fill = ConditionRearrangedNewName)) +
  geom_violin(aes(group = ConditionRearrangedNewName),linewidth = 0.1)+
  stat_summary(fun = "mean",
               geom = "point",
               color = "black",size = 0.75)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 35, vjust = 1, hjust=1, face = "bold", size = 10))+
  ylab('Normalized Expression')+ xlab('Condition') +
  ggtitle("Tuft Module Expression Across Conditions") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(axis.title.x = element_text(size = 20, face = "bold")) +
  theme(axis.title.y = element_text(size = 20, face = "bold"))+
  theme(title = element_text(size = 20, face = "bold")) +
  theme(axis.text.y = element_text(size = 10, face = "bold"))

plot + scale_fill_manual(values=c("red", "blue", "green", "yellow", "purple", "orange")) + NoLegend()
dev.off()









