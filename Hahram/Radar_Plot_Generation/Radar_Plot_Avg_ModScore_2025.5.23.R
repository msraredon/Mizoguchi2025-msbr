install.packages("devtools")
devtools::install_github("ricardo-bion/ggradar")
library(ggradar)

module_trend_df.avg <- module_trend_df %>%
  group_by(Condition_Names, Module_Type) %>%
  summarize(mean = mean(Module_Score_Value),
            sd   = sd(Module_Score_Value),
            mean_p2sd = mean + 2 * sd,
            mean_m2sd = mean - 2 * sd) %>%
  ungroup()

module_trend_df.avg <- module_trend_df.avg %>% group_by(Condition_Names)

module_trend_df.avg.filtered <- module_trend_df.avg[, c("Condition_Names", "Module_Type", "mean")]

module_trend_df.avg_unpivoted <- pivot_wider(module_trend_df.avg.filtered, names_from = Module_Type, values_from = mean) %>% group_by(Condition_Names)

module_trend_df.avg_unpivoted<- rename(module_trend_df.avg_unpivoted, 
                                       "ATI" = "ATI_Module_Values",
                                       "ATII-ATI" = "ATII_ATI_Module_Values", 
                                       "ATII" = "ATII_Module_Values",
                                       "Basal" = "Basal_Module_Values", 
                                       "BASC" = "BASC_Module_Values",
                                       "Ciliated" = "Ciliated_Module_Values",
                                       "Tuft" = "Tuft_Module_Values",
                                       "Secretory" = "Secretory_Module_Values", 
                                       "AbBasal" = "Aberrant_Basaloid_Module_Values",
                                       "Hillock" = "Hillock_Module_Values"  
                                       )
colnames(module_trend_df.avg_unpivoted)
module_trend_df.avg_unpivoted <- module_trend_df.avg_unpivoted[, c(1, 11, 10, 4, 5, 7, 2, 6, 3, 9, 8)]


module_trend_df.avg_unpivoted.Embryo <- filter(module_trend_df.avg_unpivoted, Condition_Names == "Embryo")
module_trend_df.avg_unpivoted.Start <-filter(module_trend_df.avg_unpivoted, Condition_Names == "St. Dclk1")
module_trend_df.avg_unpivoted.P0 <-  filter(module_trend_df.avg_unpivoted, Condition_Names == "P0-2D")
module_trend_df.avg_unpivoted.P3 <- filter(module_trend_df.avg_unpivoted, Condition_Names == "P3-2D")
module_trend_df.avg_unpivoted.24Hour <- filter(module_trend_df.avg_unpivoted, Condition_Names == "24 Hour")
module_trend_df.avg_unpivoted.FiveDay <- filter(module_trend_df.avg_unpivoted, Condition_Names == "Five Day")

png(file = "Avg_Module_Radar_Plot.png", width=18, height=12, units = 'in',res=300, type = "cairo")
ggradar(module_trend_df.avg_unpivoted, grid.min = -1.0, grid.max = 1.5, fill = T, group.point.size = 0, group.colours = c("red", "blue", "green", "yellow", "purple", "orange"), fill.alpha = 0.25, group.line.width = 0.5)
dev.off()


png(file = "Avg_Module_Radar_Plot.Condition.png", width=18, height=12, units = 'in',res=300, type = "cairo")
m1<-ggradar(module_trend_df.avg_unpivoted.Embryo, grid.min = -1.0, grid.max = 1.5, fill = T, group.point.size = 0, group.colours = c("red"), fill.alpha = 0.5, group.line.width = 0.5) + ggtitle("Embryo")
m2<-ggradar(module_trend_df.avg_unpivoted.Start, grid.min = -1.0, grid.max = 1.5, fill = T, group.point.size = 0, group.colours = c("blue"), fill.alpha = .5, group.line.width = 0.5)+ ggtitle("St. Dclk1")
m3<-ggradar(module_trend_df.avg_unpivoted.P0, grid.min = -1.0, grid.max = 1.5, fill = T, group.point.size = 0, group.colours = c("green"), fill.alpha = .5, group.line.width = 0.5)+ ggtitle("P0-2D")
m4<-ggradar(module_trend_df.avg_unpivoted.P3, grid.min = -1.0, grid.max = 1.5, fill = T, group.point.size = 0, group.colours = c("yellow"), fill.alpha = .5, group.line.width = 0.5)+ ggtitle("P0-3D")
m5<-ggradar(module_trend_df.avg_unpivoted.24Hour, grid.min = -1.0, grid.max = 1.5, fill = T, group.point.size = 0, group.colours = c("purple"), fill.alpha = .5, group.line.width = 0.5)+ ggtitle("24 Hour")
m6<-ggradar(module_trend_df.avg_unpivoted.FiveDay, grid.min = -1.0, grid.max = 1.5, fill = T, group.point.size = 0, group.colours = c("orange"), fill.alpha = .5, group.line.width = 0.5)+ ggtitle("Five Day")
cowplot::plot_grid(m1, m2, m3, m4, m5, m6, ncol = 3, nrow = 2)
dev.off()

numero <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
numerodos <- c(1, 2, 3)
df.names<- c("ATI", "ATII", "ATII_ATI",         
 "Basal", "BASC", "Ciliated",         
"Hillock", "Secretory","Tuft",             
"Aberrant_Basaloid")  
df.names.2 <- c("Basal", "Hillock", "Aberrant_Basaloid")
sub.df <- module_trend_plot[, c("Basal_Module_Values", "Hillock_Module_Values", "Aberrant_Basaloid_Module_Values")]
for(i in numero){
  for(j in numerodos){
    png(paste0(df.names[i], df.names.2[j], "_FeatPlot.png"), width=20, height=14, units = 'in',res=300, type = "cairo")
  p <- ggplot(module_trend_plot, aes(module_trend_plot[i], sub.df[j])) +
    geom_point(aes(color = Condition_Names), size = 3, alpha = 0.6) +
    scale_color_manual(values = c("red", "blue", "yellow", "green", "purple", "orange"))
  p + theme(axis.text=element_text(size=14), axis.title=element_text(size=20,face="bold")) + theme(legend.text = element_text(size=20), legend.key.size = unit(.5, 'cm'))  + guides(colour = guide_legend(override.aes = list(size=10)))  + theme(legend.title=element_text(size=20))  + theme(legend.title=element_text(size=20)) + theme(legend.position = "inside", legend.position.inside =  c(.90, .85))
  print(p)
  dev.off()
  }
}
