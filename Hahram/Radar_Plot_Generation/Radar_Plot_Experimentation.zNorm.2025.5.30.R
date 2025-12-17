png(file = "Avg_Module_Radar_Original_Plot.Condition.png", width=18, height=12, units = 'in',res=300, type = "cairo")
ggradar(module_trend_df.avg_unpivoted, grid.min = -1.5, grid.max = 1.5, fill = T, group.point.size = 0, group.colours = c("red", "blue", "green", "yellow", "purple", "orange"), fill.alpha = 0.25, group.line.width = 0.5,
        values.radar = c("-1.5", "0","1.5"))
dev.off()



module_trend_df.avg_unpivoted.Embryo <- filter(module_trend_df.avg_unpivoted, Condition_Names == "Embryo")
module_trend_df.avg_unpivoted.Start <-filter(module_trend_df.avg_unpivoted, Condition_Names == "St. Dclk1")
module_trend_df.avg_unpivoted.P0 <-  filter(module_trend_df.avg_unpivoted, Condition_Names == "P0-2D")
module_trend_df.avg_unpivoted.P3 <- filter(module_trend_df.avg_unpivoted, Condition_Names == "P3-2D")
module_trend_df.avg_unpivoted.24Hour <- filter(module_trend_df.avg_unpivoted, Condition_Names == "24 Hour")
module_trend_df.avg_unpivoted.FiveDay <- filter(module_trend_df.avg_unpivoted, Condition_Names == "Five Day")


png(file = "Avg_Module_Radar_Original_Plot.Condition.Split.png", width=18, height=12, units = 'in',res=300, type = "cairo")
m1<-ggradar(module_trend_df.avg_unpivoted.Embryo, grid.min = -1.5, grid.max = 1.5, fill = T, group.point.size = 0, group.colours = c("red"), fill.alpha = 0.5, group.line.width = 0.5, values.radar = c("-1.5", "0", "1.5")) + ggtitle("Embryo")
m1.5<-ggradar(module_trend_df.avg_unpivoted.Start, grid.min = -1.5, grid.max = 1.5, fill = T, group.point.size = 0, group.colours = c("blue"), fill.alpha = .5, group.line.width = 0.5, values.radar = c("-1.5", "0", "1.5"))+ ggtitle("St. Dclk1")
m3<-ggradar(module_trend_df.avg_unpivoted.P0, grid.min = -1.5, grid.max = 1.5, fill = T, group.point.size = 0, group.colours = c("green"), fill.alpha = .5, group.line.width = 0.5, values.radar = c("-1.5", "0", "1.5"))+ ggtitle("P0-1.5D")
m4<-ggradar(module_trend_df.avg_unpivoted.P3, grid.min = -1.5, grid.max = 1.5, fill = T, group.point.size = 0, group.colours = c("yellow"), fill.alpha = .5, group.line.width = 0.5, values.radar = c("-1.5", "0", "1.5"))+ ggtitle("P0-3D")
m5<-ggradar(module_trend_df.avg_unpivoted.24Hour, grid.min = -1.5, grid.max = 1.5, fill = T, group.point.size = 0, group.colours = c("purple"), fill.alpha = .5, group.line.width = 0.5, values.radar = c("-1.5", "0", "1.5"))+ ggtitle("1.54 Hour")
m6<-ggradar(module_trend_df.avg_unpivoted.FiveDay, grid.min = -1.5, grid.max = 1.5, fill = T, group.point.size = 0, group.colours = c("orange"), fill.alpha = .5, group.line.width = 0.5, values.radar = c("-1.5", "0", "1.5"))+ ggtitle("Five Day")
cowplot::plot_grid(m1, m1.5, m3, m4, m5, m6, ncol = 3, nrow = 2)
dev.off()








module_trend_df.avg_unpivoted.zNorm <- module_trend_df.avg_unpivoted

module_trend_df.avg_unpivoted.zNorm$"AbBasal" <- scale(module_trend_df.avg_unpivoted.zNorm$AbBasal, center = TRUE, scale = TRUE)
module_trend_df.avg_unpivoted.zNorm$"ATI" <- scale(module_trend_df.avg_unpivoted.zNorm$ATI, center = TRUE, scale = TRUE)
module_trend_df.avg_unpivoted.zNorm$"ATII" <- scale(module_trend_df.avg_unpivoted.zNorm$ATII, center = TRUE, scale = TRUE)
module_trend_df.avg_unpivoted.zNorm$"ATII-ATI" <- scale(module_trend_df.avg_unpivoted.zNorm$`ATII-ATI`, center = TRUE, scale = TRUE)
module_trend_df.avg_unpivoted.zNorm$"Basal" <- scale(module_trend_df.avg_unpivoted.zNorm$Basal, center = TRUE, scale = TRUE)
module_trend_df.avg_unpivoted.zNorm$"BASC" <- scale(module_trend_df.avg_unpivoted.zNorm$BASC, center = TRUE, scale = TRUE)
module_trend_df.avg_unpivoted.zNorm$"Ciliated" <- scale(module_trend_df.avg_unpivoted.zNorm$Ciliated, center = TRUE, scale = TRUE)
module_trend_df.avg_unpivoted.zNorm$"Hillock" <- scale(module_trend_df.avg_unpivoted.zNorm$Hillock, center = TRUE, scale = TRUE)
module_trend_df.avg_unpivoted.zNorm$"Secretory" <- scale(module_trend_df.avg_unpivoted.zNorm$Secretory, center = TRUE, scale = TRUE)
module_trend_df.avg_unpivoted.zNorm$"Tuft" <- scale(module_trend_df.avg_unpivoted.zNorm$Tuft, center = TRUE, scale = TRUE)



png(file = "Avg_Module_Radar_ZNorm_Plot.Condition.png", width=18, height=12, units = 'in',res=300, type = "cairo")
ggradar(module_trend_df.avg_unpivoted.zNorm, grid.min = -2, grid.max = 2,  fill = T, group.point.size = 0, group.colours = c("red", "blue", "green", "yellow", "purple", "orange"), fill.alpha = 0.25, group.line.width = 0.5, values.radar = c("-2", "0", "2"))
dev.off()

module_trend_df.avg_unpivoted.zNorm.Embryo <- filter(module_trend_df.avg_unpivoted.zNorm, Condition_Names == "Embryo")
module_trend_df.avg_unpivoted.zNorm.Start <-filter(module_trend_df.avg_unpivoted.zNorm, Condition_Names == "St. Dclk1")
module_trend_df.avg_unpivoted.zNorm.P0 <-  filter(module_trend_df.avg_unpivoted.zNorm, Condition_Names == "P0-2D")
module_trend_df.avg_unpivoted.zNorm.P3 <- filter(module_trend_df.avg_unpivoted.zNorm, Condition_Names == "P3-2D")
module_trend_df.avg_unpivoted.zNorm.24Hour <- filter(module_trend_df.avg_unpivoted.zNorm, Condition_Names == "24 Hour")
module_trend_df.avg_unpivoted.zNorm.FiveDay <- filter(module_trend_df.avg_unpivoted.zNorm, Condition_Names == "Five Day")


png(file = "Avg_Module_Radar_ZNorm_Plot.Condition.Split.png", width=18, height=12, units = 'in',res=300, type = "cairo")
m1<-ggradar(module_trend_df.avg_unpivoted.zNorm.Embryo, grid.min = -2, grid.max = 2, fill = T, group.point.size = 0, group.colours = c("red"), fill.alpha = 0.5, group.line.width = 0.5, values.radar = c("-2", "0", "2")) + ggtitle("Embryo")
m2<-ggradar(module_trend_df.avg_unpivoted.zNorm.Start, grid.min = -2, grid.max = 2, fill = T, group.point.size = 0, group.colours = c("blue"), fill.alpha = .5, group.line.width = 0.5, values.radar = c("-2", "0", "2"))+ ggtitle("St. Dclk1")
m3<-ggradar(module_trend_df.avg_unpivoted.zNorm.P0, grid.min = -2, grid.max = 2, fill = T, group.point.size = 0, group.colours = c("green"), fill.alpha = .5, group.line.width = 0.5, values.radar = c("-2", "0", "2"))+ ggtitle("P0-2D")
m4<-ggradar(module_trend_df.avg_unpivoted.zNorm.P3, grid.min = -2, grid.max = 2, fill = T, group.point.size = 0, group.colours = c("yellow"), fill.alpha = .5, group.line.width = 0.5, values.radar = c("-2", "0", "2"))+ ggtitle("P0-3D")
m5<-ggradar(module_trend_df.avg_unpivoted.zNorm.24Hour, grid.min = -2, grid.max = 2, fill = T, group.point.size = 0, group.colours = c("purple"), fill.alpha = .5, group.line.width = 0.5, values.radar = c("-2", "0", "2"))+ ggtitle("24 Hour")
m6<-ggradar(module_trend_df.avg_unpivoted.zNorm.FiveDay, grid.min = -2, grid.max = 2, fill = T, group.point.size = 0, group.colours = c("orange"), fill.alpha = .5, group.line.width = 0.5, values.radar = c("-2", "0", "2"))+ ggtitle("Five Day")
cowplot::plot_grid(m1, m2, m3, m4, m5, m6, ncol = 3, nrow = 2)
dev.off()

