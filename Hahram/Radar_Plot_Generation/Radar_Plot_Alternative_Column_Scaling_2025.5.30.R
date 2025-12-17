# import the library
library(caret)

# dataset
data = data.frame(var1 = c(120,345,145,122,596,285,211),
                  var2 = c(10,15,45,22,53,28,12),
                  var3 = c(-34,0.05,0.15,0.12,-6,0.85,0.11))

# custom function to implement min max scaling
minMax <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

#normalise data using custom function
normalisedMydata <- as.data.frame(lapply(data, minMax))
head(normalisedMydata)



module_trend_df.avg_unpivoted.altScaling <- pivot_wider(module_trend_df.avg.filtered, names_from = Condition_Names, values_from = mean) %>% group_by(Module_Type)
module_trend_df.avg_unpivoted.altScaling.zNorm <- module_trend_df.avg_unpivoted.altScaling

module_trend_df.avg_unpivoted.altScaling.zNorm$Embryo <- scale(module_trend_df.avg_unpivoted.altScaling.zNorm$Embryo, center = TRUE, scale = TRUE)
module_trend_df.avg_unpivoted.altScaling.zNorm$`St. Dclk1` <- scale(module_trend_df.avg_unpivoted.altScaling.zNorm$`St. Dclk1`, center = TRUE, scale = TRUE)
module_trend_df.avg_unpivoted.altScaling.zNorm$`P0-2D`<- scale(module_trend_df.avg_unpivoted.altScaling.zNorm$`P0-2D`, center = TRUE, scale = TRUE)
module_trend_df.avg_unpivoted.altScaling.zNorm$`P3-2D` <- scale(module_trend_df.avg_unpivoted.altScaling.zNorm$`P3-2D`)
module_trend_df.avg_unpivoted.altScaling.zNorm$`24 Hour` <- scale(module_trend_df.avg_unpivoted.altScaling.zNorm$`24 Hour`, center = TRUE, scale = TRUE)
module_trend_df.avg_unpivoted.altScaling.zNorm$`Five Day` <- scale(module_trend_df.avg_unpivoted.altScaling.zNorm$`Five Day`, center = TRUE, scale = TRUE)

 df <- module_trend_df.avg_unpivoted.altScaling.zNorm
df2 <- data.frame(t(df[-1]))

df2$Condition_Names <- c("Embryo", "St Dclk1", "P0-2D", "P3-2D", "24 Hour", "Five Day")
df2 <- df2 %>% relocate(Condition_Names)
colnames(df2) <- c("Condition_Names","ATII_ATI", "ATII", "ATI", "AbBasal", "BASC", "Basal", "Ciliated", "Hillock", "Secretory", "Tuft")


df2 <- df2[, c(1, 11, 10, 4, 5, 7, 2, 6, 3, 9, 8)]
df2 <- df2[c(1, 2, 3, 4, 5, 6),]

rownames(df2) <- NULL
png(file = "Avg_Module_Radar_ZNorm_AltScaling_Plot.Condition.png", width=18, height=12, units = 'in',res=300, type = "cairo")
ggradar(df2, grid.min = -3, grid.max = 3,  fill = T, group.point.size = 0, fill.alpha = 0.25, group.line.width = 0.5, values.radar = c("-3", "0", "3")) + scale_fill_manual(labels = c("Embryo", "St Dclk1", "P0-2D", "P3-2D", "24 Hour", "Five Day"), breaks = c("Embryo", "St Dclk1", "P0-2D", "P3-2D", "24 Hour", "Five Day"), values = c("red", "blue", "green", "yellow", "purple", "orange"))
dev.off()

