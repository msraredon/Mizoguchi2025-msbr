install.packages("fmsb")
library(fmsb)
radar.df <- module_trend_df.avg_unpivoted
radarchart(module_trend_df.avg_unpivoted)
max_min <- data.frame(
  Condition_Names =  c(1.5, -1), Tuft = c(1.5,-1), Secretory = c(1.5,-1), ATI = c(1.5,-1),
  AbBasal = c(1.5,-1), Basal = c(1.5,-1), `ATII-ATI` = c(1.5,-1),
  BASC = c(1.5,-1), ATII = c(1.5,-1), Hillock = c(1.5,-1), Ciliated = c(1.5, -1)
)



radar.df = subset(radar.df, select = -c(Condition_Names) )

rownames(max_min) <- c("Max", "Min")
colnames(max_min) <- c("Condition_Names", "Tuft", "Secretory", "ATI", "AbBasal", "Basal", "ATII-ATI", "BASC", "ATII", "Hillock", "Ciliated")
module_trend_df.avg_unpivoted <- rbind(max_min, module_trend_df.avg_unpivoted)
rownames(radar.df) <- c("Max", "Min", "Embryo", "St. Dclk1", "P0-2D", "P3-2D", "24 Hour", "Five Day")


radar.df.embryo <- radar.df[c("Max", "Min", "Embryo", "St. Dclk1", "P0-2D", "P3-2D", "24 Hour", "Five Day"), ]
radarchart(radar.df.embryo, caxislabels = c(-1, 0, .5, 1, 1.5))

           

ggradar(module_trend_df.avg_unpivoted, grid.min = -1.0, grid.max = 1.5, fill = T, group.point.size = 0, group.colours = c("red", "blue", "green", "yellow", "purple", "orange"), fill.alpha = 0.25, group.line.width = 0.5, grid.min = -1.5, grid.max = 1.5)
