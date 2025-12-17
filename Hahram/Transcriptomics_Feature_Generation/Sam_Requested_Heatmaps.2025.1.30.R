###We now have a good formula for how to get somewhat specific genes
##What sam wants to see (and what I will try and refine)
# up, down, up, up, up, down
# Down, up, down, down, down, up
# Down, down, up, up, up, up
# Up, up, down, down, down

### First: Find Features Up In Embryo, Down In Start, Up In P0, Up In P3, Up In 24Hr, Down In FiveDay
###Second: Find Features Down In Embryo, Up In Start, Down In P0, Down In P3, Down In 24Hr, Up In FiveDay
###Third: Find Features Down In Embryo, Down In Start, Up In P0, Up In P3, Up In 24Hr, Up In FiveDay 
###Fourth: Find Features Up In Embryo, Up In Start, Down In P0, Down In P3, Down In 24 Hr, Down FiveDay

##List of needed features
#####
# StartingOverEmbryoGenes.cleaned.NoCC
# StartingOverP0Genes.cleaned.NoCC
# StartingOverP3Genes.cleaned.NoCC
# StartingOver24HrGenes.cleaned.NoCC
# StartingOverFiveDayGenes.cleaned.NoCC
# 
# FiveDayOverEmbryoGenes.cleaned.NoCC
# FiveDayOverP0Genes.cleaned.NoCC
# FiveDayOverP3Genes.cleaned.NoCC
# FiveDayOver24HrGenes.cleaned.NoCC
# FiveDayOverStartGenes.cleaned.NoCC
# 
# EmbryoOverFiveDayGenes.cleaned.NoCC
# EmbryoOverP0Genes.cleaned.NoCC
# EmbryoOverP3Genes.cleaned.NoCC
# EmbryoOver24HrGenes.cleaned.NoCC
# EmbryoOverStartGenes.cleaned.NoCC
# 
# Eng24HrOverEmbryoGenes.cleaned.NoCC
# Eng24HrOverStartGenes.cleaned.NoCC
# Eng24HrOverP0Genes.cleaned.NoCC
# Eng24HrOverP3Genes.cleaned.NoCC
# Eng24HrOverFiveDayGenes.cleaned.NoCC
# 
# P0OverEmbryoGenes.cleaned.NoCC
# P0OverStartGenes.cleaned.NoCC
# P0OverP3Genes.cleaned.NoCC
# P0Over24HrGenes.cleaned.NoCC
# P0OverFiveDayGenes.cleaned.NoCC
# 
# P3OverEmbryoGenes.cleaned.NoCC
# P3OverStartGenes.cleaned.NoCC
# P3OverP0Genes.cleaned.NoCC
# P3Over24HrGenes.cleaned.NoCC
# P3OverFiveDayGenes.cleaned.NoCC
#####

###First Heatmap
#####
#We need: EmbryoOverStart, P0OverStart, P3OverStart, 24Hr Over Start, 
#         EmbryoOverFiveDay, P0OverFiveDay, P3OverFiveDay, 24HrOverFiveDay

a1 <- intersect(EmbryoOverStartGenes.cleaned.NoCC, P0OverStartGenes.cleaned.NoCC)
a2<- intersect(a1, P3OverStartGenes.cleaned.NoCC)
a3<- intersect(a2, Eng24HrOverStartGenes.cleaned.NoCC)

b1 <- intersect(EmbryoOverFiveDayGenes.cleaned.NoCC, P0OverFiveDayGenes.cleaned.NoCC)
b2 <- intersect(b1, P3OverFiveDayGenes.cleaned.NoCC)
b3 <- intersect(b2, Eng24HrOverFiveDayGenes.cleaned.NoCC)
FirstHeatmap <- intersect(a3, b3)
FirstHeatmap2 <-intersect(a3, b2)
#####

###Second Heatmap
#####
#We need: StartOverEmbryo, StartOverP0, StartOverP3, StartOver24Hr
#         FiveDayOverEmbryo, FiveDayOverP0, FiveDayOverP3, FiveDayOver24Hr
c1 <- intersect(StartingOverEmbryoGenes.cleaned.NoCC, StartingOverP0Genes.cleaned.NoCC)
c2<- intersect(c1, StartingOverP3Genes.cleaned.NoCC)
c3<- intersect(c2, StartingOver24HrGenes.cleaned.NoCC)

d1 <- intersect(FiveDayOverEmbryoGenes.cleaned.NoCC, FiveDayOverP0Genes.cleaned.NoCC)
d2 <- intersect(d1, FiveDayOverP3Genes.cleaned.NoCC)
d3 <- intersect(d2, FiveDayOver24HrGenes.cleaned.NoCC)
SecondHeatmap <- intersect(c3, d3)
#####

###Third Heatmap
#####
#We need: P0OverEmbryo, P3OverEmbryo, 24HrOverEmbryo, FiveDayOverEmbryo
#         P0OverStart, P3OverStart, 24HrOverStart, FiveDayOverStart
e1 <- intersect(P0OverEmbryoGenes.cleaned.NoCC, P3OverEmbryoGenes.cleaned.NoCC)
e2<- intersect(e1, Eng24HrOverEmbryoGenes.cleaned.NoCC)
e3<- intersect(e2, FiveDayOverEmbryoGenes.cleaned.NoCC)

f1 <- intersect(P0OverStartGenes.cleaned.NoCC, P3OverStartGenes.cleaned.NoCC)
f2 <- intersect(f1, Eng24HrOverStartGenes.cleaned.NoCC)
f3 <- intersect(f2, FiveDayOverStartGenes.cleaned.NoCC)
ThirdHeatmap <- intersect(e3, f3)
#####

###Fourth Heatmap
#####
#We need: EmbryoOverP0, EmbryoOverP3, EmbryoOver24Hr, EmbryoOverFiveDay
#         StartOverP0, StartOverP3, StartOver24Hr, StartOverFiveDay
g1 <- intersect(EmbryoOverP0Genes.cleaned.NoCC, EmbryoOverP3Genes.cleaned.NoCC)
g2<- intersect(g1, EmbryoOver24HrGenes.cleaned.NoCC)
g3<- intersect(g2, EmbryoOverFiveDayGenes.cleaned.NoCC)

h1 <- intersect(StartingOverP0Genes.cleaned.NoCC, StartingOverP3Genes.cleaned.NoCC)
h2 <- intersect(h1, StartingOver24HrGenes.cleaned.NoCC)
h3 <- intersect(h2, StartingOverFiveDayGenes.cleaned.NoCC)
FourthHeatmap <- intersect(g3, h3)

#####
setwd("/gpfs/gibbs/project/kaminski/hk738/Storage/Satoshi_Data/BASC_Analysis/BASC_Combined/2025.1.30/Heatmaps_(2000Feats)/Sam_Requested_Heatmap")

#FirstHeatmap
png(file = 'FirstHeatmap.1.Heatmap.png', width=16, height=12, units = 'in',res=300, type = "cairo")
CustomHeatmapOnly3(object = downsampled,
                   data.type = 'RNA', # sets the assay to pull data from
                   primary = 'ConditionRearranged' , # sets the first metadata slot to use (primary column grouping)
                   secondary = 'SampleRearranged' , # sets the second metadata slot to use (secondary column grouping)
                   #tertiary = 'Sample' , # sets the third metadata slot to use (tertiary column grouping)
                   #quarternary = 'orig.ident' , # not used in this specific function
                   primary.cols = col.pal$Condition, # Needs to be a named list of colors
                   secondary.cols = col.pal$SampleRearranged, # Needs to be a named list of colors
                   #tertiary.cols = col.pal$SampleRearranged, # Needs to be a named list of colors
                   #quarternary.cols = NULL, # not used in this specific function
                   features = FirstHeatmap[1:50], # defined above, customize
                   labels = c('Condition', 'Sample'), # Needs to ordered IDENTICALLY to primary, secondary, tertiary, above, or will produce inaccurately labeled plot
                   selected.row.anotations=FirstHeatmap[1:50],
                   selected.label.size = 10,
                   use.scale.data = T,
                   range.frac = 0.35)
dev.off()

png(file = 'FirstHeatmap.2.Heatmap.png', width=16, height=12, units = 'in',res=300, type = "cairo")
CustomHeatmapOnly3(object = downsampled,
                   data.type = 'RNA', # sets the assay to pull data from
                   primary = 'ConditionRearranged' , # sets the first metadata slot to use (primary column grouping)
                   secondary = 'SampleRearranged' , # sets the second metadata slot to use (secondary column grouping)
                   #tertiary = 'Sample' , # sets the third metadata slot to use (tertiary column grouping)
                   #quarternary = 'orig.ident' , # not used in this specific function
                   primary.cols = col.pal$Condition, # Needs to be a named list of colors
                   secondary.cols = col.pal$SampleRearranged, # Needs to be a named list of colors
                   #tertiary.cols = col.pal$SampleRearranged, # Needs to be a named list of colors
                   #quarternary.cols = NULL, # not used in this specific function
                   features = FirstHeatmap[51:93], # defined above, customize
                   labels = c('Condition', 'Sample'), # Needs to ordered IDENTICALLY to primary, secondary, tertiary, above, or will produce inaccurately labeled plot
                   selected.row.anotations=FirstHeatmap[51:93],
                   selected.label.size = 10,
                   use.scale.data = T,
                   range.frac = 0.35)
dev.off()
#SecondHeatmap
png(file = 'SecondHeatmap.1.Heatmap.png', width=16, height=12, units = 'in',res=300, type = "cairo")
CustomHeatmapOnly3(object = downsampled,
                   data.type = 'RNA', # sets the assay to pull data from
                   primary = 'ConditionRearranged' , # sets the first metadata slot to use (primary column grouping)
                   secondary = 'SampleRearranged' , # sets the second metadata slot to use (secondary column grouping)
                   #tertiary = 'Sample' , # sets the third metadata slot to use (tertiary column grouping)
                   #quarternary = 'orig.ident' , # not used in this specific function
                   primary.cols = col.pal$Condition, # Needs to be a named list of colors
                   secondary.cols = col.pal$SampleRearranged, # Needs to be a named list of colors
                   #tertiary.cols = col.pal$SampleRearranged, # Needs to be a named list of colors
                   #quarternary.cols = NULL, # not used in this specific function
                   features = SecondHeatmap[1:33], # defined above, customize
                   labels = c('Condition', 'Sample'), # Needs to ordered IDENTICALLY to primary, secondary, tertiary, above, or will produce inaccurately labeled plot
                   selected.row.anotations=SecondHeatmap[1:33],
                   selected.label.size = 10,
                   use.scale.data = T,
                   range.frac = 0.35)
dev.off()
#ThirdHeatmap
png(file = 'ThirdHeatmap.1.Heatmap.png', width=16, height=12, units = 'in',res=300, type = "cairo")
CustomHeatmapOnly3(object = downsampled,
                   data.type = 'RNA', # sets the assay to pull data from
                   primary = 'ConditionRearranged' , # sets the first metadata slot to use (primary column grouping)
                   secondary = 'SampleRearranged' , # sets the second metadata slot to use (secondary column grouping)
                   #tertiary = 'Sample' , # sets the third metadata slot to use (tertiary column grouping)
                   #quarternary = 'orig.ident' , # not used in this specific function
                   primary.cols = col.pal$Condition, # Needs to be a named list of colors
                   secondary.cols = col.pal$SampleRearranged, # Needs to be a named list of colors
                   #tertiary.cols = col.pal$SampleRearranged, # Needs to be a named list of colors
                   #quarternary.cols = NULL, # not used in this specific function
                   features = ThirdHeatmap[1:50], # defined above, customize
                   labels = c('Condition', 'Sample'), # Needs to ordered IDENTICALLY to primary, secondary, tertiary, above, or will produce inaccurately labeled plot
                   selected.row.anotations=ThirdHeatmap[1:50],
                   selected.label.size = 10,
                   use.scale.data = T,
                   range.frac = 0.35)
dev.off()
png(file = 'ThirdHeatmap.2.Heatmap.png', width=16, height=12, units = 'in',res=300, type = "cairo")
CustomHeatmapOnly3(object = downsampled,
                   data.type = 'RNA', # sets the assay to pull data from
                   primary = 'ConditionRearranged' , # sets the first metadata slot to use (primary column grouping)
                   secondary = 'SampleRearranged' , # sets the second metadata slot to use (secondary column grouping)
                   #tertiary = 'Sample' , # sets the third metadata slot to use (tertiary column grouping)
                   #quarternary = 'orig.ident' , # not used in this specific function
                   primary.cols = col.pal$Condition, # Needs to be a named list of colors
                   secondary.cols = col.pal$SampleRearranged, # Needs to be a named list of colors
                   #tertiary.cols = col.pal$SampleRearranged, # Needs to be a named list of colors
                   #quarternary.cols = NULL, # not used in this specific function
                   features = ThirdHeatmap[51:100], # defined above, customize
                   labels = c('Condition', 'Sample'), # Needs to ordered IDENTICALLY to primary, secondary, tertiary, above, or will produce inaccurately labeled plot
                   selected.row.anotations=ThirdHeatmap[51:100],
                   selected.label.size = 10,
                   use.scale.data = T,
                   range.frac = 0.35)
dev.off()

png(file = 'ThirdHeatmap.3.Heatmap.png', width=16, height=12, units = 'in',res=300, type = "cairo")
CustomHeatmapOnly3(object = downsampled,
                   data.type = 'RNA', # sets the assay to pull data from
                   primary = 'ConditionRearranged' , # sets the first metadata slot to use (primary column grouping)
                   secondary = 'SampleRearranged' , # sets the second metadata slot to use (secondary column grouping)
                   #tertiary = 'Sample' , # sets the third metadata slot to use (tertiary column grouping)
                   #quarternary = 'orig.ident' , # not used in this specific function
                   primary.cols = col.pal$Condition, # Needs to be a named list of colors
                   secondary.cols = col.pal$SampleRearranged, # Needs to be a named list of colors
                   #tertiary.cols = col.pal$SampleRearranged, # Needs to be a named list of colors
                   #quarternary.cols = NULL, # not used in this specific function
                   features = ThirdHeatmap[101:150], # defined above, customize
                   labels = c('Condition', 'Sample'), # Needs to ordered IDENTICALLY to primary, secondary, tertiary, above, or will produce inaccurately labeled plot
                   selected.row.anotations=ThirdHeatmap[101:150],
                   selected.label.size = 10,
                   use.scale.data = T,
                   range.frac = 0.35)
dev.off()

png(file = 'ThirdHeatmap.4.Heatmap.png', width=16, height=12, units = 'in',res=300, type = "cairo")
CustomHeatmapOnly3(object = downsampled,
                   data.type = 'RNA', # sets the assay to pull data from
                   primary = 'ConditionRearranged' , # sets the first metadata slot to use (primary column grouping)
                   secondary = 'SampleRearranged' , # sets the second metadata slot to use (secondary column grouping)
                   #tertiary = 'Sample' , # sets the third metadata slot to use (tertiary column grouping)
                   #quarternary = 'orig.ident' , # not used in this specific function
                   primary.cols = col.pal$Condition, # Needs to be a named list of colors
                   secondary.cols = col.pal$SampleRearranged, # Needs to be a named list of colors
                   #tertiary.cols = col.pal$SampleRearranged, # Needs to be a named list of colors
                   #quarternary.cols = NULL, # not used in this specific function
                   features = ThirdHeatmap[151:191], # defined above, customize
                   labels = c('Condition', 'Sample'), # Needs to ordered IDENTICALLY to primary, secondary, tertiary, above, or will produce inaccurately labeled plot
                   selected.row.anotations=ThirdHeatmap[151:191],
                   selected.label.size = 10,
                   use.scale.data = T,
                   range.frac = 0.35)
dev.off()

#FourthHeatmap

png(file = 'FourthHeatmap.1.Heatmap.png', width=16, height=12, units = 'in',res=300, type = "cairo")
CustomHeatmapOnly3(object = downsampled,
                   data.type = 'RNA', # sets the assay to pull data from
                   primary = 'ConditionRearranged' , # sets the first metadata slot to use (primary column grouping)
                   secondary = 'SampleRearranged' , # sets the second metadata slot to use (secondary column grouping)
                   #tertiary = 'Sample' , # sets the third metadata slot to use (tertiary column grouping)
                   #quarternary = 'orig.ident' , # not used in this specific function
                   primary.cols = col.pal$Condition, # Needs to be a named list of colors
                   secondary.cols = col.pal$SampleRearranged, # Needs to be a named list of colors
                   #tertiary.cols = col.pal$SampleRearranged, # Needs to be a named list of colors
                   #quarternary.cols = NULL, # not used in this specific function
                   features = FourthHeatmap[1:50], # defined above, customize
                   labels = c('Condition', 'Sample'), # Needs to ordered IDENTICALLY to primary, secondary, tertiary, above, or will produce inaccurately labeled plot
                   selected.row.anotations=FourthHeatmap[1:50],
                   selected.label.size = 10,
                   use.scale.data = T,
                   range.frac = 0.35)
dev.off()

png(file = 'FourthHeatmap.2.Heatmap.png', width=16, height=12, units = 'in',res=300, type = "cairo")
CustomHeatmapOnly3(object = downsampled,
                   data.type = 'RNA', # sets the assay to pull data from
                   primary = 'ConditionRearranged' , # sets the first metadata slot to use (primary column grouping)
                   secondary = 'SampleRearranged' , # sets the second metadata slot to use (secondary column grouping)
                   #tertiary = 'Sample' , # sets the third metadata slot to use (tertiary column grouping)
                   #quarternary = 'orig.ident' , # not used in this specific function
                   primary.cols = col.pal$Condition, # Needs to be a named list of colors
                   secondary.cols = col.pal$SampleRearranged, # Needs to be a named list of colors
                   #tertiary.cols = col.pal$SampleRearranged, # Needs to be a named list of colors
                   #quarternary.cols = NULL, # not used in this specific function
                   features = FourthHeatmap[51:100], # defined above, customize
                   labels = c('Condition', 'Sample'), # Needs to ordered IDENTICALLY to primary, secondary, tertiary, above, or will produce inaccurately labeled plot
                   selected.row.anotations=FourthHeatmap[51:100],
                   selected.label.size = 10,
                   use.scale.data = T,
                   range.frac = 0.35)
dev.off()

png(file = 'FourthHeatmap.3.Heatmap.png', width=16, height=12, units = 'in',res=300, type = "cairo")
CustomHeatmapOnly3(object = downsampled,
                   data.type = 'RNA', # sets the assay to pull data from
                   primary = 'ConditionRearranged' , # sets the first metadata slot to use (primary column grouping)
                   secondary = 'SampleRearranged' , # sets the second metadata slot to use (secondary column grouping)
                   #tertiary = 'Sample' , # sets the third metadata slot to use (tertiary column grouping)
                   #quarternary = 'orig.ident' , # not used in this specific function
                   primary.cols = col.pal$Condition, # Needs to be a named list of colors
                   secondary.cols = col.pal$SampleRearranged, # Needs to be a named list of colors
                   #tertiary.cols = col.pal$SampleRearranged, # Needs to be a named list of colors
                   #quarternary.cols = NULL, # not used in this specific function
                   features = FourthHeatmap[101:150], # defined above, customize
                   labels = c('Condition', 'Sample'), # Needs to ordered IDENTICALLY to primary, secondary, tertiary, above, or will produce inaccurately labeled plot
                   selected.row.anotations=FourthHeatmap[101:150],
                   selected.label.size = 10,
                   use.scale.data = T,
                   range.frac = 0.35)
dev.off()

png(file = 'FourthHeatmap.4.Heatmap.png', width=16, height=12, units = 'in',res=300, type = "cairo")
CustomHeatmapOnly3(object = downsampled,
                   data.type = 'RNA', # sets the assay to pull data from
                   primary = 'ConditionRearranged' , # sets the first metadata slot to use (primary column grouping)
                   secondary = 'SampleRearranged' , # sets the second metadata slot to use (secondary column grouping)
                   #tertiary = 'Sample' , # sets the third metadata slot to use (tertiary column grouping)
                   #quarternary = 'orig.ident' , # not used in this specific function
                   primary.cols = col.pal$Condition, # Needs to be a named list of colors
                   secondary.cols = col.pal$SampleRearranged, # Needs to be a named list of colors
                   #tertiary.cols = col.pal$SampleRearranged, # Needs to be a named list of colors
                   #quarternary.cols = NULL, # not used in this specific function
                   features = FourthHeatmap[151:181], # defined above, customize
                   labels = c('Condition', 'Sample'), # Needs to ordered IDENTICALLY to primary, secondary, tertiary, above, or will produce inaccurately labeled plot
                   selected.row.anotations=FourthHeatmap[151:181],
                   selected.label.size = 10,
                   use.scale.data = T,
                   range.frac = 0.35)
dev.off()
