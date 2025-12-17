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



##EMBRYO AND P0, Over Start, P3, BDL, FiveDay
a1 <- intersect(EmbryoOverStartGenes.cleaned.NoCC, EmbryoOverP3Genes.cleaned.NoCC)
a2 <- intersect(a1, EmbryoOver24HrGenes.cleaned.NoCC)
a3 <- intersect(a2, EmbryoOverFiveDayGenes.cleaned.NoCC)

b1 <- intersect(P0OverStartGenes.cleaned.NoCC, P0OverP3Genes.cleaned.NoCC)
b2 <- intersect(b1, P0Over24HrGenes.cleaned.NoCC)
b3 <- intersect(b2, P0OverFiveDayGenes.cleaned.NoCC)
EmbryoP0.OverRest <- intersect(a3, b3)


##EMBRYO AND P3, over Start, P0, BDL, FiveDay
c1 <- intersect(EmbryoOverStartGenes.cleaned.NoCC, EmbryoOverP0Genes.cleaned.NoCC)
c2 <- intersect(c1, EmbryoOver24HrGenes.cleaned.NoCC)
c3 <- intersect(c2, EmbryoOverFiveDayGenes.cleaned.NoCC)

d1 <- intersect(P3OverStartGenes.cleaned.NoCC, P3OverP0Genes.cleaned.NoCC)
d2 <- intersect(d1, P3Over24HrGenes.cleaned.NoCC)
d3 <- intersect(d2, P3OverFiveDayGenes.cleaned.NoCC)
EmbryoP3.OverRest <- intersect(c3, d3)


##EMBRYO AND BDL, over Start, P0, P3, FiveDaY
e1 <- intersect(EmbryoOverStartGenes.cleaned.NoCC, EmbryoOverP0Genes.cleaned.NoCC)
e2 <- intersect(e1, EmbryoOverP3Genes.cleaned.NoCC)
e3 <- intersect(e2, EmbryoOverFiveDayGenes.cleaned.NoCC)

f1 <- intersect(Eng24HrOverStartGenes.cleaned.NoCC, Eng24HrOverP0Genes.cleaned.NoCC)
f2 <- intersect(f1, Eng24HrOverP3Genes.cleaned.NoCC)
f3 <- intersect(f2, Eng24HrOverFiveDayGenes.cleaned.NoCC)
Embryo24Hr.OverRest <- intersect(e3, f3)

##EMBRYO AND FiveDay, over Start, P0, P3, BDL
g1 <- intersect(EmbryoOverStartGenes.cleaned.NoCC, EmbryoOverP0Genes.cleaned.NoCC)
g2 <- intersect(e1, EmbryoOverP3Genes.cleaned.NoCC)
g3 <- intersect(e2, EmbryoOver24HrGenes.cleaned.NoCC)

h1 <- intersect(FiveDayOverStartGenes.cleaned.NoCC, FiveDayOverP0Genes.cleaned.NoCC)
h2 <- intersect(f1, FiveDayOverP3Genes.cleaned.NoCC)
h3 <- intersect(f2, FiveDayOver24HrGenes.cleaned.NoCC)
EmbryoFiveDay.OverRest <- intersect(g3, h3)



##Start and FiveDay, Over Embryo, P0
m1 <- intersect(StartingOverEmbryoGenes.cleaned.NoCC, StartingOverP0Genes.cleaned.NoCC)
n1 <- intersect(FiveDayOverStartGenes.cleaned.NoCC, FiveDayOverP0Genes.cleaned.NoCC)

StartFiveDay.OverEmbryoP0 <- intersect(m1, n1)

##Start and FiveDay, Over Embryo, P3
o1 <- intersect(StartingOverEmbryoGenes.cleaned.NoCC, StartingOverP3Genes.cleaned.NoCC)
p1 <- intersect(FiveDayOverStartGenes.cleaned.NoCC, FiveDayOverP3Genes.cleaned.NoCC)

StartFiveDay.OverEmbryoP3 <- intersect(o1, p1)

##Start and FiveDay, Over Embryo, BDL
q1 <- intersect(StartingOverEmbryoGenes.cleaned.NoCC, StartingOver24HrGenes.cleaned.NoCC)
r1 <- intersect(FiveDayOverStartGenes.cleaned.NoCC, FiveDayOver24HrGenes.cleaned.NoCC) 

StartFiveDay.OverEmbryo24Hr <- intersect(q1, r1)

###########################Heatmaps


png(file = 'EmbryoP0.OverRest.Heatmap.png', width=16, height=12, units = 'in',res=300, type = "cairo")
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
                   features = EmbryoP0.OverRest[1:29], # defined above, customize
                   labels = c('Condition', 'Sample'), # Needs to ordered IDENTICALLY to primary, secondary, tertiary, above, or will produce inaccurately labeled plot
                   selected.row.anotations=EmbryoP0.OverRest[1:29],
                   selected.label.size = 10,
                   use.scale.data = T,
                   range.frac = 0.35)
dev.off()


png(file = 'EmbryoP3.OverRest.Heatmap.png', width=16, height=12, units = 'in',res=300, type = "cairo")
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
                   features = EmbryoP3.OverRest, # defined above, customize
                   labels = c('Condition', 'Sample'), # Needs to ordered IDENTICALLY to primary, secondary, tertiary, above, or will produce inaccurately labeled plot
                   selected.row.anotations=EmbryoP3.OverRest,
                   selected.label.size = 10,
                   use.scale.data = T,
                   range.frac = 0.35)
dev.off()

png(file = 'Embryo24Hr.OverRest.1.Heatmap.png', width=16, height=12, units = 'in',res=300, type = "cairo")
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
                   features = Embryo24Hr.OverRest[1:50], # defined above, customize
                   labels = c('Condition', 'Sample'), # Needs to ordered IDENTICALLY to primary, secondary, tertiary, above, or will produce inaccurately labeled plot
                   selected.row.anotations=Embryo24Hr.OverRest[1:50],
                   selected.label.size = 10,
                   use.scale.data = T,
                   range.frac = 0.35)
dev.off()

png(file = 'Embryo24Hr.OverRest.2.Heatmap.png', width=16, height=12, units = 'in',res=300, type = "cairo")
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
                   features = Embryo24Hr.OverRest[51:95], # defined above, customize
                   labels = c('Condition', 'Sample'), # Needs to ordered IDENTICALLY to primary, secondary, tertiary, above, or will produce inaccurately labeled plot
                   selected.row.anotations=Embryo24Hr.OverRest[51:95],
                   selected.label.size = 10,
                   use.scale.data = T,
                   range.frac = 0.35)
dev.off()


png(file = 'StartFiveDay.OverEmbryoP0.Heatmap.png', width=16, height=12, units = 'in',res=300, type = "cairo")
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
                   features = StartFiveDay.OverEmbryoP0, # defined above, customize
                   labels = c('Condition', 'Sample'), # Needs to ordered IDENTICALLY to primary, secondary, tertiary, above, or will produce inaccurately labeled plot
                   selected.row.anotations=StartFiveDay.OverEmbryoP0,
                   selected.label.size = 10,
                   use.scale.data = T,
                   range.frac = 0.35)
dev.off()

png(file = 'StartFiveDay.OverEmbryoP3.Heatmap.png', width=16, height=12, units = 'in',res=300, type = "cairo")
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
                   features = StartFiveDay.OverEmbryoP3, # defined above, customize
                   labels = c('Condition', 'Sample'), # Needs to ordered IDENTICALLY to primary, secondary, tertiary, above, or will produce inaccurately labeled plot
                   selected.row.anotations=StartFiveDay.OverEmbryoP3,
                   selected.label.size = 10,
                   use.scale.data = T,
                   range.frac = 0.35)
dev.off()


png(file = 'StartFiveDay.OverEmbryo24Hr.Heatmap.png', width=16, height=12, units = 'in',res=300, type = "cairo")
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
                   features = StartFiveDay.OverEmbryo24Hr, # defined above, customize
                   labels = c('Condition', 'Sample'), # Needs to ordered IDENTICALLY to primary, secondary, tertiary, above, or will produce inaccurately labeled plot
                   selected.row.anotations=StartFiveDay.OverEmbryo24Hr,
                   selected.label.size = 10,
                   use.scale.data = T,
                   range.frac = 0.35)
dev.off()