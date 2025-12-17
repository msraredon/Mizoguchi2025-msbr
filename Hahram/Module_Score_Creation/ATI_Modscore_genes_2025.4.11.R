
png(file = 'ATI_genes.png', width=16, height=12, units = 'in',res=300, type = "cairo")
CustomHeatmapOnly3(object = downsampled,
                   data.type = 'RNA', # sets the assay to pull data from
                   primary = 'ConditionRename' , # sets the first metadata slot to use (primary column grouping)
                   secondary = 'SampleRearranged' , # sets the second metadata slot to use (secondary column grouping)
                   #tertiary = 'Sample' , # sets the third metadata slot to use (tertiary column grouping)
                   #quarternary = 'orig.ident' , # not used in this specific function
                   primary.cols = col.pal$Condition, # Needs to be a named list of colors
                   secondary.cols = col.pal$SampleRearranged, # Needs to be a named list of colors
                   #tertiary.cols = col.pal$SampleRearranged, # Needs to be a named list of colors
                   #quarternary.cols = NULL, # not used in this specific function
                   features = ATI_genes.revised.2[1:50], # defined above, customize
                   labels = c('Condition', 'Sample'), # Needs to ordered IDENTICALLY to primary, secondary, tertiary, above, or will produce inaccurately labeled plot
                   selected.row.anotations=ATI_genes.revised.2[1:50],
                   selected.label.size = 10,
                   use.scale.data = T,
                   range.frac = 0.25)
dev.off()

png(file = 'ATI_genes.2.png', width=16, height=12, units = 'in',res=300, type = "cairo")
CustomHeatmapOnly3(object = downsampled,
                   data.type = 'RNA', # sets the assay to pull data from
                   primary = 'ConditionRename' , # sets the first metadata slot to use (primary column grouping)
                   secondary = 'SampleRearranged' , # sets the second metadata slot to use (secondary column grouping)
                   #tertiary = 'Sample' , # sets the third metadata slot to use (tertiary column grouping)
                   #quarternary = 'orig.ident' , # not used in this specific function
                   primary.cols = col.pal$Condition, # Needs to be a named list of colors
                   secondary.cols = col.pal$SampleRearranged, # Needs to be a named list of colors
                   #tertiary.cols = col.pal$SampleRearranged, # Needs to be a named list of colors
                   #quarternary.cols = NULL, # not used in this specific function
                   features = ATI_genes.revised.2[51:84], # defined above, customize
                   labels = c('Condition', 'Sample'), # Needs to ordered IDENTICALLY to primary, secondary, tertiary, above, or will produce inaccurately labeled plot
                   selected.row.anotations=ATI_genes.revised.2[51:84],
                   selected.label.size = 10,
                   use.scale.data = T,
                   range.frac = 0.25)
dev.off()
Idents(downsampled) <- downsampled$ConditionRename
VlnPlot(downsampled, features = ATI_genes.revised[41:84])

taking.out <- c("Cyp2b1", "Ccn2")
ATI_genes.revised.2 <- setdiff(ATI_genes.revised, taking.out)
ATI_genes.revised.2 <- c(ATI_genes.revised.2, "Hopx", "Igfbp2")
VlnPlot(downsampled, features = "Pdpn")



shared.genes.ATI <- intersect(Figure_7_All_Genes, ATI_genes.revised.2)
