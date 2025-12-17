setwd("/gpfs/gibbs/project/kaminski/hk738/Storage/Satoshi_Data/BASC_Analysis/BASC_Combined/2025.5.19")
####Figure_8 Update
load("/gpfs/gibbs/project/kaminski/hk738/Storage/Satoshi_Data/BASC_Analysis/BASC_Combined/2025.2.4/AllEpiObjectMerge.Updated.2025.3.2.Robj")

table(merged.obj$Condition)
Idents(merged.obj) <- merged.obj$Condition
downsampled <- subset(merged.obj, downsample = 2315)
downsampled <- ScaleData(downsampled,features = rownames(downsampled))

Idents(downsampled) <- downsampled$Condition
downsampled <- RenameIdents(downsampled,
                            "Embryogenesis" = "Embryo",
                            "Engineered-24Hr" = "24 Hour",
                            "Engineered-FiveDay" = "Five Day")
downsampled$ConditionRename <- Idents(downsampled)
downsampled$ConditionRename <- factor(downsampled$ConditionRename, levels = c("Embryo", "St. Dclk1", "P0-2D", "P3-2D", "24 Hour", "Five Day"))

genes_from_first <- figure_8
genes_from_second<-c( "Dclk1", "Pdpn", "Aqp5", "Clic5", "Akap5", "Cxcl2",  "Il1a", "Cxcl3", "Cxcl1", "Cxcl6", "Lamc2", "Ptgs2", "Areg", "Ereg",  "Krt5",  "Cdkn2a", "Il24", "Foxq1", "Wnt4", "Dsg2", "Lamb3", "Krt17", "Nrg1", "Jag2", "Vtcn1",  "Mmp7", "Ebf1", "Ephb2"
)


common_features <- intersect(genes_from_first, genes_from_second)


#"Ereg"  "Areg"  "Ptgs2" "Cxcl1" "Mmp7"  "Ebf1" , “Mfge8”, “Efnb2”, “Krt5


fig_8_newlist<- c("Mdk", "Igf2", "Tead2", "Ror2", "Col2a1", "Emid1", "Zfp521", "Cbln2", "Cdkn1c", "Igfbp4", "Sox6", "Sox5", "Bmp5", "Fgf12",  "Muc5b", "Cfap99", "Cfap65", "Drc1", "Dcdc5", "Socs2", "Dmc1", "Cxcl3", "Scgb3a1", "Sftpd", "Sftpb", "Akap5", "Pdpn", "Aqp5", "Clic5", "Mt2A", "Slc7a5", "Cpe", "Plat", "Mt1", "Pamr1", "Sord", "Them6", "Aldoc",    "Chac1", "Il1a", "Tnfaip3", "Ereg", "Areg", "Ptgs2", "Plaur", "Fosl1", "Nppb", "Cxcl1", "Ripk3", "Col17a1",  "**Fst", "Ptgr1", "Aldh3a1","Cst6",  "Krt16", "Hbegf", "Sox15", "Tp63", "Snai2", "Serpinb2", "Wnt4", "Lamb3", "Krt17", "Cdkn2a", "Krt5", "Efnb2", "Mfge8", "Lamc2", "Itgb6", "Egfr",  "Cryba4", "Vim", "Serpine1", "Wfdc1", "Tgfbi", "Bmp6",  "Mmp7",  "Vtcn1", "Ephb2", "Ptprb", "Ebf1", "Birc3")

fig_8_foi <- c("Mdk", "Sox6", "Igf2", "Col2a1", "Sox5", "Bmp5", "Muc5b","Fgf12","Tead2", "Il1a", "Cxcl3", "Cxcl1", "Hbegf", "Sox15", "Tp63", "Snai2", "Egfr", "Wnt4", "Lamb3", "Krt17", "Vtcn1", "Ephb2", "Ereg", "Areg", "Mmp7", "Ebf1", "Ptgs2", "Pdpn", "Aqp5", "Clic5", "Akap5")

fig_8_foi_v2 <- c("Mdk", "Sox6", "Igf2", "Col2a1", "Emid1", "Tgfbp4", "Sox5", "Bmp5", "Fgf12", "Muc5b", "Cfap99", "Cfap65", "Drc1", "Tead2", "Il1a", "Cxcl3", "Cxcl1", "Scgb3a1", "Sftpd", "Sftpb", "Plat", "Plaur", "Fosl1", "Ripk3", "Col17a1", "Ptgr1", "Aldh3a1", "Cst6", "Krt16", "Serpinb2", "Tnfaip3", "Serpine1", "Wfdc1", "Tgfbi", "Bmp6", "Ptprb", "Ebf1", "Birc3", "Hbegf", "Sox15", "Tp63", "Snai2", "Egfr", "Wnt4", "Lamb3", "Krt17", "Vtcn1", "Ephb2", "Ereg", "Areg", "Mmp7", "Ebf1", "Ptgs2", "Pdpn", "Aqp5", "Clic5", "Akap5", "Vim", "Lamc2", "Cdkn2a", "Itgb6","Mfge8", "Efnb2", "Krt5")

abBasal.foi <- c("Vim", "Lamc2", "Mmp7", "Cdkn2a", "Itgb6", "Ephb2", "Mfge8", "Efnb2", "Krt5", "Runx1")


VlnPlot(downsampled, features = "Hopx", group.by = "ConditionRename")

png(file = 'Figure8.Heatmap.Updated.2025.5.29.png', width=15.5, height=13, units = 'in',res=300, type = "cairo")
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
                   features = fig_8_newlist, # defined above, customize
                   labels = c('Condition', 'Sample'), # Needs to ordered IDENTICALLY to primary, secondary, tertiary, above, or will produce inaccurately labeled plot
                   selected.row.anotations=fig_8_foi_v2,
                   selected.label.size = 12,
                   use.scale.data = T,
                   range.frac = 0.25)
dev.off()

png(file = 'Figure8.Heatmap.Updated.Full.2025.5.19.png', width=16, height=12, units = 'in',res=300, type = "cairo")
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
                   features = abBasal.foi, # defined above, customize
                   labels = c('Condition', 'Sample'), # Needs to ordered IDENTICALLY to primary, secondary, tertiary, above, or will produce inaccurately labeled plot
                   selected.row.anotations=abBasal.foi,
                   selected.label.size = 10,
                   use.scale.data = T,
                   range.frac = 0.25)
dev.off()
