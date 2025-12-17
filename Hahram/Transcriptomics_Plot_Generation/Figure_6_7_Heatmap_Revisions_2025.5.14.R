Figure_6_All_Genes <- c("Cxcl2", "Cxcl3", "Cxcl6", "Cxcl1", "Il1a",	"Ereg", "Areg", "Hbegf", "Asf1b",	"Cish", "Selenoi", "Mastl", "Ska3", "Trip13", "Rrm2",	"Gen1",	"Spc25", "Pbk", "Ptgs2", "Angptl4", "Adtrp", "Pde4c.1", "Odc1", "Ptgr1", "Lyar", "Depdc1b",	  "Zc3h12a", "Cnfn", "Melk", "Acot7", "Nfkbia", "Tnfaip3", "Mms22l", "Fen1", "Ckap2l", "Ska1", "Rpa2",	"Arhgap11a",	"Ncapg2", "Neto2", "Exo1", "Trmt61a", "Taf4b", "Atf3", "Mybl2",	"Gins1",	"Mcm5",	"Ube2c", "Prc1", "Igf1r", "Ldlr", "Lrp8", "F2rl2", "Aurka", "Plaur", "Tnfrsf12a", "Itga6",	"Hras", "Il24", "Lamb3",	"Nrg1", "Rtn4", "Vasp",	"Nptn")

Figure_6_All_Genes_Select <- c("Cxcl2", "Cxcl3", "Cxcl6", "Cxcl1", "Il1a", "Ereg", "Areg", "Hbegf", "Ptgs2", "Angptl4", "Adtrp", "Pde4c.1", "Odc1", "Ptgr1", "Lyar", "Depdc1b", "Zc3h12a", "Cnfn", "Melk", "Acot7", "Nfkbia",	  "Tnfaip3", "Taf4b", "Atf3", "Mybl2", "Igf1r", "Ldlr", "Lrp8", "F2rl2", "Itga6", "Nptn"
)

Fig6.diff <- setdiff(Figure_6_All_Genes, Figure_6_All_Genes_Select)
downsampled.2dto3d <- RenameIdents(downsampled.2dto3d,
                                   "Engineered-24Hr" = "24 Hour")
unique(Idents(downsampled.2dto3d))
downsampled.2dto3d$ConditionRearranged <- Idents(downsampled.2dto3d)

downsampled.2dto3d$ConditionRearranged <- factor(downsampled.2dto3d$ConditionRearranged, levels = c("P0-2D", "P3-2D", "24 Hour"))

png(file = 'Figure_6_Heatmap.png', width=8.5, height=6, units = 'in',res=300, type = "cairo")

CustomHeatmapOnly3(object = downsampled.2dto3d,
                   data.type = 'RNA', # sets the assay to pull data from
                   primary = 'ConditionRearranged' , # sets the first metadata slot to use (primary column grouping)
                   secondary = 'SampleRearranged' , # sets the second metadata slot to use (secondary column grouping)
                   #tertiary = 'Sample' , # sets the third metadata slot to use (tertiary column grouping)
                   #quarternary = 'orig.ident' , # not used in this specific function
                   primary.cols = c("green", "yellow", "purple"), # Needs to be a named list of colors
                   secondary.cols = col.pal$SampleRearranged, # Needs to be a named list of colors
                   #tertiary.cols = col.pal$SampleRearranged, # Needs to be a named list of colors
                   #quarternary.cols = NULL, # not used in this specific function
                   features = Figure_6_All_Genes, # defined above, customize
                   labels = c('Condition','Sample'), # Needs to ordered IDENTICALLY to primary, secondary, tertiary, above, or will produce inaccurately labeled plot
                   selected.row.anotations=Figure_6_All_Genes_Select,
                   selected.label.size = 9,
                   use.scale.data = T,
                   range.frac = 0.2)
dev.off()

Figure_7_All_Genes<- c("Ebf1",	"Bcl6",	"Bhlhe40",	"Litaf",	"Camta1",	"Zbtb20",	"Runx1",	"Epas1",	"Elf3",	"Foxo1",	"Serpine1",	"Acap1",	"Ldlrad4",	"Pmepa1",	"Entpd4", "Ppp1r14a",	"Sema4g",	"Itpr1",	"Cxcr4",	"Ackr3",	"Plxna2",	
                       "Bnip3",	"Ltbp3",	"Fzd2",	"Ephb2",	"Barx2",	"Parm1",	"Car9",	"Tmprss2",	"Cmbl", "Wfdc1",	"Mmp7",	"Lamb2",	"Tgfbi",	"Ankub1",	"Zc3h6",	
                       "Rhobtb3",		"Vill",	"Scgb1a1",	"Rtn4rl1",	"Bmp6",	"Vegfa",	"Efnb2", "Mfge8",	"Vtcn1",	"Igfbp3",	
                       "Dlk2",	"Fosl2",	"Rere",		
                       "Tsc22d3",	"Aff1", 	"Ahr", 	"Itgb6",	"Itgav",	"Cux1",	"Dbp",	"Nfat5",	"Atrn",	"F3","Sdc4",	"Tspan1",	"Stra6",	"St6gal1",	"Camk2n1",	
                       "Smad7",	"Dlc1",	"Wwc2",	"Zfp36l2",	"Vim",	"Pik3cb", "Prom1",	"Krt7",	"Pik3c2a",	"Furin",	"Sulf1", "Aqp5", "Pdpn")


Figure_7_All_Genes_Select <- c("Ebf1",	"Bcl6",	"Bhlhe40",	"Litaf",	"Camta1",	"Zbtb20",	"Runx1",	"Epas1",	"Elf3",	"Foxo1","Cxcr4",	"Ackr3",	"Plxna2", "Bnip3",	"Ltbp3",	"Fzd2",	"Ephb2", "Wfdc1",	"Mmp7",	"Lamb2",	"Tgfbi", "Bmp6",	"Vegfa",	"Efnb2",	"Mfge8",	"Vtcn1",	"Igfbp3",	"Dlk2",	"Itgb6", "Itgav", "St6gal1",	"Camk2n1", "Smad7",	"Dlc1",	"Wwc2",	"Zfp36l2",	"Vim",	"Pik3cb", "Prom1",	"Krt7",	"Pik3c2a",	"Furin",	"Sulf1", "Aqp5", "Pdpn")

Idents(downsampled.3dto5d) <- downsampled.3dto5d$ConditionRearranged
downsampled.3dto5d <- RenameIdents(downsampled.3dto5d,
                                   "Engineered-24Hr" = "24 Hour",
                                   "Engineered-FiveDay" = "Five Day"
)
downsampled.3dto5d$ConditionRearranged <- Idents(downsampled.3dto5d)
downsampled.3dto5d$ConditionRearranged <- factor(downsampled.3dto5d$ConditionRearranged, levels = c("P3-2D", "24 Hour", "Five Day"))
png(file = 'Figure_7_Heatmap.png', width=8.5, height=6, units = 'in',res=300, type = "cairo")
CustomHeatmapOnly3(object = downsampled.3dto5d,
                   data.type = 'RNA', # sets the assay to pull data from
                   primary = 'ConditionRearranged' , # sets the first metadata slot to use (primary column grouping)
                   secondary = 'SampleRearranged' , # sets the second metadata slot to use (secondary column grouping)
                   #tertiary = 'Sample' , # sets the third metadata slot to use (tertiary column grouping)
                   #quarternary = 'orig.ident' , # not used in this specific function
                   primary.cols = c("yellow", "purple", "orange"), # Needs to be a named list of colors
                   secondary.cols = col.pal$SampleRearranged, # Needs to be a named list of colors
                   #tertiary.cols = col.pal$SampleRearranged, # Needs to be a named list of colors
                   #quarternary.cols = NULL, # not used in this specific function
                   features = Figure_7_All_Genes, # defined above, customize
                   labels = c('Condition','Sample'), # Needs to ordered IDENTICALLY to primary, secondary, tertiary, above, or will produce inaccurately labeled plot
                   selected.row.anotations=Figure_7_All_Genes_Select,
                   selected.label.size = 8,
                   use.scale.data = T,
                   range.frac = 0.2)
dev.off()
