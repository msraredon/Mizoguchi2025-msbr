Embryogenesis.remove <- c("Kcnc3", "Nnt", "Rbp1", "P3h4", "Ung", 
                          "Fgfrl1", "Bicc1", "Nek6", "Zfp428", "Dscc1",
                          "Fbln1", "P3h4", "Ctnnd2", "Mdfi", "H1fx", "Rasl11b",
                          "Fam111a", "Hba-a1", "Gpc6", "Dcakd", "Tmed1", "Crlf1", "Rps4y2",
                          "Etv5", "Pcpb4", "Pclaf", "Pcbd1", "Rfc3", "Pafah1b3", "Tmem9", "Hist1h1e",
                          "Tead2", "RGD1562378.1", "Acvr2b", "Tmem35b", "Bclaf3", "Fanci", "Tgif2",
                          "Thumpd2", "Zfp955a", "Ncs1", "Nsun6", "Mcm8", "Psrc1", "E2f7", "Mum1l1",
                          "Fiz1", "Pkdcc", "Haus5", "Logals1", "RGD1306227", "Twf2", "Hba-a2", "Nusap1", "Bend5",
                          "Osgep", "Skp2", "Blm", "Hdc", "Hif3a"
                          )

Starting.remove <- c("RGD156355", "5330417C22Rik", "Alas1", "Sftpb", "Pde4d", "Zfp36", "Pla2g16", "Npas2", "Rnasek",
                     "Cxcl17", "Ctsc", "Osbpl6", "AC128059.1", "Tc2n", "S100a1", "Pon3", "Wfdc2", "Scgb1a1", "Gda", "Rmdn2",
                     "Ank3", "Mt-nd4l", "Scp2", "St3gal1", "Vamp8", "Ica1", "RT1-S3", "Hist2h2aa3", "Crip1",
                     "Cyba", "Tyrobp", "Ptpn18", "Lgals3bp", "Tgoln2", "Shroom3", "Ilr1", "Atp2c2", "Ly6i", "Lrrc26",
                     "Irak2", "Pdpn", "Hspb1", "Epas1", "Ogfrl1", "Cebpd", "Acss1", "Ramp2", "Nebl", "Dnaja4", "Il13ra1", "Ilr1", "Sult1a1"
                     )
P0.remove <- c("Abo", "Akr7a3", "Ifi27l2b", "Hsd17b2", "Erg28", "Fuca2", "Nhp2")

FiveDay.remove <- c("Aff1", "RGD1562420", "Col8a1")
Heatmap1.remove <- c("Hist1h1b", "Pclaf", "Spc25", "Gar1", "Pa2g4","Hpf1", "Qdpr",
                     "Bud23", "Prim1", "Rrm1", "Apex1", "Fam83d", "Arhgap11a",
                     "Zfp57", "Mettl18"
                     )

Heatmap3.remove <- c( "Pinlyp", "Sncg", "Plaur", "Cxcl6", "Card19", "Car9", "Ptgr1",
                      "Iffo2", "Plat", "Adh7", "Dbi", "Ybx1-ps3",
                      "Tes", "Porcn", "Krt5", "Krt15", "Pgam1", 
                      "Akr1b1", "Med24", "Cald1", "Dusp5", "Sesn2",
                      "Serpinb2", "Rps6ka4", "Plcxd2", "Abcc3", "Rbp2",
                      "Mybl1", "Areg", "Cldn4", "Car13", "Fam83b", "Dennd2c"
                      )


Heatmap4.remove <- c("Mdk", "Pde3a", "St6galnac3", "Mia", "Adgrl3", "Foxp2", "H1fx", "Msi1",
                     "Cadm1", "Col9a2", "Hdc"," Srgap1", "Map3k21", "Rapgef5", "Nrp1", "Ralgapa1",
                     "Eln", "Napsa", "Plcb4", "Sftpa1", "Scl15a2", "Lrrk2", "Slc4a8", "Tcte1",
                     "Ntn1", "Slc34a2", "Pclo", "Mme", "Sidt1", "Kcne2", ""
)

Embryo.24Hr.remove <- c("Pclaf", "Hist1h2bk", "Mest", "Dtl", "Rfc2",
                        "Haus8", "Tcof1", "Spdl1", "Neto2", "Mcm10",
                        "Chaf1b", "Pole")
Embryo.P0.remove <- c("Tmem59l", "Tmem213", "Spats2l", "B4galnt3", "Cacna2d3", "Aif1l",
                      "Cmtm3", "Sn5a")

FiveDay.Starting.remove <- c("Ikbke", "Tns1", "Arl4a")

EmbryoInt4.adjusted <- EmbryoInt4[1:200]
EmbryoInt4.adjusted <- EmbryoInt4.adjusted[! EmbryoInt4.adjusted %in% Embryogenesis.remove]

StartInt4.adjusted <- StartInt4[1:200]
StartInt4.adjusted <- StartInt4.adjusted[! StartInt4.adjusted %in% Starting.remove]

P0Int4.adjusted <- P0Int4[! P0Int4 %in% P0.remove]

FiveDayInt4.adjusted <- FiveDayInt4[1:200]
FiveDayInt4.adjusted <- FiveDayInt4.adjusted[! FiveDayInt4.adjusted %in% FiveDay.remove]

FirstHeatmap.adjusted <- FirstHeatmap[! FirstHeatmap %in% Heatmap1.remove]
ThirdHeatmap.adjusted <- ThirdHeatmap[! ThirdHeatmap %in% Heatmap3.remove]
FourthHeatmap.adjusted <- FourthHeatmap[! FourthHeatmap %in% Heatmap4.remove]
Embryo24Hr.OverRest.adjusted <- Embryo24Hr.OverRest[! Embryo24Hr.OverRest %in% Embryo.24Hr.remove]
EmbryoP0.OverRest.adjusted <- EmbryoP0.OverRest[! EmbryoP0.OverRest %in% Embryo.P0.remove]
StartFiveDay.OverEmbryoP0.adjusted <- StartFiveDay.OverEmbryoP0[!StartFiveDay.OverEmbryoP0 %in% FiveDay.Starting.remove]
