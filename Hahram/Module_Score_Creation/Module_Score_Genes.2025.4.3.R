####We are going to attempt to generate a series of module scores for our engineered data based off of existing native cell cluster markers.
####We begin by loading in the combined rat lung atlas. 

###Our approach will generate four markerlists for each module using FindMarkers; one comparing the target cluster against the entire object, another comparing the target cluster against the epithelium, another comparing the target cluster against the most spatially proximal epithelial cells, and finally one comparing the target cluster against the cluster most similar to itself. We then use intersect to find what common genes exist across all four markerlists, and use that final product as our module score feature list.

####Starting with ATI 
####Misc Testing
DimPlot(lung.combined, group.by = "CellClass_Final")

DimPlot(lung.combined, split.by = "CellClass_Final", label = T)
FeaturePlot(lung.combined, features = c("Scgb3a2", "Bpifa5", "Bpifb1", "Scgb1a1"), label = T)
FeaturePlot(lung.combined, features = c("Avil"), label = T)
###ATI
################
###ATI Over All
ATIOverAll <- FindMarkers(lung.combined, 
                          ident.1 = "ATI", 
                          ident.2 = NULL, 
                          only.pos = T, 
                          logfc.threshold = 0.1)
ATIOverAll$ratio <- ATIOverAll$pct.1/ATIOverAll$pct.2
ATIOverAll$power <-ATIOverAll$ratio*ATIOverAll$avg_log2FC

ATIOverAll <- ATIOverAll %>% 
  filter (p_val_adj < 0.05) %>% 
  filter (pct.1 >= 0.1) 

ATIOverAll <- ATIOverAll %>%  top_n(200, power)
ATIOverAll <- rownames(ATIOverAll)


ATIOverAll.toRemove <- grep(paste(c("AABR0", "LOC"), collapse = "|"), 
                            ATIOverAll, value = TRUE)
ATIOverAll <- ATIOverAll[! ATIOverAll %in% ATIOverAll.toRemove]

###ATII Over Epithelium
ATIOverEpithelium<- FindMarkers(lung.combined, 
                                only.pos = T, 
                                ident.1 = "ATI", 
                                ident.2 = c("Tuft", "BASC", "ATII-ATI", "Secretory", "Ciliated", "ATII"), 
                                logfc.threshold = 0.1)
ATIOverEpithelium$ratio <- ATIOverEpithelium$pct.1/ATIOverEpithelium$pct.2
ATIOverEpithelium$power <-ATIOverEpithelium$ratio*ATIOverEpithelium$avg_log2FC

ATIOverEpithelium <- ATIOverEpithelium %>% 
  filter (p_val_adj < 0.05) %>% 
  filter (pct.1 >= 0.1) 
ATIOverEpithelium <- ATIOverEpithelium %>%  top_n(200, power)
ATIOverEpithelium <- rownames(ATIOverEpithelium)


ATIOverEpithelium.toRemove <- grep(paste(c("AABR0", "LOC"), collapse = "|"), 
                                   ATIOverEpithelium, value = TRUE)
ATIOverEpithelium <- ATIOverEpithelium[! ATIOverEpithelium %in% ATIOverEpithelium.toRemove]

###ATI Over Proximal
ATIOverProximal<- FindMarkers(lung.combined, 
                              only.pos = T,
                              ident.1 = "ATI", 
                              ident.2 = c("Secretory", "Ciliated", "Tuft", "BASC"), 
                              logfc.threshold = 0.1)

ATIOverProximal$ratio <- ATIOverProximal$pct.1/ATIOverProximal$pct.2

ATIOverProximal$power <-ATIOverProximal$ratio*ATIOverProximal$avg_log2FC

ATIOverProximal <- ATIOverProximal %>% 
  filter (p_val_adj < 0.05) %>% 
  filter (pct.1 >= 0.1) 
ATIOverProximal <- ATIOverProximal %>%  top_n(200, power)
ATIOverProximal <- rownames(ATIOverProximal)


ATIOverProximal.toRemove <- grep(paste(c("AABR0", "LOC"), collapse = "|"), 
                                 ATIOverProximal, value = TRUE)
ATIOverProximal <- ATIOverProximal[! ATIOverProximal %in% ATIOverProximal.toRemove]

###ATI Over ATII

ATIOverATII<- FindMarkers(lung.combined, 
                          only.pos = T,
                          ident.1 = "ATI", 
                          ident.2 = "ATII", 
                          logfc.threshold = 0.1)
ATIOverATII$ratio <- ATIOverATII$pct.1/ATIOverATII$pct.2

ATIOverATII$power <-ATIOverATII$ratio*ATIOverATII$avg_log2FC

ATIOverATII <- ATIOverATII %>% 
  filter (p_val_adj < 0.05) %>% 
  filter (pct.1 >= 0.1) 

ATIOverATII <- ATIOverATII %>% 
  filter (p_val_adj < 0.05) %>% 
  filter (pct.1 >= 0.1) 
ATIOverATII <- ATIOverATII %>%  top_n(200, power)
ATIOverATII <- rownames(ATIOverATII)


ATIOverATII.toRemove <- grep(paste(c("AABR0", "LOC"), collapse = "|"), 
                             ATIOverATII, value = TRUE)
ATIOverATII <- ATIOverATII[! ATIOverATII %in% ATIOverATII.toRemove]

ATI_genes <- intersect(intersect(intersect(ATIOverAll, 
                                           ATIOverEpithelium),
                                 ATIOverProximal), 
                       ATIOverATII)

###For printing genes vertically
cat(ATIOverAll, sep = "\n")
cat(ATIOverEpithelium, sep = "\n")
cat(ATIOverProximal, sep = "\n")
cat(ATIOverATII, sep = "\n")
cat(ATI_genes, sep = "\n")
#######



###ATII
################
###ATII Over All
ATIIOverAll<- FindMarkers(lung.combined, 
                          ident.1 = "ATII", 
                          ident.2 = NULL, 
                          only.pos = T, 
                          logfc.threshold = 0.1)
ATIIOverAll$ratio <- ATIIOverAll$pct.1/ATIIOverAll$pct.2
ATIIOverAll$power <-ATIIOverAll$ratio*ATIIOverAll$avg_log2FC

ATIIOverAll <- ATIIOverAll %>% 
  filter (p_val_adj < 0.05) %>% 
  filter (pct.1 >= 0.1) 

ATIIOverAll <- ATIIOverAll %>%  top_n(200, power)
ATIIOverAll <- rownames(ATIIOverAll)


ATIIOverAll.toRemove <- grep(paste(c("AABR0", "LOC"), collapse = "|"), 
                             ATIIOverAll, value = TRUE)
ATIIOverAll <- ATIIOverAll[! ATIIOverAll %in% ATIIOverAll.toRemove]

###ATII Over Epithelium
ATIIOverEpithelium<- FindMarkers(lung.combined, 
                                 only.pos = T, 
                                 ident.1 = "ATII", 
                                 ident.2 = c("Tuft", "BASC", "ATII-ATI", "Secretory", "Ciliated", "ATI"), 
                                 logfc.threshold = 0.1)
ATIIOverEpithelium$ratio <- ATIIOverEpithelium$pct.1/ATIIOverEpithelium$pct.2
ATIIOverEpithelium$power <-ATIIOverEpithelium$ratio*ATIIOverEpithelium$avg_log2FC

ATIIOverEpithelium <- ATIIOverEpithelium %>% 
  filter (p_val_adj < 0.05) %>% 
  filter (pct.1 >= 0.1) 
ATIIOverEpithelium <- ATIIOverEpithelium %>%  top_n(200, power)
ATIIOverEpithelium <- rownames(ATIIOverEpithelium)


ATIIOverEpithelium.toRemove <- grep(paste(c("AABR0", "LOC"), collapse = "|"), 
                                    ATIIOverEpithelium, value = TRUE)
ATIIOverEpithelium <- ATIIOverEpithelium[! ATIIOverEpithelium %in% ATIIOverEpithelium.toRemove]

###ATII Over Proximal
ATIIOverProximal<- FindMarkers(lung.combined, 
                               only.pos = T,
                               ident.1 = "ATII", 
                               ident.2 = c("Secretory", "Ciliated", "Tuft", "BASC"), 
                               logfc.threshold = 0.1)

ATIIOverProximal$ratio <- ATIIOverProximal$pct.1/ATIIOverProximal$pct.2

ATIIOverProximal$power <-ATIIOverProximal$ratio*ATIIOverProximal$avg_log2FC

ATIIOverProximal <- ATIIOverProximal %>% 
  filter (p_val_adj < 0.05) %>% 
  filter (pct.1 >= 0.1) 
ATIIOverProximal <- ATIIOverProximal %>%  top_n(200, power)
ATIIOverProximal <- rownames(ATIIOverProximal)


ATIIOverProximal.toRemove <- grep(paste(c("AABR0", "LOC"), collapse = "|"), 
                                  ATIIOverProximal, value = TRUE)
ATIIOverProximal <- ATIIOverProximal[! ATIIOverProximal %in% ATIIOverProximal.toRemove]

###ATII Over ATI
ATIIOverATI<- FindMarkers(lung.combined, 
                          only.pos = T,
                          ident.1 = "ATII", 
                          ident.2 = "ATI", 
                          logfc.threshold = 0.1)
ATIIOverATI$ratio <- ATIIOverATI$pct.1/ATIIOverATI$pct.2

ATIIOverATI$power <-ATIIOverATI$ratio*ATIIOverATI$avg_log2FC

ATIIOverATI <- ATIIOverATI %>% 
  filter (p_val_adj < 0.05) %>% 
  filter (pct.1 >= 0.1) 

ATIIOverATI <- ATIIOverATI %>% 
  filter (p_val_adj < 0.05) %>% 
  filter (pct.1 >= 0.1) 
ATIIOverATI <- ATIIOverATI %>%  top_n(200, power)
ATIIOverATI <- rownames(ATIIOverATI)


ATIIOverATI.toRemove <- grep(paste(c("AABR0", "LOC"), collapse = "|"), 
                             ATIIOverATI, value = TRUE)
ATIIOverATI <- ATIIOverATI[! ATIIOverATI %in% ATIIOverATI.toRemove]

ATII_genes <- intersect(intersect(intersect(ATIIOverAll, 
                                            ATIIOverEpithelium),
                                  ATIIOverProximal), 
                        ATIIOverATI)

cat(ATIIOverAll, sep = "\n")
cat(ATIIOverEpithelium, sep = "\n")
cat(ATIIOverProximal, sep = "\n")
cat(ATIIOverATI, sep = "\n")
cat(ATII_genes, sep = "\n")
#######

#Ciliated
##############
CiliatedOverAll<- FindMarkers(lung.combined, 
                              ident.1 = "Ciliated", 
                              ident.2 = NULL, 
                              only.pos = T, 
                              logfc.threshold = 0.1)
CiliatedOverAll$ratio <- CiliatedOverAll$pct.1/CiliatedOverAll$pct.2
CiliatedOverAll$power <-CiliatedOverAll$ratio*CiliatedOverAll$avg_log2FC

CiliatedOverAll <- CiliatedOverAll %>% 
  filter (p_val_adj < 0.05) %>% 
  filter (pct.1 >= 0.1) 

CiliatedOverAll <- CiliatedOverAll %>%  top_n(200, power)
CiliatedOverAll <- rownames(CiliatedOverAll)


CiliatedOverAll.toRemove <- grep(paste(c("AABR0", "LOC"), collapse = "|"), 
                                 CiliatedOverAll, value = TRUE)
CiliatedOverAll <- CiliatedOverAll[! CiliatedOverAll %in% CiliatedOverAll.toRemove]

###Ciliated Over Epithelium
CiliatedOverEpithelium<- FindMarkers(lung.combined, 
                                     only.pos = T, 
                                     ident.1 = "Ciliated", 
                                     ident.2 = c("Tuft", "BASC", "ATII-ATI", "Secretory", "ATII", "ATI"), 
                                     logfc.threshold = 0.1)
CiliatedOverEpithelium$ratio <- CiliatedOverEpithelium$pct.1/CiliatedOverEpithelium$pct.2

CiliatedOverEpithelium$power <-CiliatedOverEpithelium$ratio*CiliatedOverEpithelium$avg_log2FC

CiliatedOverEpithelium <- CiliatedOverEpithelium %>% 
  filter (p_val_adj < 0.05) %>% 
  filter (pct.1 >= 0.1) 

CiliatedOverEpithelium <- CiliatedOverEpithelium %>%  top_n(200, power)
CiliatedOverEpithelium <- rownames(CiliatedOverEpithelium)


CiliatedOverEpithelium.toRemove <- grep(paste(c("AABR0", "LOC"), collapse = "|"), 
                                 CiliatedOverEpithelium, value = TRUE)
CiliatedOverEpithelium <- CiliatedOverEpithelium[! CiliatedOverEpithelium %in% CiliatedOverEpithelium.toRemove]

###Ciliated Over ATI and ATII
CiliatedOverProximal<- FindMarkers(lung.combined, 
                                   only.pos = T,
                                   ident.1 = "Ciliated", 
                                   ident.2 = c("ATI", "ATII-ATI", "ATII"), 
                                   logfc.threshold = 0.1)
CiliatedOverProximal$ratio <- CiliatedOverProximal$pct.1/CiliatedOverProximal$pct.2

CiliatedOverProximal$power <-CiliatedOverProximal$ratio*CiliatedOverProximal$avg_log2FC

CiliatedOverProximal <- CiliatedOverProximal %>% 
  filter (p_val_adj < 0.05) %>% 
  filter (pct.1 >= 0.1) 


CiliatedOverProximal <- CiliatedOverProximal %>%  top_n(200, power)
CiliatedOverProximal <- rownames(CiliatedOverProximal)


CiliatedOverProximal.toRemove <- grep(paste(c("AABR0", "LOC"), collapse = "|"), 
                                      CiliatedOverProximal, value = TRUE)
CiliatedOverProximal <- CiliatedOverProximal[! CiliatedOverProximal %in% CiliatedOverProximal.toRemove]


###Ciliated Over Secretory
CiliatedOverSecretory<- FindMarkers(lung.combined, 
                                    only.pos = T,
                                    ident.1 = "Ciliated", 
                                    ident.2 = "Secretory", 
                                    logfc.threshold = 0.1)

CiliatedOverSecretory$ratio <- CiliatedOverSecretory$pct.1/CiliatedOverSecretory$pct.2

CiliatedOverSecretory$power <-CiliatedOverSecretory$ratio*CiliatedOverSecretory$avg_log2FC

CiliatedOverSecretory <- CiliatedOverSecretory %>% 
  filter (p_val_adj < 0.05) %>% 
  filter (pct.1 >= 0.1) 

CiliatedOverSecretory <- CiliatedOverSecretory %>%  top_n(200, power)
CiliatedOverSecretory <- rownames(CiliatedOverSecretory)


CiliatedOverSecretory.toRemove <- grep(paste(c("AABR0", "LOC"), collapse = "|"), 
                                       CiliatedOverSecretory, value = TRUE)
CiliatedOverSecretory <- CiliatedOverSecretory[! CiliatedOverSecretory %in% CiliatedOverSecretory.toRemove]


Ciliated_genes <- intersect(intersect(intersect(CiliatedOverAll, CiliatedOverEpithelium), CiliatedOverProximal), CiliatedOverSecretory)

cat(CiliatedOverAll, sep = "\n")
cat(CiliatedOverEpithelium, sep = "\n")
cat(CiliatedOverProximal, sep = "\n")
cat(CiliatedOverSecretory, sep = "\n")
cat(Ciliated_genes, sep = "\n")

#######

#####Secretory
###################
SecretoryOverAll<- FindMarkers(lung.combined, 
                              ident.1 = "Secretory", 
                              ident.2 = NULL, 
                              only.pos = T, 
                              logfc.threshold = 0.1)
SecretoryOverAll$ratio <- SecretoryOverAll$pct.1/SecretoryOverAll$pct.2
SecretoryOverAll$power <-SecretoryOverAll$ratio*SecretoryOverAll$avg_log2FC

SecretoryOverAll <- SecretoryOverAll %>% 
  filter (p_val_adj < 0.05) %>% 
  filter (pct.1 >= 0.1) 

SecretoryOverAll <- SecretoryOverAll %>%  top_n(200, power)
SecretoryOverAll <- rownames(SecretoryOverAll)


SecretoryOverAll.toRemove <- grep(paste(c("AABR0", "LOC"), collapse = "|"), 
                                 SecretoryOverAll, value = TRUE)
SecretoryOverAll <- SecretoryOverAll[! SecretoryOverAll %in% SecretoryOverAll.toRemove]

###Secretory Over Epithelium
SecretoryOverEpithelium<- FindMarkers(lung.combined, 
                                     only.pos = T, 
                                     ident.1 = "Secretory", 
                                     ident.2 = c("Tuft", "BASC", "ATII-ATI", "Ciliated", "ATII", "ATI"), 
                                     logfc.threshold = 0.1)
SecretoryOverEpithelium$ratio <- SecretoryOverEpithelium$pct.1/SecretoryOverEpithelium$pct.2

SecretoryOverEpithelium$power <-SecretoryOverEpithelium$ratio*SecretoryOverEpithelium$avg_log2FC

SecretoryOverEpithelium <- SecretoryOverEpithelium %>% 
  filter (p_val_adj < 0.05) %>% 
  filter (pct.1 >= 0.1) 

SecretoryOverEpithelium <- SecretoryOverEpithelium %>%  top_n(200, power)
SecretoryOverEpithelium <- rownames(SecretoryOverEpithelium)


SecretoryOverEpithelium.toRemove <- grep(paste(c("AABR0", "LOC"), collapse = "|"), 
                                        SecretoryOverEpithelium, value = TRUE)
SecretoryOverEpithelium <- SecretoryOverEpithelium[! SecretoryOverEpithelium %in% SecretoryOverEpithelium.toRemove]

###Secretory Over ATI and ATII
SecretoryOverProximal<- FindMarkers(lung.combined, 
                                   only.pos = T,
                                   ident.1 = "Secretory", 
                                   ident.2 = c("ATI", "ATII-ATI", "ATII"), 
                                   logfc.threshold = 0.1)
SecretoryOverProximal$ratio <- SecretoryOverProximal$pct.1/SecretoryOverProximal$pct.2

SecretoryOverProximal$power <-SecretoryOverProximal$ratio*SecretoryOverProximal$avg_log2FC

SecretoryOverProximal <- SecretoryOverProximal %>% 
  filter (p_val_adj < 0.05) %>% 
  filter (pct.1 >= 0.1) 


SecretoryOverProximal <- SecretoryOverProximal %>%  top_n(200, power)
SecretoryOverProximal <- rownames(SecretoryOverProximal)


SecretoryOverProximal.toRemove <- grep(paste(c("AABR0", "LOC"), collapse = "|"), 
                                      SecretoryOverProximal, value = TRUE)
SecretoryOverProximal <- SecretoryOverProximal[! SecretoryOverProximal %in% SecretoryOverProximal.toRemove]


###Secretory Over Ciliated
SecretoryOverCiliated<- FindMarkers(lung.combined, 
                                    only.pos = T,
                                    ident.1 = "Secretory", 
                                    ident.2 = "Ciliated", 
                                    logfc.threshold = 0.1)

SecretoryOverCiliated$ratio <- SecretoryOverCiliated$pct.1/SecretoryOverCiliated$pct.2

SecretoryOverCiliated$power <-SecretoryOverCiliated$ratio*SecretoryOverCiliated$avg_log2FC

SecretoryOverCiliated <- SecretoryOverCiliated %>% 
  filter (p_val_adj < 0.05) %>% 
  filter (pct.1 >= 0.1) 

SecretoryOverCiliated <- SecretoryOverCiliated %>%  top_n(200, power)
SecretoryOverCiliated <- rownames(SecretoryOverCiliated)


SecretoryOverCiliated.toRemove <- grep(paste(c("AABR0", "LOC"), collapse = "|"), 
                                       SecretoryOverCiliated, value = TRUE)
SecretoryOverCiliated <- SecretoryOverCiliated[! SecretoryOverCiliated %in% SecretoryOverCiliated.toRemove]


Secretory_genes <- intersect(intersect(intersect(SecretoryOverAll, SecretoryOverEpithelium), SecretoryOverProximal), SecretoryOverCiliated)

cat(SecretoryOverAll, sep = "\n")
cat(SecretoryOverEpithelium, sep = "\n")
cat(SecretoryOverProximal, sep = "\n")
cat(SecretoryOverCiliated, sep = "\n")
cat(Secretory_genes, sep = "\n")

########


#ATII-ATI
#####################
ATII_ATIOverAll<- FindMarkers(lung.combined, 
                              ident.1 = "ATII-ATI", 
                              ident.2 = NULL, 
                              only.pos = T, 
                              logfc.threshold = 0.1)
ATII_ATIOverAll$ratio <- ATII_ATIOverAll$pct.1/ATII_ATIOverAll$pct.2
ATII_ATIOverAll$power <-ATII_ATIOverAll$ratio*ATII_ATIOverAll$avg_log2FC

ATII_ATIOverAll <- ATII_ATIOverAll %>% 
  filter (p_val_adj < 0.05) %>% 
  filter (pct.1 >= 0.1) 

ATII_ATIOverAll <- ATII_ATIOverAll %>%  top_n(200, power)
ATII_ATIOverAll <- rownames(ATII_ATIOverAll)


ATII_ATIOverAll.toRemove <- grep(paste(c("AABR0", "LOC"), collapse = "|"), 
                                 ATII_ATIOverAll, value = TRUE)
ATII_ATIOverAll <- ATII_ATIOverAll[! ATII_ATIOverAll %in% ATII_ATIOverAll.toRemove]

###ATII-ATI Over Epithelium
ATII_ATIOverEpithelium<- FindMarkers(lung.combined, 
                                     only.pos = T, 
                                     ident.1 = "ATII-ATI", 
                                     ident.2 = c("Tuft", "BASC", "Secretory", "Ciliated", "ATII", "ATI"), 
                                     logfc.threshold = 0.1)
ATII_ATIOverEpithelium$ratio <- ATII_ATIOverEpithelium$pct.1/ATII_ATIOverEpithelium$pct.2

ATII_ATIOverEpithelium$power <-ATII_ATIOverEpithelium$ratio*ATII_ATIOverEpithelium$avg_log2FC

ATII_ATIOverEpithelium <- ATII_ATIOverEpithelium %>% 
  filter (p_val_adj < 0.05) %>% 
  filter (pct.1 >= 0.1) 

ATII_ATIOverEpithelium <- ATII_ATIOverEpithelium %>%  top_n(200, power)
ATII_ATIOverEpithelium <- rownames(ATII_ATIOverEpithelium)


ATII_ATIOverEpithelium.toRemove <- grep(paste(c("AABR0", "LOC"), collapse = "|"), 
                                        ATII_ATIOverEpithelium, value = TRUE)
ATII_ATIOverEpithelium <- ATII_ATIOverEpithelium[! ATII_ATIOverEpithelium %in% ATII_ATIOverEpithelium.toRemove]

###ATII-ATI Over Proximal
ATII_ATIOverProximal<- FindMarkers(lung.combined, 
                                   only.pos = T,
                                   ident.1 = "ATII-ATI", 
                                   ident.2 = c("Tuft", "BASC", "Secretory", "Ciliated"), 
                                   logfc.threshold = 0.1)
ATII_ATIOverProximal$ratio <- ATII_ATIOverProximal$pct.1/ATII_ATIOverProximal$pct.2

ATII_ATIOverProximal$power <-ATII_ATIOverProximal$ratio*ATII_ATIOverProximal$avg_log2FC

ATII_ATIOverProximal <- ATII_ATIOverProximal %>% 
  filter (p_val_adj < 0.05) %>% 
  filter (pct.1 >= 0.1) 


ATII_ATIOverProximal <- ATII_ATIOverProximal %>%  top_n(200, power)
ATII_ATIOverProximal <- rownames(ATII_ATIOverProximal)


ATII_ATIOverProximal.toRemove <- grep(paste(c("AABR0", "LOC"), collapse = "|"), 
                                      ATII_ATIOverProximal, value = TRUE)
ATII_ATIOverProximal <- ATII_ATIOverProximal[! ATII_ATIOverProximal %in% ATII_ATIOverProximal.toRemove]


###ATII-ATI Over ATI
ATII_ATIOverATI<- FindMarkers(lung.combined, 
                              only.pos = T,
                              ident.1 = "ATII-ATI", 
                              ident.2 = c("ATI"),
                              logfc.threshold = 0.1)

ATII_ATIOverATI$ratio <- ATII_ATIOverATI$pct.1/ATII_ATIOverATI$pct.2

ATII_ATIOverATI$power <-ATII_ATIOverATI$ratio*ATII_ATIOverATI$avg_log2FC

ATII_ATIOverATI <- ATII_ATIOverATI %>% 
  filter (p_val_adj < 0.05) %>% 
  filter (pct.1 >= 0.1) 

ATII_ATIOverATI <- ATII_ATIOverATI %>%  top_n(200, power)
ATII_ATIOverATI <- rownames(ATII_ATIOverATI)


ATII_ATIOverATI.toRemove <- grep(paste(c("AABR0", "LOC"), collapse = "|"), 
                                 ATII_ATIOverATI, value = TRUE)
ATII_ATIOverATI <- ATII_ATIOverATI[! ATII_ATIOverATI %in% ATII_ATIOverATI.toRemove]


###ATII-ATI Over ATII
ATII_ATIOverATII<- FindMarkers(lung.combined, 
                               only.pos = T,
                               ident.1 = "ATII-ATI", 
                               ident.2 = c("ATII"),
                               logfc.threshold = 0.1)

ATII_ATIOverATII$ratio <- ATII_ATIOverATII$pct.1/ATII_ATIOverATII$pct.2

ATII_ATIOverATII$power <-ATII_ATIOverATII$ratio*ATII_ATIOverATII$avg_log2FC

ATII_ATIOverATII <- ATII_ATIOverATII %>% 
  filter (p_val_adj < 0.05) %>% 
  filter (pct.1 >= 0.1) 

ATII_ATIOverATII <- ATII_ATIOverATII %>%  top_n(200, power)
ATII_ATIOverATII <- rownames(ATII_ATIOverATII)


ATII_ATIOverATII.toRemove <- grep(paste(c("AABR0", "LOC"), collapse = "|"), 
                                  ATII_ATIOverATII, value = TRUE)
ATII_ATIOverATII <- ATII_ATIOverATII[! ATII_ATIOverATII %in% ATII_ATIOverATII.toRemove]


ATII_ATI_genes <- intersect(intersect(intersect(
  ATII_ATIOverAll, ATII_ATIOverEpithelium), 
  ATII_ATIOverProximal),
  ATII_ATIOverATII)

cat(ATII_ATIOverAll, sep = "\n")
cat(ATII_ATIOverEpithelium, sep = "\n")
cat(ATII_ATIOverProximal, sep = "\n")
cat(ATII_ATIOverATI, sep = "\n")
cat(ATII_ATIOverATII, sep = "\n")
cat(ATII_ATI_genes, sep = "\n")
#############


#BASC
#########################
BASCOverAll<- FindMarkers(lung.combined, 
                          ident.1 = "BASC", 
                          ident.2 = NULL, 
                          only.pos = T, 
                          logfc.threshold = 0.1)
BASCOverAll$ratio <- BASCOverAll$pct.1/BASCOverAll$pct.2
BASCOverAll$power <-BASCOverAll$ratio*BASCOverAll$avg_log2FC

BASCOverAll <- BASCOverAll %>% 
  filter (p_val_adj < 0.05) %>% 
  filter (pct.1 >= 0.1) 

BASCOverAll <- BASCOverAll %>%  top_n(200, power)
BASCOverAll <- rownames(BASCOverAll)


BASCOverAll.toRemove <- grep(paste(c("AABR0", "LOC"), collapse = "|"), 
                             BASCOverAll, value = TRUE)
BASCOverAll <- BASCOverAll[! BASCOverAll %in% BASCOverAll.toRemove]

###BASC Over Epithelium
BASCOverEpithelium<- FindMarkers(lung.combined, 
                                 only.pos = T, 
                                 ident.1 = "BASC", 
                                 ident.2 = c("Tuft", "ATII-ATI", "Secretory", "Ciliated", "ATII", "ATI"), 
                                 logfc.threshold = 0.1)
BASCOverEpithelium$ratio <- BASCOverEpithelium$pct.1/BASCOverEpithelium$pct.2

BASCOverEpithelium$power <-BASCOverEpithelium$ratio*BASCOverEpithelium$avg_log2FC

BASCOverEpithelium <- BASCOverEpithelium %>% 
  filter (p_val_adj < 0.05) %>% 
  filter (pct.1 >= 0.1) 

BASCOverEpithelium <- BASCOverEpithelium %>%  top_n(200, power)
BASCOverEpithelium <- rownames(BASCOverEpithelium)


BASCOverEpithelium.toRemove <- grep(paste(c("AABR0", "LOC"), collapse = "|"), 
                                    BASCOverEpithelium, value = TRUE)
BASCOverEpithelium <- BASCOverEpithelium[! BASCOverEpithelium %in% BASCOverEpithelium.toRemove]

###BASC Over Proximal
BASCOverProximal<- FindMarkers(lung.combined, 
                               only.pos = T,
                               ident.1 = "BASC", 
                               ident.2 = c("ATII", "ATI", "ATII-ATI"), 
                               logfc.threshold = 0.1)
BASCOverProximal$ratio <- BASCOverProximal$pct.1/BASCOverProximal$pct.2

BASCOverProximal$power <-BASCOverProximal$ratio*BASCOverProximal$avg_log2FC

BASCOverProximal <- BASCOverProximal %>% 
  filter (p_val_adj < 0.05) %>% 
  filter (pct.1 >= 0.1) 


BASCOverProximal <- BASCOverProximal %>%  top_n(200, power)
BASCOverProximal <- rownames(BASCOverProximal)


BASCOverProximal.toRemove <- grep(paste(c("AABR0", "LOC"), collapse = "|"), 
                                  BASCOverProximal, value = TRUE)
BASCOverProximal <- BASCOverProximal[! BASCOverProximal %in% BASCOverProximal.toRemove]


###BASC Over "Tuft", "Secretory", "Ciliated"
BASCOverTuftSecretoryCiliated<- FindMarkers(lung.combined, 
                                            only.pos = T,
                                            ident.1 = "BASC", 
                                            ident.2 = c("Tuft", "Secretory", "Ciliated"),
                                            logfc.threshold = 0.1)

BASCOverTuftSecretoryCiliated$ratio <- BASCOverTuftSecretoryCiliated$pct.1/BASCOverTuftSecretoryCiliated$pct.2

BASCOverTuftSecretoryCiliated$power <-BASCOverTuftSecretoryCiliated$ratio*BASCOverTuftSecretoryCiliated$avg_log2FC

BASCOverTuftSecretoryCiliated <- BASCOverTuftSecretoryCiliated %>% 
  filter (p_val_adj < 0.05) %>% 
  filter (pct.1 >= 0.1) 

BASCOverTuftSecretoryCiliated <- BASCOverTuftSecretoryCiliated %>%  top_n(200, power)
BASCOverTuftSecretoryCiliated <- rownames(BASCOverTuftSecretoryCiliated)


BASCOverTuftSecretoryCiliated.toRemove <- grep(paste(c("AABR0", "LOC"), collapse = "|"), 
                                               BASCOverTuftSecretoryCiliated, value = TRUE)
BASCOverTuftSecretoryCiliated <- BASCOverTuftSecretoryCiliated[! BASCOverTuftSecretoryCiliated %in% BASCOverTuftSecretoryCiliated.toRemove]

BASC_genes <- intersect(intersect(intersect(
  BASCOverAll, BASCOverEpithelium), 
  BASCOverProximal),
  BASCOverTuftSecretoryCiliated)
BASC_genes <- intersect(
  BASCOverAll, BASCOverEpithelium) 

cat(BASCOverAll, sep = "\n")
cat(BASCOverEpithelium, sep = "\n")
cat(BASCOverProximal, sep = "\n")
cat(BASCOverTuftSecretoryCiliated, sep = "\n")
cat(BASC_genes, sep = "\n")

#############


##Tuft
#######################
TuftOverAll<- FindMarkers(lung.combined, 
                          ident.1 = "Tuft", 
                          ident.2 = NULL, 
                          only.pos = T, 
                          logfc.threshold = 0.1)
TuftOverAll$ratio <- TuftOverAll$pct.1/TuftOverAll$pct.2
TuftOverAll$power <-TuftOverAll$ratio*TuftOverAll$avg_log2FC

TuftOverAll <- TuftOverAll %>% 
  filter (p_val_adj < 0.05) %>% 
  filter (pct.1 >= 0.1) 

TuftOverAll <- TuftOverAll %>%  top_n(200, power)
TuftOverAll <- rownames(TuftOverAll)


TuftOverAll.toRemove <- grep(paste(c("AABR0", "LOC"), collapse = "|"), 
                             TuftOverAll, value = TRUE)
TuftOverAll <- TuftOverAll[! TuftOverAll %in% TuftOverAll.toRemove]

###Tufft Over Epithelium
TuftOverEpithelium<- FindMarkers(lung.combined, 
                                 only.pos = T, 
                                 ident.1 = "Tuft", 
                                 ident.2 = c("BASC", "ATII-ATI", "Secretory", "Ciliated", "ATII", "ATI"), 
                                 logfc.threshold = 0.1)
TuftOverEpithelium$ratio <- TuftOverEpithelium$pct.1/TuftOverEpithelium$pct.2

TuftOverEpithelium$power <-TuftOverEpithelium$ratio*TuftOverEpithelium$avg_log2FC

TuftOverEpithelium <- TuftOverEpithelium %>% 
  filter (p_val_adj < 0.05) %>% 
  filter (pct.1 >= 0.1) 

TuftOverEpithelium <- TuftOverEpithelium %>%  top_n(200, power)
TuftOverEpithelium <- rownames(TuftOverEpithelium)


TuftOverEpithelium.toRemove <- grep(paste(c("AABR0", "LOC"), collapse = "|"), 
                                    TuftOverEpithelium, value = TRUE)
TuftOverEpithelium <- TuftOverEpithelium[! TuftOverEpithelium %in% TuftOverEpithelium.toRemove]

###Tufft Over Proximal
TuftOverProximal<- FindMarkers(lung.combined, 
                               only.pos = T,
                               ident.1 = "Tuft", 
                               ident.2 = c("ATII", "ATI", "ATII-ATI"), 
                               logfc.threshold = 0.1)
TuftOverProximal$ratio <- TuftOverProximal$pct.1/TuftOverProximal$pct.2

TuftOverProximal$power <-TuftOverProximal$ratio*TuftOverProximal$avg_log2FC

TuftOverProximal <- TuftOverProximal %>% 
  filter (p_val_adj < 0.05) %>% 
  filter (pct.1 >= 0.1) 


TuftOverProximal <- TuftOverProximal %>%  top_n(200, power)
TuftOverProximal <- rownames(TuftOverProximal)


TuftOverProximal.toRemove <- grep(paste(c("AABR0", "LOC"), collapse = "|"), 
                                  TuftOverProximal, value = TRUE)
TuftOverProximal <- TuftOverProximal[! TuftOverProximal %in% TuftOverProximal.toRemove]


###Tuft Over BASC
TuftOverBASC<- FindMarkers(lung.combined, 
                           only.pos = T,
                           ident.1 = "Tuft", 
                           ident.2 = "BASC",
                           logfc.threshold = 0.1)

TuftOverBASC$ratio <- TuftOverBASC$pct.1/TuftOverBASC$pct.2

TuftOverBASC$power <-TuftOverBASC$ratio*TuftOverBASC$avg_log2FC

TuftOverBASC <- TuftOverBASC %>% 
  filter (p_val_adj < 0.05) %>% 
  filter (pct.1 >= 0.1) 

TuftOverBASC <- TuftOverBASC %>%  top_n(200, power)
TuftOverBASC <- rownames(TuftOverBASC)


TuftOverBASC.toRemove <- grep(paste(c("AABR0", "LOC"), collapse = "|"), 
                              TuftOverBASC, value = TRUE)
TuftOverBASC <- TuftOverBASC[! TuftOverBASC %in% TuftOverBASC.toRemove]

Tuft_genes <- intersect(intersect(intersect(
  TuftOverAll, TuftOverEpithelium), 
  TuftOverProximal),
  TuftOverBASC)


cat(TuftOverAll, sep = "\n")
cat(TuftOverEpithelium, sep = "\n")
cat(TuftOverProximal, sep = "\n")
cat(TuftOverBASC, sep = "\n")
cat(Tuft_genes, sep = "\n")

##############