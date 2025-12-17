## IMPROVED SEGMENTATION PARAMETERS and GLOBAL NORMALIZATION for ASSAY 1##
# SET WD
setwd('/gpfs/gibbs/project/raredon/vl325/CytospinQuant/ImageQuant.git/Vi Sandbox/Image Quantification/SimpleSeg')

# Load Packages
library(simpleSeg)
library(EBImage)
library(cytomapper)
library(BiocParallel)

BiocParallel::register(BiocParallel::MulticoreParam(40))

#### RUN IF IMAGES HAVE NOT BEEN LOADED AND NORMALIZED ####
# Get path to image directory
pathToImages <- "/gpfs/gibbs/project/raredon/vl325/CytospinQuant/ImageQuant.git/Vi Sandbox/Image Quantification/Cytospin images/EPCAM.Krt5.Sox9"
# Get directories of images
imageDirs <- dir(pathToImages, "P", full.names = TRUE)
names(imageDirs) <- dir(pathToImages, "P", full.names = FALSE)

# Get files in each directory
files <- files <- lapply(
  imageDirs,
  list.files,
  pattern = ".Raw",
  full.names = TRUE
)

# Read files with readImage from EBImage
images <- lapply(files, EBImage::readImage, as.is = TRUE)

# Convert to cytoImageList
images <- cytomapper::CytoImageList(images)
imagesnorm <- normalize(images, separateChannels = TRUE, separateImages = FALSE, ft = c(0,1), overwrite = TRUE)

mcols(images)$imageID <- names(images)
mcols(imagesnorm)$imageID <- names(images)


#### PROCESSING IMAGES ####
# Load previously normalized images
load("Assay 1 Norm Images File.Robj")
imagesnorm
### Segment ####
for (i in 5:9){
  masks[i] <- simpleSeg::simpleSeg(imagesnorm[i],
                                nucleus = "DAPI.Raw",
                                cellBody = c("EPCAM.Raw","Krt5.Raw","Ager.Raw"),
                                # sizeSelection = 250,
                                transform = "norm99",
                                watershed = "intensity",
                                tolerance = 0.0225,
                                smooth = 1.8,
                                # # ext = 3,
                                # # pca = TRUE,
                                cores = 40
                                )

save(masks, file = "NEW Assay 1 Masks.Robj")
}

for (i in 7:9){
  png(filename = paste(names(imagesnorm[i]), "NEW Image Segmentation with Outlines.png"), width = 8.5, height = 11 , units = "in", res = 450)
      cytomapper::plotPixels(
        image = imagesnorm[i],
        mask = masks[i],
        img_id = "imageID",
        colour_by = c("DAPI.Raw", "EPCAM.Raw", "Ager.Raw", "Krt5.Raw"),
        display = "single",
        colour = list(
          DAPI.Raw = c("black", "blue"),
          EPCAM.Raw = c("black", "green"),
          Ager.Raw = c("black", "red"),
          Krt5.Raw = c("black", "cyan")
                      ),
        image_title = list(
          cex = 1.25),
        legend = NULL
              )
      dev.off()
}      
    
    # Plot and check

    png(filename = paste(names(sing), "NEW Image Segmentation with Outlines.png", width = 8.5, height = 11 , units = "in", res = 450)
      cytomapper::plotPixels(
        image = sing,
        mask = masks,
        img_id = "imageID",
        colour_by = c("DAPI.Raw", "EPCAM.Raw", "Ager.Raw", "Krt5.Raw"),
        display = "single",
        colour = list(
          DAPI.Raw = c("black", "blue"),
          EPCAM.Raw = c("black", "green"),
          Ager.Raw = c("black", "red"),
          Krt5.Raw = c("black", "cyan")
        )
        ,
        # bcg = list(
        #   DAPI.Raw = c(0, 1.5, 1),
        #   #   #CD3 = c(0, 1, 1.5),
        #   EPCAM.Raw = c(0, 1, 1),
        #   Vimentin.Raw = c(0, 1, 1),
        #   Krt5.Raw = c(0, 1, 1)
        # )
        image_title = list(# text = c("image1", "image2", "image3"),
          # position = "topleft",
          # colour = "blue",
          # margin = c(0,5),
          # font = 2,
          cex = 1.25),
        scale_bar = NULL,
        legend = NULL
                            )
    dev.off()
    }
      
#### SUMMARIZE CELL FEATURES ####
cellSCE <- cytomapper::measureObjects(masks, imagesnorm, img_id = "imageID")
#cellInt <- cytomapper::measureObjects(masks, images, img_id = "imageID", feature_types = "basic")
save(cellSCE, file = "NEW Assay 1 Norm Cell Features.Robj")

# Extract marker data and bind with information about images
df.new.1 <- as.data.frame(cbind(colData(cellSCE), t(assay(cellSCE, "counts"))))
save(df.new.1, file = "NEW Assay 1 Norm Data Frame.Robj")


#### THRESHOLD ####
library(dplyr)
library(ggplot2)
library(EBImage)
library(ggthemes)
library(patchwork)

#Initialize Thresholds
dapi.thresh = 0
epcam.thresh = 0
ager.thresh = 0
krt5.thresh = 0
abca3.thresh = 0
sox2.thresh = 0
sox9.thresh = 0

#EXTRACT SPECIFIC MARKER DATA
#ASSAY 1
DAPI <- df %>%
  select(imageID:m.eccentricity, DAPI.Raw) %>%
  mutate(y.flip = 0 - m.cy,
         `Pos/Neg` = (if_else(DAPI.Raw >= dapi.thresh, true = "Pos", false = "Neg"))
  )

EpCAM <- df %>%
  select(imageID:m.eccentricity, EPCAM.Raw) %>%
  mutate(y.flip = 0 - m.cy,
         `Pos/Neg` = (if_else(EPCAM.Raw >= epcam.thresh, true = "Pos", false = "Neg"))
  )

Ager <- df %>%
  select(imageID:m.eccentricity, Ager.Raw) %>%
  mutate(y.flip = 0 - m.cy,
         `Pos/Neg` = (if_else(Ager.Raw >= ager.thresh, true = "Pos", false = "Neg"))
  )

Krt5 <- df %>%
  select(imageID:m.eccentricity, Krt5.Raw) %>%
  mutate(y.flip = 0 - m.cy,
         `Pos/Neg` = (if_else(Krt5.Raw >= krt5.thresh, true = "Pos", false = "Neg"))
  )

#ASSAY 2
ABCA3 <- df %>%
  filter(Feature == 'ABCA3',
         `Norm Intensity` >= abca3.thresh)

Sox2 <- df %>%
  filter(Feature == 'Sox2',
         `Norm Intensity` >= sox2.thresh)

Sox9 <- df %>%
  filter(Feature == 'Sox9',
         `Norm Intensity` >= sox9.thresh)

### PLOT POINTS ####
#Assay 1
DAPI.pt <- ggplot(data = DAPI, aes(x = m.cx, y = y.flip)) +
  geom_point(size = 2, aes(colour = `DAPI.Raw`)) +
  scale_colour_gradient(low = "lightgray", high = "darkblue") +
  # facet_grid(~Passage ~ ~Replicate) +
  ggtitle("DAPI") +
  coord_fixed() +
  theme_bw() +
  theme(legend.position = "bottom",
        # legend.key.size = unit(1.25, 'cm'), #change legend key size
        legend.key.height = unit(0.5, 'cm'), #change legend key height
        legend.key.width = unit(0.85, 'cm'), #change legend key width
        legend.title = element_blank(),
        legend.text = element_text(size=8)) #change legend text font size)
DAPI.pt

EpCAM.pt <- ggplot(data = EpCAM, aes(x = m.cx, y = y.flip)) +
  geom_point(size = 2, aes(colour = `EPCAM.Raw`)) +
  scale_colour_gradient(low = "lightgray", high = "green") +
  # facet_grid(~Passage ~ ~Replicate) +
  ggtitle("EpCAM") +
  coord_fixed() +
  theme_bw() +
  theme(legend.position = "bottom",
        # legend.key.size = unit(1.25, 'cm'), #change legend key size
        legend.key.height = unit(0.5, 'cm'), #change legend key height
        legend.key.width = unit(0.85, 'cm'), #change legend key width
        legend.title = element_blank(),
        legend.text = element_text(size=8)) #change legend text font size)
EpCAM.pt

Ager.pt <- ggplot(data = Ager, aes(x = m.cx, y = y.flip)) +
  geom_point(size = 2, aes(colour = `Ager.Raw`)) +
  scale_colour_gradient(low = "lightgray", high = "red") +
  # facet_grid(~Passage ~ ~Replicate) +
  ggtitle("Ager") +
  coord_fixed() +
  theme_bw() +
  theme(legend.position = "bottom",
        # legend.key.size = unit(1.25, 'cm'), #change legend key size
        legend.key.height = unit(0.5, 'cm'), #change legend key height
        legend.key.width = unit(0.85, 'cm'), #change legend key width
        legend.title = element_blank(),
        legend.text = element_text(size=8)) #change legend text font size)
Ager.pt

Krt5.pt <- ggplot(data = Krt5, aes(x = m.cx, y = y.flip)) +
  geom_point(size = 2, aes(colour = `Krt5.Raw`)) +
  scale_colour_gradient(low = "lightgray", high = "cyan") +
  # facet_grid(~Passage ~ ~Replicate) +
  ggtitle("Krt5") +
  coord_fixed() +
  theme_bw() +
  theme(legend.position = "bottom",
        # legend.key.size = unit(1.25, 'cm'), #change legend key size
        legend.key.height = unit(0.5, 'cm'), #change legend key height
        legend.key.width = unit(0.85, 'cm'), #change legend key width
        legend.title = element_blank(),
        legend.text = element_text(size=8)) #change legend text font size)
Krt5.pt

png(filename = "P1 B12 A1 Small Image Color Point Plots.png", width = 15, height = 5, units = "in", res = 450)
DAPI.pt|EpCAM.pt|Ager.pt|Krt5.pt
dev.off()


### Plot threshold points ####

#Set Thresholds
dapi.thresh = 0.1
epcam.thresh = 0.05
ager.thresh = 0.06
krt5.thresh = 0.0275
abca3.thresh = 0
sox2.thresh = 0
sox9.thresh = 0

# Update dataframes
DAPI <- DAPI %>%
  mutate(`Pos/Neg` = (if_else(DAPI.Raw >= dapi.thresh, true = "Pos", false = "Neg")))

EpCAM <- EpCAM %>%
  mutate(`Pos/Neg` = (if_else(EPCAM.Raw >= epcam.thresh, true = "Pos", false = "Neg")))

Ager <- Ager %>%
  mutate(`Pos/Neg` = (if_else(Ager.Raw >= ager.thresh, true = "Pos", false = "Neg")))

Krt5 <- Krt5 %>%
  mutate(`Pos/Neg` = (if_else(Krt5.Raw >= krt5.thresh, true = "Pos", false = "Neg")))

# Plot thresholded points
DAPI.t <- ggplot(DAPI, aes(x = m.cx, y = y.flip, colour = `Pos/Neg`)) +
  geom_point() +
  scale_colour_manual(values = c("gray", "darkblue")) +
  ggtitle("DAPI") +
  coord_fixed() +
  theme_bw() +
  theme(legend.position = "bottom",
        # legend.key.size = unit(0, 'cm'), #change legend key size
        legend.key.height = unit(0, 'cm'), #change legend key height
        legend.key.width = unit(0.5, 'cm'), #change legend key width
        legend.title = element_blank(),
        legend.text = element_text(size = 10)) + #change legend text font size)
  annotate("text", x = 750, y = 50, fontsize = 6,
           label = paste("Thresh = ", dapi.thresh), color = "red")


EpCAM.t <- ggplot(EpCAM, aes(x = m.cx, y = y.flip, colour = `Pos/Neg`)) +
  geom_point() +
  scale_colour_manual(values = c("gray", "green")) +
  ggtitle("EpCAM") +
  coord_fixed() +
  theme_bw()  +
  theme(legend.position = "bottom",
        # legend.key.size = unit(0, 'cm'), #change legend key size
        legend.key.height = unit(0, 'cm'), #change legend key height
        legend.key.width = unit(0.5, 'cm'), #change legend key width
        legend.title = element_blank(),
        legend.text = element_text(size = 10))+ #change legend text font size)
  annotate("text", x = 750, y = 50, fontsize = 6,
           label = paste("Thresh = ", epcam.thresh), color = "red")


Ager.t <- ggplot(Ager, aes(x = m.cx, y = y.flip, colour = `Pos/Neg`)) +
  geom_point() +
  scale_colour_manual(values = c("gray", "red")) +
  ggtitle("Ager") +
  coord_fixed() +
  theme_bw() +
  theme(legend.position = "bottom",
        # legend.key.size = unit(0, 'cm'), #change legend key size
        legend.key.height = unit(0, 'cm'), #change legend key height
        legend.key.width = unit(0.5, 'cm'), #change legend key width
        legend.title = element_blank(),
        legend.text = element_text(size = 10))+ #change legend text font size)
  annotate("text", x = 750, y = 50, fontsize = 6,
           label = paste("Thresh = ", ager.thresh), color = "red")


Krt5.t <- ggplot(Krt5, aes(x = m.cx, y = y.flip, colour = `Pos/Neg`)) +
  geom_point() +
  scale_colour_manual(values = c("gray", "cyan")) +
  ggtitle("Krt5") +
  coord_fixed() +
  theme_bw() +
  theme(legend.position = "bottom",
        # legend.key.size = unit(0, 'cm'), #change legend key size
        legend.key.height = unit(0, 'cm'), #change legend key height
        legend.key.width = unit(0.5, 'cm'), #change legend key width
        legend.title = element_blank(),
        legend.text = element_text(size = 10)) + #change legend text font size)
  annotate("text", x = 750, y = 50, fontsize = 6,
           label = paste("Thresh = ", krt5.thresh), color = "red")

# DAPI.t|EpCAM.t|Ager.t|Krt5.t
png(filename = "P1 B12 A1 Stacked Color w Higher Threshold Plots.png", width = 15, height = 8, units = "in", res = 450)
(DAPI.pt/DAPI.t)|(EpCAM.pt/EpCAM.t)|(Ager.pt/Ager.t)|(Krt5.pt/Krt5.t)
dev.off()