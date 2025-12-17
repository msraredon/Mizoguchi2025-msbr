## IMPROVED SEGMENTATION PARAMETERS and GLOBAL NORMALIZATION for ASSAY 2 ##

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
pathToImages <- "/gpfs/gibbs/project/raredon/vl325/CytospinQuant/ImageQuant.git/Vi Sandbox/Image Quantification/Cytospin images/ABCA3.SOX2.SOX9"
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
load("Assay 2 Norm Images File.Robj")
imagesnorm
### Segment ####
for (i in 5:9){
  masks[i] <- simpleSeg::simpleSeg(imagesnorm[i],
                                   nucleus = "DAPI.Raw",
                                   cellBody = c("ABCA3.Raw", "Sox2.Raw"),
                                   # sizeSelection = 250,
                                   transform = "norm99",
                                   watershed = "intensity",
                                   tolerance = 0.005,
                                   smooth = 5,
                                   # ext = 20,
                                   pca = TRUE,
                                   cores = 20
                                  )
  save(masks, file = "NEW Assay 2 Masks.Robj")
}

for (i in 1:9){
png(filename = paste(names(imagesnorm[i]), "NEW Image Segmentation with Outlines.png"), width = 8.5, height = 11 , units = "in", res = 450)
cytomapper::plotPixels(
  image = imagesnorm[i],
  mask = masks[i],
  img_id = "imageID",
  colour_by = c("DAPI.Raw", "ABCA3.Raw", "Sox2.Raw", "Sox9.Raw"),
  display = "single",
  colour = list(
    DAPI.Raw = c("black", "blue"),
    #CD3 = c("black", "purple"),
    ABCA3.Raw = c("black", "green"),
    Sox2.Raw = c("black", "red"),
    Sox9.Raw = c("black", "cyan")
  ),
  image_title = list(
    cex = 1.25),
  legend = NULL
)
dev.off()
}



png(filename = paste(names(imagesnorm[3]), "NEW Image Segmentation with Outlines.png"), width = 8.5, height = 11 , units = "in", res = 450)
cytomapper::plotPixels(
  image = imagesnorm[3],
  mask = masks[3],
  img_id = "imageID",
  colour_by = c("DAPI.Raw", "ABCA3.Raw", "Sox2.Raw", "Sox9.Raw"),
  display = "single",
  colour = list(
    DAPI.Raw = c("black", "blue"),
    #CD3 = c("black", "purple"),
    ABCA3.Raw = c("black", "green"),
    Sox2.Raw = c("black", "red"),
    Sox9.Raw = c("black", "cyan")
  ),
  image_title = list(
    cex = 1.25),
  legend = NULL
)
dev.off()




# Plot and check
{
  png(filename = paste(names(sing), "NEW Image Segmentation with Outlines.png", width = 8.5, height = 11 , units = "in", res = 450)
      cytomapper::plotPixels(
        image = sing,
        mask = masks,
        img_id = "imageID",
        colour_by = c("DAPI.Raw", "ABCA3.Raw", "Sox2.Raw", "Sox9.Raw"),
        display = "single",
        colour = list(
          DAPI.Raw = c("black", "blue"),
          #CD3 = c("black", "purple"),
          ABCA3.Raw = c("black", "green"),
          Sox2.Raw = c("black", "red"),
          Sox9.Raw = c("black", "cyan")
        )
        ,
        # bcg = list(
        #   DAPI.Raw = c(0, 1.5, 1),
        #   #   #CD3 = c(0, 1, 1.5),
        #   ABCA3.Raw = c(0, 1, 1),
        #   Vimentin.Raw = c(0, 1, 1),
        #   Sox9.Raw = c(0, 1, 1)
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
save(cellSCE, file = "NEW Assay 2 Norm Cell Features.Robj")

# Extract marker data and bind with information about images
df.new.2 <- as.data.frame(cbind(colData(cellSCE), t(assay(cellSCE, "counts"))))
save(df.new.2, file = "NEW Assay 2 Norm Data Frame.Robj")


