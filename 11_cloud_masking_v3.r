# CLOUD MASKING

rm(list=ls(all=TRUE))

require(rgdal)
require(gdalUtils)
#require(raster)

Server <- TRUE
WriteLog <- TRUE

# IMAGE INPUT
#################################################################
args <- commandArgs(trailingOnly = TRUE)
input <- Sys.getenv('input') # input="2014-05-28T104747_RE5_1B-NAC_17275562_208298"

# DEFINE CLOUD COVERAGE FROM METADATA
image <- paste(input, "_P001_MUL", ".tif", sep = "")
imagebase <- unlist(strsplit(image, "[.]"))[1]

if(Server==TRUE){ 
	root <- Sys.getenv('baseDG')
} else{
	root <- "S:/derived/DG_v8"
}
	
# Set input & output directories
pathin <- paste(root, "6_image_matching", input, sep="/")
pathout <- paste(root, "7_cloud_mask", input, sep="/")
pathmeta <- paste(root, "0_categ", input, imagebase, sep = "/")

if(!file.exists(pathout))dir.create(pathout,showWarnings=FALSE, recursive=TRUE)

# GET METADATA
#################################################################
setwd(pathmeta)
metadata <- read.table(paste("metadata_", input, ".txt", sep = ""))
Band_nr <- metadata[3,2]
Study_area <- as.character(metadata[9,2])
Cloud_cover <- as.numeric(as.character(metadata[13,2]))
#y <- as.numeric(levels(metadata[13,2]))[metadata[13,2]]

# Set blue band
if(WriteLog){
  if(Band_nr==4)
  {blue <- 1
  } else if(Band_nr==8) 
  {blue <- 2
  } else if(Band_nr==5) 
  {blue <- 1
  } else 
  {blue <- "Unknown"
  }
}

setwd(pathin)
Files <- list.files(pattern = paste(input, ".*tif$", sep = ""), ignore.case = TRUE)
filename <- unlist(strsplit(Files[1], "[.]"))[1]

if(length(Files)!=1) next
# GDAL
x <- Files[1]
D <- readGDAL(x, band = blue)
threshold <- quantile(D@data$band1, 1-Cloud_cover, na.rm = TRUE)
if(threshold>10000){threshold <- 10000}
maxrefl <- max(D@data$band1, na.rm = TRUE)
pixelnu_image <- length(D@data$band1)

# DERIVE CLOUD MASK - REMOVE HIGH REFLECTANCE VALUES ACCORDING TO THE METADATA (0.052 - <CLOUDCOVER> - Cloud cover percentage)

# # Raster package
# R <- raster(D)
# #quantile(R@data@values, 1-Cloud_cover)
# maxrefl <- max(R@data@values, na.rm = TRUE)# R@data@max # maximum value
# pixelnu_image <- length(R@data@values) #length(R@data@values) # number of pixels
# #hist(R@data@values, breaks = 100, xlim = range(0, 15000), plot=TRUE, freq = FALSE, main = "Mali", xlab = "Reflectance*1000")
# #abline(v = max(R@data@values), col = "red", lwd = 3)
# 
# # mask non-cloudy values 
# maskR <- R
# maskR@data@values[maskR@data@values>quantile(maskR@data@values, 1-Cloud_cover, na.rm = TRUE)] <- NA
# #quantile(maskR@data@values, 1-Cloud_cover, na.rm = TRUE) # quantile 
# threshold <- max(maskR@data@values, na.rm = TRUE) # threshold for cloud detection
# pixelnu_clouds <- length(maskR@data@values[maskR@data@values>threshold]) # number of cloudy pixels
# percentage <- pixelnu_clouds/pixelnu_image
# 
# # export mask
# setwd(pathout)
# mask(R, maskR, filename=paste(filename, "_cloud_mask", sep =""), inverse=TRUE, format="GTiff", datatype='INT2S', overwrite=TRUE)

# create cloud mask
#mask <- SpatialGridDataFrame(D@grid, data.frame(band=array(0,nrow(D@data))), D@proj4string)
mask <- D
mask@data$band1[mask@data$band1<threshold] <- 1 # non-cloudy pixels
mask@data$band1[mask@data$band1>threshold] <- 0 # cloud pixels
#mask@data$band1[is.na(mask@data$band1)] <- 15 # NA pixels
#thres_check <- max(mask@data$band1[mask@data$band1 < 10000], na.rm = TRUE) # threshold for cloud detection
pixelnu_clouds <- length(mask@data$band1[mask@data$band1==1]) # number of cloudy pixels
percentage <- pixelnu_clouds/pixelnu_image

# export mask
setwd(pathout)
mask.tif<-create2GDAL(mask, drivername="GTiff",type="Byte", mvFlag = 255L)
saveDataset(mask.tif, paste(input, "_cloud_mask", ".tif", sep = ""))
GDAL.close(mask.tif)
#gdal_translate(paste(pathin, x, sep="/"), paste(filename, ".tif", sep = ""), of="GTiff", ot = "UInt16", verbose=TRUE)
#

if(WriteLog){
  #export report
  report <- rbind(filename, pixelnu_image, pixelnu_clouds, Cloud_cover*100, threshold/10000)
  rownames(report) <- c("Filename", "Input image pixels #", "Masked cloud pixels #", "Masked cloud pixels %", "Reflectance threshold band 2")
  write.table(report, file = paste(filename, "_report", ".txt", sep = ""), sep="\t", row.names = TRUE, col.names = FALSE)
  
  # HISTOGRAM
  png(filename = paste(filename, "_histogram", ".png", sep = ""), width = 1000, height = 500, units = "px", bg = "white", type = "cairo")
  #par(mar=c(0,0,0,0), mai=c(0,0,0,0))
  hist(D@data$band1, breaks = 300, xlim = range(0, 10000), plot=TRUE, freq = FALSE, main = filename, xlab = "Reflectance*1000")
  abline(v = threshold, col = "red", lwd = 3)
  legend("right", paste("Masked cloud pixels: ", as.character(report[4]), " %", sep =""), cex = 2)
  #curve(dnorm(x, mean=mean(D@data$band1, sd=sqrt(var(D@data$band1)))), col="darkblue", lwd=2, add=TRUE, yaxt="n")
  dev.off()
}

# CREATE MASKED-OUT FILE
#create
setwd(pathin)
#maskedout <- SpatialGridDataFrame(D@grid, data.frame(band=array(0,nrow(D@data))), D@proj4string)
maskedout <- readGDAL(x)

# for(i in 1:length(maskedout))
# {
#   if(maskedout@data$band1[i]>threshold) {
#     for(j in 1:length(maskedout@data))
#     {
#     maskedout@data$band1[j] <- 45001
#     }
#   }
# }

#maskedout@data$band1[maskedout@data$band1>threshold] <- 45001 # mask out cloud pixels
# maskremoval <- D
# maskremoval@data$band1[maskremoval@data$band1<threshold] <- 1 # non-cloudy pixels
# maskremoval@data$band1[maskremoval@data$band1>threshold] <- 0 # cloud pixels
# maskremoval_m <- as.matrix(maskremoval@data)

# maskedout_m <- as.matrix(maskedout@data$band1)
#maskedout@data$band1 <- mask@data$band1*maskedout@data$band1
setwd(pathout)
for(i in 1:length(maskedout@data)) {
  maskedout@data[i] <- mask@data$band1*maskedout@data[i]
  }

# export
setwd(pathout)
maskedout.tif<-create2GDAL(maskedout, drivername="GTiff",type="UInt16", mvFlag = 255L)
saveDataset(maskedout.tif, paste(filename, "_cloud", ".tif", sep = ""))
GDAL.close(maskedout.tif)

# convert to 0-255 (8bit)
#smin=0
#smax=255
#stre8b <- as.integer((E - min(E) ) * smax / ( max(E) - min(E) ) + smin)
#F <- as.matrix(stre8b)

# Otsu method - returns a threshold value to reduce a greyscale image to a binary file 
#F <- as.matrix(E) # convert to array
#grayimage <- channel(F,"gray")
#threshold <- otsu(F, range = c(0, 10000), levels = 16) # levels = 65536
#
