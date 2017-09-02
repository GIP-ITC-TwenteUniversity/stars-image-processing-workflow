### STATISTICS VISUALIZATION
rm(list=ls(all=TRUE))

#install and load R packages
ipak <- function(pkg){ # check to see if packages are installed. Install them if they are not, then load them into the R session.
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}
packages <- c("lattice", "ggplot2", "reshape", "rgdal", "car")
ipak(packages)

# SET WORKING DIRECTORIES
Server <- FALSE

if(Server){
  root <- Sys.getenv('baseDG')
  pathin <- paste(root, "9_cssl_statistics", sep="/")
  pathout <- paste(root, "10_statistics_analysis", sep="/")
  pathin_shp <- paste(root, "0_shapefiles_fmu", sep="/")
} else{
    root <- "C:/DELL_backup_15_June_2016/STARS/statistics_visualization" #  D:/STARS/statistics_visualization
    pathin <- paste(root, "9_cssl_statistics", sep="/")
    pathout <- paste(root, "10_statistics_analysis", sep="/")
    pathin_shp <- paste(root, "0_shapefiles_fmu", sep="/")
  }

# SET CRITERIA
study_area <- "ML_Sukumba" # "ML_Sukumba" | "NG_Kofa" !!!!! Mali with 2014 data, Nigeria with 2015 data!!!!
year <- 2014 # 2014 | 2015
# colnames(df)[6:23] 
vegindices <- c("NDVI", "sd_NDVI", "NDVI_green", "sd_NDVI_green", "EVI", "sd_EVI", "TCARI", "sd_TCARI", "NIR.R", "sd_NIR.R", "SARVI", "sd_SARVI", "SAVI", "sd_SAVI", "MSAVI", "sd_MSAVI", "WV2VI", "sd_WV2VI")
vegindex <- "TCARI"
#oid <- 1015
#crop <- "all" # all | 

# READ STATISTIC FILES
setwd(pathin)
files <- list.files(pattern = ".*csv$", ignore.case = TRUE)
files <- subset(files, grepl("spectral", files))
files <- subset(files, grepl(study_area, files))
files <- subset(files, grepl(year, files))

#  assign(substr(files[i], 28, 38), read.csv(files[i]))

# READ IN-SITU INFO FROM SHAPEFILE
mali2015full_buffer2m <- readOGR(dsn = pathin_shp, layer = 'mali2015full_buffer2m')
mali2015full_buffer2m <- mali2015full_buffer2m@data
mali_farmingpractice <- read.csv(paste(pathin_shp, "mali_farmingpractice_comma_delimited.csv", sep = "/"), sep = ",", quote = "", colClasses = "character")
insitu <- merge(mali2015full_buffer2m, mali_farmingpractice, by.x = "FIELD_ID", by.y = "numero_champ", all = TRUE)   
insitu  <- insitu[, c("ID", "FIELD_ID", "CROPTYPE", "date_labour", "date_semis")]
colnames(insitu) <- c("OID", "Field_id", "Crop_type", "Ploughing_date", "Sowing_date")
croptype <- as.character(unique(insitu[["Crop_type"]]))


################################################################################################
# EXTRACT VEGETATION INDICES TEMPORAL CURVES (PER FMU)
################################################################################################

nrow <- length(files)
ncol <- length(mali2015full_buffer2m[,1])

for (k in 1:length(vegindices)) {
vegindex <- vegindices[k]
  
  fmu_stats <- data.frame(array(NA,c(nrow,ncol)),stringsAsFactors = FALSE)
  datedf  <- data.frame(array(NA,c(nrow,1)),stringsAsFactors = FALSE)
  fmudf <- data.frame(array(NA,c(ncol,2)),stringsAsFactors = FALSE)
  names(datedf) <- "date"
  names(fmudf) <- c("variable", "sowing_date")
  for (j in 1:ncol) {
    for (i in 1:nrow) {
      csvfile <- read.csv(files[i])
      csvfile <- merge(insitu, csvfile, by.x = "OID", by.y = "target_of_interest_oid", all = TRUE)  # MATCH FIELD_ID (field data ID) AND OID (database ID)   
      df <- csvfile[j,]
      date <- strsplit(files[i], split = "_")[[1]][5]
      #date <- substr(files[i], 28, 38)
      date <- paste(substring(date, 1, 4), substring(date, 6, 7), substring(date, 9, 10), sep = "/")
      # convert date to julian date
      date <- strptime(date, "%Y/%m/%d")$yday+1
      if (j==1) {
        if (date %in% datedf$date) {
          date <- date+1000
          row.names(fmu_stats)[i] <- date
        }
      }
      datedf[i,1] <- date
      fmu <- paste(df$OID, df$Crop_type, sep =" - ")
      fmudf[j,] <- c(fmu, strptime(insitu$Sowing_date[j], "%d/%m/%Y")$yday+1)
      vegindex_value <- df[[vegindex]]
      
      names(fmu_stats)[j] <- fmu
      fmu_stats[i,j] <- vegindex_value    
    }
  }

  fmudf$sowing_date <- as.numeric(fmudf$sowing_date)
  fmu <- cbind(datedf, fmu_stats)

  # PLOT WITH SOWING DATES
  dfmelt <- melt(fmu,  id.vars = 'date', variable.name = 'variable')
  #dfmelt <- dfmelt[grepl("Cotton", dfmelt$variable), ]
  fmu_names <- unique(substring(dfmelt$variable, 1, 4))
  g <- ggplot(dfmelt, aes(date,value, group =1)) + geom_smooth() + geom_point() + geom_vline(aes(xintercept=sowing_date), data=fmudf, colour = "red")
  g + facet_wrap(~variable) + 
    ggtitle(paste(vegindex, "for individual fmu in", study_area, year, sep = " ")) +
    xlab("Julian day") +
    ylab(vegindex)
  ggsave(paste(pathout, "/", "per_fmu", "/", vegindex, "_for_individual_fmu_", study_area, year, ".png", sep = ""), width = 18, height = 9)

}

################################################################################################
# EXTRACT VEGETATION INDICES TEMPORAL CURVES (PER CROP_TYPE)
################################################################################################

nrow <- length(files)
ncol <- length(croptype)

for (k in 1:length(vegindices)) {
vegindex <- vegindices[k]
  
  fmu_stats <- data.frame(array(NA,c(nrow,ncol)),stringsAsFactors = FALSE)
  datedf  <- data.frame(array(NA,c(nrow,1)),stringsAsFactors = FALSE)
  names(datedf) <- "date"
  for (j in 1:ncol) {
    for (i in 1:nrow) {
      csvfile <- read.csv(files[i])
      csvfile[["weighted_value"]] <- csvfile[[vegindex]] * csvfile$number_of_contributing_pixels
      csvfile <- merge(insitu, csvfile, by.x = "OID", by.y = "target_of_interest_oid", all = TRUE)  # MATCH FIELD_ID (field data ID) AND OID (database ID)   
      csvfile_croptype <- csvfile[csvfile["Crop_type"] == croptype[j], ]
      vegindex_value <- sum(csvfile_croptype$weighted_value, na.rm=TRUE)/sum(csvfile_croptype$number_of_contributing_pixels, na.rm=TRUE)
      #sd <- aggregate(weighted_value ~ Crop_type, csvfile, sd, na.rm = TRUE)
      #crop_stats <- cbind(meanval, sd$weighted_value)
      #colnames(crop_stats) <- c("Croptype", "mean", "sd")
   
      date <- strsplit(files[i], split = "_")[[1]][5]
      #date <- substr(files[i], 28, 38)
      date <- paste(substring(date, 1, 4), substring(date, 6, 7), substring(date, 9, 10), sep = "/")
      # convert date to julian date
      date <- strptime(date, "%Y/%m/%d")$yday+1
      if (j==1) {
        if (date %in% datedf$date) {
          date <- date+1000
          row.names(fmu_stats)[i] <- date
        }
      }
      datedf[i,1] <- date
      names(fmu_stats)[j] <- croptype[j]
      fmu_stats[i,j] <- vegindex_value 
    }
  }
  
  fmu <- cbind(datedf, fmu_stats)
  
  # INDIVIDUAL PLOTS 
  dfmelt <- melt(fmu,  id.vars = 'date', variable.name = 'variable')
  #levels(dfmelt$variable) <- c("S", "Ve", "Vi")
  #ggplot(dfmelt, aes(date,value, group =1)) + geom_point(aes(colour = variable))
  g <- ggplot(dfmelt[which(dfmelt$value != 0),], aes(date,value, group =1)) + geom_smooth() + geom_point()
  g + facet_wrap(~variable) + 
    ggtitle(paste(vegindex, "per crop type in", study_area, year, sep = " ")) +
    xlab("Julian day") +
    ylab(vegindex)
  ggsave(paste(pathout, "/", "per_crop_type", "/", vegindex, "_for_crop_type_", study_area, year, ".png", sep = ""), width = 18, height = 9)

  # # CUMMULATIVE PLOTS
  # library("reshape2")
  # dfmelt <- melt(fmu, id.vars="date", value.name="value", variable.name="Year")
  # g <- ggplot(dfmelt, aes(date,value, group = variable, colour = variable)) + geom_smooth(se = FALSE) + geom_point()
  # g + 
  #   ggtitle(paste(vegindex, "per crop type in", study_area, year, sep = " ")) +
  #   xlab("Julian day") +
  #   ylab(vegindex) +
  #   theme(legend.text=element_text(size=22)) +
  #   theme(legend.title=element_text(size=28))
  # ggsave(paste(pathout, "/", "per_crop_type", "/", vegindex, "_for_crop_type_", "cummulative_", study_area, year, ".png", sep = ""), width = 18, height = 9)
  # 
}

################################################################################################
# TEST 1 ------------ EXTRACT SPECTRO-TEMPORAL SURFACES PER CROP TYPE
################################################################################################

nrow <- length(files)
ncol <- 8

fmu <- data.frame(array(NA,c(nrow,ncol+1)),stringsAsFactors = FALSE)
datedf  <- data.frame(array(NA,c(nrow,1)),stringsAsFactors = FALSE)
names(datedf) <- "date"
wavelength <- c(300, 400, 500, 600, 700, 750, 800, 900)
for (j in 1:ncol) {
  for (i in 1:nrow) {
    csvfile <- read.csv(files[i])
    csvfile <- csvfile[c("target_of_interest_oid", "number_of_contributing_pixels", "mean1", "mean2", "mean3", "mean4", "mean5", "mean6", "mean7", "mean8")]
    for (n in 3:(ncol+2)) {    
      csvfile[[n]] <- csvfile[[n]] * csvfile$number_of_contributing_pixels
    }
    #csvfile[["weighted_value"]] <- csvfile[[vegindex]] * csvfile$number_of_contributing_pixels
    csvfile <- merge(insitu, csvfile, by.x = "OID", by.y = "target_of_interest_oid", all = TRUE)  # MATCH FIELD_ID (field data ID) AND OID (database ID)   
    csvfile_croptype <- csvfile[csvfile["Crop_type"] == croptype[j], ]
    
    band1 <- sum(csvfile_croptype[[7]], na.rm=TRUE)/sum(csvfile_croptype$number_of_contributing_pixels, na.rm=TRUE)
    band2 <- sum(csvfile_croptype[[8]], na.rm=TRUE)/sum(csvfile_croptype$number_of_contributing_pixels, na.rm=TRUE)
    band3 <- sum(csvfile_croptype[[9]], na.rm=TRUE)/sum(csvfile_croptype$number_of_contributing_pixels, na.rm=TRUE)
    band4 <- sum(csvfile_croptype[[10]], na.rm=TRUE)/sum(csvfile_croptype$number_of_contributing_pixels, na.rm=TRUE)
    band5 <- sum(csvfile_croptype[[11]], na.rm=TRUE)/sum(csvfile_croptype$number_of_contributing_pixels, na.rm=TRUE)
    band6 <- sum(csvfile_croptype[[12]], na.rm=TRUE)/sum(csvfile_croptype$number_of_contributing_pixels, na.rm=TRUE)
    band7 <- sum(csvfile_croptype[[13]], na.rm=TRUE)/sum(csvfile_croptype$number_of_contributing_pixels, na.rm=TRUE)
    band8 <- sum(csvfile_croptype[[14]], na.rm=TRUE)/sum(csvfile_croptype$number_of_contributing_pixels, na.rm=TRUE)

    #sd <- aggregate(weighted_value ~ Crop_type, csvfile, sd, na.rm = TRUE)
    #crop_stats <- cbind(meanval, sd$weighted_value)
    #colnames(crop_stats) <- c("Croptype", "mean", "sd")
    
    date <- strsplit(files[i], split = "_")[[1]][5]
    #date <- substr(files[i], 28, 38)
    date <- paste(substring(date, 1, 4), substring(date, 6, 7), substring(date, 9, 10), sep = "/")
    # convert date to julian date
    date <- strptime(date, "%Y/%m/%d")$yday+1
    if (j==1) {
      if (date %in% datedf$date) {
        date <- date+1000
        row.names(fmu_stats)[i] <- date
      }
    }
    
    
    fmu[i,] <- cbind(date, band1, band2, band3, band4, band5, band6, band7, band8)

  }
}

dfmelt <- melt(fmu,  id.vars = 'date', variable.name = 'variable')
colnames(dfmelt) <- c("date", "crop type", vegindex)

data_x <- dfmelt$date
data_y <- dfmelt$`crop type`
data_z <- dfmelt$NDVI

scatter3d(x = data_x, y = data_y, z = data_z, point.col = "blue", grid = FALSE, fit = "smooth")
rgl.snapshot(filename = "3d_spectro_temporal_plane.png")
setwd(paste(output, "per_crop_spectro_temporal_planes", sep = "/"))
rgl.postscript("plot.pdf",fmt="pdf")
  


################################################################################################
# EXTRACT SPECTRO-TEMPORAL SURFACES PER CROP TYPE
################################################################################################

nrow <- length(files)*8 #add test for 4 or 8 bands
ncol <- 3

fmu <- data.frame(array(NA,c(nrow,ncol)),stringsAsFactors = FALSE)
datedf  <- data.frame(array(NA,c(nrow/8,1)),stringsAsFactors = FALSE)
names(datedf) <- "date"
bands <- c("band1", "band2", "band3", "band4", "band5", "band6", "band7", "band8")
wavelength <- c(300, 400, 500, 600, 700, 750, 800, 900)
  for (i in 1:nrow) {
    csvfile <- read.csv(files[i])
    csvfile <- csvfile[c("target_of_interest_oid", "number_of_contributing_pixels", "mean1", "mean2", "mean3", "mean4", "mean5", "mean6", "mean7", "mean8")]
    for (n in 3:length(csvfile)) {    
      csvfile[[n]] <- csvfile[[n]] * csvfile$number_of_contributing_pixels
    }
    #csvfile[["weighted_value"]] <- csvfile[[vegindex]] * csvfile$number_of_contributing_pixels
    csvfile <- merge(insitu, csvfile, by.x = "OID", by.y = "target_of_interest_oid", all = TRUE)  # MATCH FIELD_ID (field data ID) AND OID (database ID)   
    csvfile_croptype <- csvfile[csvfile["Crop_type"] == croptype[j], ] #merge per crop type
    
    for (k in 1:8) {
      bands[k] <- sum(csvfile_croptype[[(k+6)]], na.rm=TRUE)/sum(csvfile_croptype$number_of_contributing_pixels, na.rm=TRUE)
      sample <- cbind(date, wavelength[k], band1) # populate the dataframe here!!!!!!!!!
    }
    
    band2 <- sum(csvfile_croptype[[8]], na.rm=TRUE)/sum(csvfile_croptype$number_of_contributing_pixels, na.rm=TRUE)
    band3 <- sum(csvfile_croptype[[9]], na.rm=TRUE)/sum(csvfile_croptype$number_of_contributing_pixels, na.rm=TRUE)
    band4 <- sum(csvfile_croptype[[10]], na.rm=TRUE)/sum(csvfile_croptype$number_of_contributing_pixels, na.rm=TRUE)
    band5 <- sum(csvfile_croptype[[11]], na.rm=TRUE)/sum(csvfile_croptype$number_of_contributing_pixels, na.rm=TRUE)
    band6 <- sum(csvfile_croptype[[12]], na.rm=TRUE)/sum(csvfile_croptype$number_of_contributing_pixels, na.rm=TRUE)
    band7 <- sum(csvfile_croptype[[13]], na.rm=TRUE)/sum(csvfile_croptype$number_of_contributing_pixels, na.rm=TRUE)
    band8 <- sum(csvfile_croptype[[14]], na.rm=TRUE)/sum(csvfile_croptype$number_of_contributing_pixels, na.rm=TRUE)
    
    #sd <- aggregate(weighted_value ~ Crop_type, csvfile, sd, na.rm = TRUE)
    #crop_stats <- cbind(meanval, sd$weighted_value)
    #colnames(crop_stats) <- c("Croptype", "mean", "sd")
    
    date <- strsplit(files[i], split = "_")[[1]][5]
    #date <- substr(files[i], 28, 38)
    date <- paste(substring(date, 1, 4), substring(date, 6, 7), substring(date, 9, 10), sep = "/")
    # convert date to julian date
    date <- strptime(date, "%Y/%m/%d")$yday+1
    if (j==1) {
      if (date %in% datedf$date) {
        date <- date+1000
        row.names(fmu_stats)[i] <- date
      }
    }
    
    
    fmu[i,] <- cbind(date, band1, band2, band3, band4, band5, band6, band7, band8)
    
  }
}

dfmelt <- melt(fmu,  id.vars = 'date', variable.name = 'variable')
colnames(dfmelt) <- c("date", "crop type", vegindex)

data_x <- dfmelt$date
data_y <- dfmelt$`crop type`
data_z <- dfmelt$NDVI

scatter3d(x = data_x, y = data_y, z = data_z, point.col = "blue", grid = FALSE, fit = "smooth")
rgl.snapshot(filename = "3d_spectro_temporal_plane.png")
setwd(paste(output, "per_crop_spectro_temporal_planes", sep = "/"))

# expected results:
#rice monocroping 
#maize similar to sorghum (the latter is delayed)
#Millet
#Spice crops
#Soybeans