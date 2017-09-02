# RAPIDEYE IMAGE ARCHIVE METADATA REPORT
rm(list=ls(all=TRUE))

# INSTALL AND LOAD R PACKAGES
################################################################
ipak <- function(pkg){ # check to see if packages are installed. Install them if they are not, then load them into the R session.
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, repos="http://cran.rstudio.com/", dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}
packages <- c("XML", "rgdal", "proj4", "raster")
ipak(packages)

# DEFINE WORKING DIRECTORIES
#################################################################
root <- Sys.getenv('baseDG')
#root <- "/home/stratouliasd/stars/derived/RE_v8"
pathin <- paste(root, "0_categ", sep = "/")
#pathin <- "/home/stratouliasd/stars/acquired/RapidEye"
pathin_shp <- paste(root, "0_shapefiles", sep = "/")
pathout <- root
setwd(pathin)

# DEFINE PROJECTION SYSTEMS
#################################################################
lat_long <- "+proj=longlat +ellps=WGS84"
#UTM30N <- "+proj=utm +zone=30 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"


# SHAPEFILES OF STUDY AREAS - CONVERT TO LAT/LONG
#################################################################

shp_latlong <- function(shp){ 
  setwd(pathin_shp)
  readshp <- readOGR(dsn = path.expand(pathin_shp), layer = shp)
  readshp@proj4string@projargs
  transf_latlong <- spTransform(readshp, CRS(lat_long))
}

# READ SHAPEFILES
MALI_latlong <- shp_latlong("MALI10_10_NEW_BOX")
NIGERIA_latlong <- shp_latlong("KKMa3_Kofa_box")
TANZANIA_NJOMBE_latlong <- shp_latlong("TZ_Njombe20x20km")
TANZANIA_KILOSA_latlong <- shp_latlong("TZ_Kilosa20x20km")
TANZANIA_SAME_latlong <- shp_latlong("TZ_Same20x20km")
UGANDA_latlong <- shp_latlong("UG_Moroto20x20km")
BANGLADESH_latlong <- shp_latlong("BD_6_IrMASaT_Sites_DG_20150128_10by10km")
MALI_BAMAKO_CAMPUS_latlong <- shp_latlong("ML_Bamako_campus")
TANZANIA_MOROGORO_CAMPUS_latlong <- shp_latlong("TZ_Morogoro_campus")

# GET EXTENT OF SHAPEFILES IN LAT LONG
NIGERIA_latlong_extent <- extent(NIGERIA_latlong@bbox[1,1], NIGERIA_latlong@bbox[1,2], NIGERIA_latlong@bbox[2,1], NIGERIA_latlong@bbox[2,2]) 
MALI_latlong_extent <- extent(MALI_latlong@bbox[1,1], MALI_latlong@bbox[1,2], MALI_latlong@bbox[2,1], MALI_latlong@bbox[2,2])
BANGLADESH_latlong_extent <- extent(BANGLADESH_latlong@bbox[1,1], BANGLADESH_latlong@bbox[1,2], BANGLADESH_latlong@bbox[2,1], BANGLADESH_latlong@bbox[2,2]) 
TANZANIA_KILOSA_latlong_extent <- extent(TANZANIA_KILOSA_latlong@bbox[1,1], TANZANIA_KILOSA_latlong@bbox[1,2], TANZANIA_KILOSA_latlong@bbox[2,1], TANZANIA_KILOSA_latlong@bbox[2,2]) 
TANZANIA_NJOMBE_latlong_extent <- extent(TANZANIA_NJOMBE_latlong@bbox[1,1], TANZANIA_NJOMBE_latlong@bbox[1,2], TANZANIA_NJOMBE_latlong@bbox[2,1], TANZANIA_NJOMBE_latlong@bbox[2,2]) 
TANZANIA_SAME_latlong_extent <- extent(TANZANIA_SAME_latlong@bbox[1,1], TANZANIA_SAME_latlong@bbox[1,2], TANZANIA_SAME_latlong@bbox[2,1], TANZANIA_SAME_latlong@bbox[2,2]) 
UGANDA_latlong_extent <- extent(UGANDA_latlong@bbox[1,1], UGANDA_latlong@bbox[1,2], UGANDA_latlong@bbox[2,1], UGANDA_latlong@bbox[2,2]) 
MALI_BAMAKO_CAMPUS_latlong_extent <- extent(MALI_BAMAKO_CAMPUS_latlong@bbox[1,1], MALI_BAMAKO_CAMPUS_latlong@bbox[1,2], MALI_BAMAKO_CAMPUS_latlong@bbox[2,1], MALI_BAMAKO_CAMPUS_latlong@bbox[2,2]) 
TANZANIA_MOROGORO_CAMPUS_latlong_extent <- extent(TANZANIA_MOROGORO_CAMPUS_latlong@bbox[1,1], TANZANIA_MOROGORO_CAMPUS_latlong@bbox[1,2], TANZANIA_MOROGORO_CAMPUS_latlong@bbox[2,1], TANZANIA_MOROGORO_CAMPUS_latlong@bbox[2,2]) 

areadf <- data.frame("Class" = character(1), stringsAsFactors=FALSE)
dirdf <- data.frame("Class" = character(1), stringsAsFactors=FALSE)
cat_iddf <- data.frame("Class" = character(1), stringsAsFactors=FALSE)
study_areadf <- data.frame("Class" = character(1), stringsAsFactors=FALSE)
order_nudf <- data.frame("Class" = character(1), stringsAsFactors=FALSE)
datedf <- data.frame("Class" = character(1), stringsAsFactors=FALSE)
timedf <- data.frame("Class" = character(1), stringsAsFactors=FALSE)
satdf <- data.frame("Class" = character(1), stringsAsFactors=FALSE)
nadirdf <- data.frame("Class" = character(1), stringsAsFactors=FALSE)
sunazdf <- data.frame("Class" = character(1), stringsAsFactors=FALSE)
suneldf <- data.frame("Class" = character(1), stringsAsFactors=FALSE)
satazdf <- data.frame("Class" = character(1), stringsAsFactors=FALSE)
sateldf <- data.frame("Class" = character(1), stringsAsFactors=FALSE)
incidenceAngledf <- data.frame("Class" = character(1), stringsAsFactors=FALSE)
proc_leveldf <- data.frame("Class" = character(1), stringsAsFactors=FALSE)
cloudcoverdf <- data.frame("Class" = character(1), stringsAsFactors=FALSE)
demdf <- data.frame("Class" = character(1), stringsAsFactors=FALSE)
spresdf <- data.frame("Class" = character(1), stringsAsFactors=FALSE)
bandnudf <- data.frame("Class" = character(1), stringsAsFactors=FALSE)
ULXdf <- data.frame("Class" = character(1), stringsAsFactors=FALSE)
ULYdf <- data.frame("Class" = character(1), stringsAsFactors=FALSE)
URXdf <- data.frame("Class" = character(1), stringsAsFactors=FALSE)
URYdf <- data.frame("Class" = character(1), stringsAsFactors=FALSE)
LRXdf <- data.frame("Class" = character(1), stringsAsFactors=FALSE)
LRYdf <- data.frame("Class" = character(1), stringsAsFactors=FALSE)
LLXdf <- data.frame("Class" = character(1), stringsAsFactors=FALSE)
LLYdf <- data.frame("Class" = character(1), stringsAsFactors=FALSE)

# READ XML FILES
setwd(pathin)
RapidEye <- data.frame("Class" = character(0), stringsAsFactors=FALSE) 

xml <- list.files(pattern = '*metadata.xml', recursive = T, ignore.case = TRUE)

for (i in 1:length(xml)){

xml_split_slash <- unlist(strsplit(xml[i], "[/]"))
dir <- xml_split_slash[1]
l_slash <- length(xml_split_slash)
xml_split <- unlist(strsplit(xml_split_slash[l_slash], "[_]"))

if (xml_split[3] == "1B-NAC") {
  date <- substring(xml_split[1], 1, 10)
  substring(date, 5) <- "Y"
  substring(date, 8) <- "M"
  date <- paste(date, "D", sep = "")
  time <- substring(xml_split[1], 11, 17)
  time <- paste(substr(time, 2, 3), "hh", substr(time, 4, 5), "mm", substr(time, 6, 7), "ss", sep = "")
  sat <- xml_split[2]
  proc_level <- "1B"
  cat_id <- xml_split[4]
  order_nu <- xml_split[5]
 
  filename <- paste(proc_level, sat, date, order_nu, sep = "_")
  
  }	else if (xml_split[4] == "3A") {
    tile_id <- xml_split[1]
    date <- substring(xml_split[2], 1, 10)
    substring(date, 5) <- "Y"
    substring(date, 8) <- "M"
    date <- paste(date, "D", sep = "")
    time <- "not_defined"
    sat <- xml_split[3]
    proc_level <- "3A"
	cat_id <- "no_catalogue_id"
    order_nu <- xml_split[5]
    
    filename <- paste(proc_level, sat, date, order_nu, sep = "_")
    
    } else {

      sat <- "UNKNOWN"
      date <- "UNKNOWN"
	  proc_level <- "UNKNOWN"
	  cat_id <- "UNKNOWN"
      order_nu <- "UNKNOWN"
      
      filename <- paste(proc_level, sat, date, order_nu, sep = "_")
     }


# GET INFO FROM XML FILE
xmlfile <- xmlParse(xml[i])
xmlvalues <- xmlToList(xmlfile)
bandnu <- xmlvalues$resultOf$EarthObservationResult$product$ProductInformation$numBands
spresrow <- xmlvalues$resultOf$EarthObservationResult$product$ProductInformation$rowGsd
sprescol <- xmlvalues$resultOf$EarthObservationResult$product$ProductInformation$columnGsd
if (spresrow == sprescol) {
	spres <- spresrow
   }  else {
   	  spres <- "row and column GSD not equal"
	  }
dem <- xmlvalues$resultOf$EarthObservationResult$product$ProductInformation$elevationCorrectionApplied
sunaz <- xmlvalues$using$EarthObservationEquipment$acquisitionParameters$Acquisition$illuminationAzimuthAngle$text
sunel <- xmlvalues$using$EarthObservationEquipment$acquisitionParameters$Acquisition$illuminationElevationAngle$text
sataz <- xmlvalues$using$EarthObservationEquipment$acquisitionParameters$Acquisition$azimuthAngle$text

# satellite elevation is assumed 90degrees - incidence angle
incidenceAngle <- xmlvalues$using$EarthObservationEquipment$acquisitionParameters$Acquisition$incidenceAngle$text
satel <- as.character(90 - as.numeric(incidenceAngle))

nadir <- xmlvalues$using$EarthObservationEquipment$acquisitionParameters$Acquisition$spaceCraftViewAngle$text
cloudcover_reported <- xmlvalues$resultOf$EarthObservationResult$cloudCoverPercentage$text
cloudcover <- (as.numeric(cloudcover_reported))/100

ULX <- as.numeric(xmlvalues$target$Footprint$geographicLocation$topLeft$longitude)
ULY <- as.numeric(xmlvalues$target$Footprint$geographicLocation$topLeft$latitude)
URX <- as.numeric(xmlvalues$target$Footprint$geographicLocation$topRight$longitude)
URY <- as.numeric(xmlvalues$target$Footprint$geographicLocation$topRight$latitude)
LRX <- as.numeric(xmlvalues$target$Footprint$geographicLocation$bottomRight$longitude)
LRY <- as.numeric(xmlvalues$target$Footprint$geographicLocation$bottomRight$latitude)
LLX <- as.numeric(xmlvalues$target$Footprint$geographicLocation$bottomLeft$longitude)
LLY <- as.numeric(xmlvalues$target$Footprint$geographicLocation$bottomLeft$latitude)

# GET STUDY AREA VIA POINT IN POLYGON OPERATION
#################################################################

# FIND EXTENT OF THE RASTER (DELIVERY) IN LAT LONG
NWLAT <- as.numeric(xmlvalues$target$Footprint$geographicLocation$topLeft$latitude)
NWLONG <- as.numeric(xmlvalues$target$Footprint$geographicLocation$topLeft$longitude)
SELAT <- as.numeric(xmlvalues$target$Footprint$geographicLocation$bottomRight$latitude)
SELONG <- as.numeric(xmlvalues$target$Footprint$geographicLocation$bottomRight$longitude)
raster_latlong_extent <- extent(NWLONG, SELONG, SELAT, NWLAT)

# FIND WHICH SHAPEFILE INTERSECTS WITH THE EXTENT OF THE RASTER LAYER
if(!is.null(intersect(BANGLADESH_latlong_extent, raster_latlong_extent))){
  study_area <<- "BG"
  study_area_oid <<- 0000
} else if(!is.null(intersect(NIGERIA_latlong_extent, raster_latlong_extent))){
  study_area <<- "NG_Kofa"
  study_area_oid <<- 1005
} else if(!is.null(intersect(MALI_latlong_extent, raster_latlong_extent))){
  study_area <<- "ML_Sukumba"
  study_area_oid <<- 1000
} else if(!is.null(intersect(TANZANIA_KILOSA_latlong_extent, raster_latlong_extent))){
  study_area <<- "TZ_Kilosa"
  study_area_oid <<- 1015 
} else if(!is.null(intersect(TANZANIA_NJOMBE_latlong_extent, raster_latlong_extent))){
  study_area <<- "TZ_Njombe"
  study_area_oid <<- 1010 
} else if(!is.null(intersect(TANZANIA_SAME_latlong_extent, raster_latlong_extent))){
  study_area <<- "TZ_Same"
  study_area_oid <<- 1020 
} else if(!is.null(intersect(UGANDA_latlong_extent, raster_latlong_extent))){
  study_area <<- "UG_Moroto"
  study_area_oid <<- 1025 
} else if(!is.null(intersect(MALI_BAMAKO_CAMPUS_latlong_extent, raster_latlong_extent))){
  study_area <<- "ML_Bamako_campus"
  study_area_oid <<- 1045 
} else if(!is.null(intersect(TANZANIA_MOROGORO_CAMPUS_latlong_extent, raster_latlong_extent))){
  study_area <<- "TZ_Morogoro_campus"
  study_area_oid <<- 1050 
} else {
  study_area <<- "Unknown"
  study_area_oid <<- 9999
}

dirdf <- rbind(dirdf, dir)
cat_iddf <- rbind(cat_iddf, cat_id)
study_areadf <- rbind(study_areadf, study_area)
order_nudf <- rbind(order_nudf, order_nu)
datedf <- rbind(datedf, date)
timedf <- rbind(timedf, time)
satdf <- rbind(satdf, sat)
nadirdf <- rbind(nadirdf, nadir)
sunazdf <- rbind(sunazdf, sunaz)
suneldf <- rbind(suneldf,sunel)
satazdf <- rbind(satazdf, sataz)
sateldf <- rbind(sateldf, satel)
incidenceAngledf <- rbind(incidenceAngledf, incidenceAngle)
proc_leveldf <- rbind(proc_leveldf, proc_level)
cloudcoverdf <- rbind(cloudcoverdf, cloudcover)
demdf <- rbind(demdf, dem)
spresdf <- rbind(spresdf, spres)
bandnudf <- rbind(bandnudf, bandnu)
ULXdf <- rbind(ULXdf, ULX)
ULYdf <- rbind(ULYdf, ULY)
URXdf <- rbind(URXdf, URX)
URYdf <- rbind(URYdf, URY)
LRXdf <- rbind(LRXdf, LRX)
LRYdf <- rbind(LRYdf, LRY)
LLXdf <- rbind(LLXdf, LLX)
LLYdf <- rbind(LLYdf, LLY)

colbind <- cbind(dir, cat_id, study_area, order_nu, date, time, sat, "mul", nadir, sunaz, sunel, sataz, satel, incidenceAngle, "resampling", proc_level, cloudcover, "descriptor", dem, spres, bandnu, ULX, ULY, URX, URY, LRX, LRY, LLX, LLY)

RapidEye <- rbind(RapidEye, colbind)

}

colnames(RapidEye) <- c("Folder", "Catalog_id", "Area_description", "Sales_Order", "Acquisition_date", "Acquisition_start_time", "Satellite", "Band_id", "Nadir_angle", "Sun_azimuth", "Sun_elevation", "Satellite_azimuth", "Satellite_elevation", "Incidence_angle", "Resampling_method", "Level", "Cloud_cover", "Descriptor", "Elevation_Correction_Applied", "Ground_Sample_Distance", "Band_number", "UL_easting", "UL_northing", "UR_easting", "UR_northing", "LR_easting", "LR_northing", "LL_easting", "LL_northing")

setwd(pathout)
write.csv(RapidEye, file = paste("RapidEye_metadata_report_", format(Sys.time(), "%d-%b-%Y"), ".csv", sep = ""))

# 1B: RapidEye Basic Product - Radiometric and sensor corrections applied to the data. On-board spacecraft attitude and ephemeris applied to the data. 
# 3A: RapidEye Ortho Product - Radiometric, sensor and geometric corrections applied to the data. The product accuracy depends on the quality of the ground control and DEMs used. Product is processed as an individual 25 km by 25 km tile
# 3B: RapidEye Ortho Take Product - Large-scale orthorectified product based on RapidEye Image Takes. Multiple images over an AOI will be bundle adjusted together for accuracy purposes.