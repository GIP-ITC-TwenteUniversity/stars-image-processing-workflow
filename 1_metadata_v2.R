#================================================================
# IMAGE CRITICAL INFO TO FEED IN SEQUENCE OF R SCRIPTS
#================================================================
require("XML")
require("proj4")
require("rgdal")
require("rgeos")

require("RPostgreSQL")
require("wkb")

args <- commandArgs(trailingOnly = TRUE)
input <- args[1]
nr_strips <- args[2]
# input <- Sys.getenv("input")
# nr_strips <- as.numeric(Sys.getenv("nr_strips"))

stars_dir <- Sys.getenv("stars")
root <- Sys.getenv("baseDG")
base <- Sys.getenv("base")
#================================================================
# DEFINE WORKING DIRECTORIES
#================================================================
path_origin <- paste(stars_dir,"acquired/DG", input, sep="/")
pathin_readme <- paste(root, "0_categ", input, sep = "/")
pathin_shp <- paste(base, "shapefiles", sep = "/")
path_product <- paste(root,"0_categ",input,"GIS_FILES", sep="/")
#================================================================
# DEFINE PROJECTION SYSTEMS
#================================================================
# lat_long <- "+proj=longlat +ellps=WGS84"
# UTM30N <- "+proj=utm +zone=30 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"

# EPSG 4326
target_projection <- CRS("+init=epsg:4326")
#================================================================
# SHAPEFILES OF STUDY AREAS - CONVERT TO LAT/LONG
#================================================================
for(strip_id in 1:nr_strips){
	strip_str <- paste("P00", strip_id, sep="")
	subprod_str <- paste(input,strip_str,"MUL",sep="_")
	pathin <- paste(root, "0_categ", input, subprod_str, sep = "/")
	#================================================================
	# DG IMAGES GET BAND NUMBER AND STUDY AREA VIA POINT IN POLYGON OPERATION
	#================================================================
	setwd(path_origin)
	file <- list.files(pattern = paste("*P00",strip_id,".TIF$",sep=""), , recursive = TRUE)[1]
	image.info <- GDALinfo(file, returnStats=FALSE)
	projection <- attributes(image.info)$projection
	n_bands <- image.info[["bands"]]

	setwd(pathin_readme)  
	readmexml <- list.files(pattern = "*README.XML", recursive = FALSE)
	listoffiles <- lapply(readmexml, read.table, header=TRUE, colClasses="character", sep="\t")

	xmlfiles <- xmlTreeParse(readmexml)
	xmltop <- xmlRoot(xmlfiles)
	xmlvalues <- xmlSApply(xmltop, function(x) xmlSApply(x, xmlValue))
	df <- data.frame(t(xmlvalues),row.names=NULL)

	NWLAT  <- as.numeric(df$NWLAT)
	NWLONG <- as.numeric(df$NWLONG)
	SELAT  <- as.numeric(df$SELAT)
	SELONG <- as.numeric(df$SELONG)

	collection <- as.character(df$COLLECTIONSTART) 
	date <- substring(collection, 1, 10)
	substring(date, 5) <- "Y"
	substring(date, 8) <- "M"
	date <- paste(date, "D", sep = "")
	year <- as.numeric(substr(date, 1, 4))

	timestart <- substring(collection, 12, 19)
	timestart <- paste(substr(timestart, 1, 2), "hh", substr(timestart, 4, 5), "mm", substr(timestart, 7, 8), "ss", sep = "")

	setwd(pathin) 
	readimagexml <- list.files(pattern = paste("*P00",strip_id,".XML",sep=""), recursive = TRUE)
	readimagexml <- Filter(function(x) grepl("-M", x), readimagexml)
	listoffiles <- lapply(readimagexml, read.table, header=TRUE, colClasses="character", comment.char = "\t", sep="\t")

	xmlfiles <- xmlTreeParse(readimagexml)
	xmltop = xmlRoot(xmlfiles)
	xmlvalues <- xmlSApply(xmltop, function(x) xmlSApply(x, xmlValue))
	df <- data.frame(t(xmlvalues),row.names=NULL)

	sat_str <- as.character(xmlValue(xmlRoot(xmlfiles)[["IMD"]][["IMAGE"]][["SATID"]]))

	sat <- switch(sat_str,
		"WV03" = "WorldView-3",
		"WV02" = "WorldView-2",
		"WV01" = "WorldView-1",
		"QB02" = "QuickBird",
		"GE01" = "GeoEye-1",
		"Unknown"
	)

	sunaz <- as.character(xmlValue(xmlRoot(xmlfiles)[["IMD"]][["IMAGE"]][["MEANSUNAZ"]]))
	sunel <- as.character(xmlValue(xmlRoot(xmlfiles)[["IMD"]][["IMAGE"]][["MEANSUNEL"]]))
	sataz <- as.character(xmlValue(xmlRoot(xmlfiles)[["IMD"]][["IMAGE"]][["MEANSATAZ"]]))
	satel <- as.character(xmlValue(xmlRoot(xmlfiles)[["IMD"]][["IMAGE"]][["MEANSATEL"]]))
	nadir <- as.character(xmlValue(xmlRoot(xmlfiles)[["IMD"]][["IMAGE"]][["MEANOFFNADIRVIEWANGLE"]]))
	cloudcover <- as.character(xmlValue(xmlRoot(xmlfiles)[["IMD"]][["IMAGE"]][["CLOUDCOVER"]]))

	level <- as.character(xmlvalues[["IMD"]][["IMAGEDESCRIPTOR"]])

	# Identify, read and reproject the strip product footprint shapefile
	files <- list.files(path_product,pattern=paste(strip_str,"_PIXEL_SHAPE.SHP",sep=""),ignore.case=TRUE)
	filename <- grep("M[23][AD]S", files, ignore.case=TRUE,value=TRUE)
	# Delete the file extension
	tmp <- unlist(strsplit(filename,".",fixed=TRUE))
	filename <- tmp[-length(tmp)]
	strip_shape <- readOGR(dsn=path_product, layer=filename)
	strip_shape <- spTransform(strip_shape, target_projection)

	# Convert image footprint to text string
	xy <- strip_shape@polygons[[1]]@Polygons[[1]]@coords

	geom_str <- 'POLYGON(('
	for(i in 1:nrow(xy)){
		tmp <- paste(xy[i,1],' ',xy[i,2],sep='')
		if(i<nrow(xy)) tmp <- paste(tmp,', ',sep='') else tmp <- paste(tmp,'))',sep='')
		geom_str <- paste(geom_str,tmp,sep='')
	}

	# WITH RPostgreSQL
	drv <- dbDriver("PostgreSQL") # load PostgreSQL driver
	con <- dbConnect(drv, host='linux352.itc.utwente.nl', port='5432', dbname='cssl', user='dimitris')
	#dbListTables(con) # list existing tables

	query <- postgresqlExecStatement(con,"select s.oid as id, s.name as name, ST_AsBinary(st_transform(s.geometry,4326)) as geom, common.cs_utmzone_from_polygon(s.geometry) as utmzonecode from common.study_area as s where st_intersects(s.geometry,st_setsrid(st_geomfromtext($1),4326)) and $2 between s.year_start and s.year_end",c(geom_str,year))

	st_area_info_db <- fetch(query,n=-1)

	dbDisconnect(con)     # closes the connection
	dbUnloadDriver(drv)   # free all the resources on the driver

	# How many study areas intersect with the image?
	nr_rel_areas <- nrow(st_area_info_db)
	oids		 <- st_area_info_db$id

	# st_area_polygon <- SpatialPolygonsDataFrame(readWKB(hex2raw(st_area_info_db$geom), 1:nr_rel_areas, target_projection), data.frame(oid=oids))
	# writeOGR(st_area_polygon, ".", "test", driver="ESRI Shapefile")
	
	# Number of fixed columns
	nr_fixed <- 14+1
	metadata <- data.frame(array(NA,c(1,nr_fixed+nr_rel_areas)))

	colnames(metadata)[1:nr_fixed] <- c("Folder", "Strip", "Satellite", "Number_of_bands","Polygon_as_text", "Sun_azimuth", "Sun_elevation", "Satellite_azimuth", "Satellite_elevation", "Off_nadir_angle", "Processing_level", "Acquisition_date", "Acquisition_time", "Cloud_cover", "Nr_study_areas")

	if(nr_rel_areas>0) colnames(metadata)[nr_fixed+1:nr_rel_areas] <- paste("OID_study_area_",1:nr_rel_areas,sep="")

	metadata$Folder					<- input
	metadata$Strip					<- strip_str
	metadata$Satellite				<- sat
	metadata$Number_of_bands		<- n_bands
	metadata$Polygon_as_text		<- paste("\"",geom_str,"\"",sep="")
	metadata$Sun_azimuth			<- sunaz
	metadata$Sun_elevation			<- sunel
	metadata$Satellite_azimuth		<- sataz
	metadata$Satellite_elevation	<- satel
	metadata$Off_nadir_angle		<- nadir
	metadata$Processing_level		<- level
	metadata$Acquisition_date		<- date
	metadata$Acquisition_time		<- timestart
	metadata$Cloud_cover			<- cloudcover
	metadata$Nr_study_areas			<- nr_rel_areas

	if(nr_rel_areas>0) metadata[1,nr_fixed+1:nr_rel_areas] <- oids[1:nr_rel_areas]

	write.csv(metadata, file = paste("metadata_", input, ".csv", sep = ""),row.names = FALSE, quote = FALSE)
}
