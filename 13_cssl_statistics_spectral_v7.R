# THIS FILE PRODUCES STATISTICS FROM THE IMAGES BASED ON THE FMUs

# REMOVE ALL VARIABLES FROM THE WORKSPACE
rm(list=ls(all=TRUE))

##### INSTALL AND LOAD R PACKAGES
ipak <- function(pkg){ # check to see if packages are installed. Install them if they are not, then load them into the R session.
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, repos="http://cran.rstudio.com/", dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}
packages <- c("rgdal", "raster", "rgeos", "sp", "RPostgreSQL", "moments", "wkb")
ipak(packages)

Server <- TRUE 						# run on server
EXPORT_CSV <- TRUE 					# export csv files
EXPORT_POSTGRESQL_table <- FALSE 	# export to POSTGRESQL as table

# DEFINE WORKING DIRECTORIES
if(Server==TRUE){
	root <- Sys.getenv('baseDG')
	#root <- "/home/stratouliasd/stars/derived/DG_v8" # FROM linux352.itc.utwente.nl
	Path_lib <- Sys.getenv('basescript')
}else{
	root <- "S:/derived/DG_v9"
	Path_lib <- "S:/derived/scripts/v9"
}

setwd(Path_lib)
source("scale_space_lib_v7.R")
# source("tree_measurement_lib_v9.R")
# sourceCpp("matcher_v5.cpp")
# sourceCpp("blobs_lib_v4.cpp")
# sourceCpp("tree_mask_from_blobs_v3.cpp")
#========================================================
# IMAGE INPUT
#========================================================
args <- commandArgs(trailingOnly = TRUE)
input <- Sys.getenv('input')
ms.imagefn <- paste(input, "_P001_MUL", ".tif", sep = "")
imagebase <- unlist(strsplit(ms.imagefn, "[.]"))[1]

pathin <- paste(root, "6_image_matching", input, sep = "/")
pathout <- paste(root, "9_cssl_statistics", sep = "/")
pathin_shp <- paste(root, "0_shapefiles_fmu", sep = "/")
pathmeta <- paste(root, "0_categ", input, imagebase, sep = "/")
path_tree_mask <- paste(root, "8_tree_mask", input, sep = "/")
path_cloud_mask <- paste(root, "7_cloud_mask", input, sep = "/")

# GET IMAGE METADATA
setwd(pathmeta)
metadata <- read.table(paste("metadata_", input, ".txt", sep = ""), stringsAsFactors=FALSE)
study_area <- as.character(unlist(metadata[9,2]))
date <- as.character(unlist(metadata[11,2]))
study_area_OID <- as.integer(unlist(metadata[14,2]))
year <- as.integer(substr(date, 1, 4))
erosion <- as.numeric(2.0) # it does not work at the moment, erosion is set in the query line 

# For the database: select years that the polygons are active
  if(study_area_OID == "BG"){
  year_OID <- 2000
  } else if(study_area_OID == 1005){ #NG_Kofa
  year_OID <- 2015
  } else if(study_area_OID == 1000){ # ML_Sukumba
  year_OID <- 2014
  } else if(study_area_OID == "TZ_Kilosa"){
  year_OID <- 2000
  } else if(study_area_OID == "TZ_Njombe"){
  year_OID <- 2000
  } else if(study_area_OID == "TZ_Same"){
  year_OID <- 2000
  } else if(study_area_OID == "UG_Moroto"){
  year_OID <- 2000
  } else {
  year_OID <- "Unknown"
  }

# if(study_area != "NG_Kofa")  quit(save = "no", status = 0, runLast = TRUE)

# RULES TO CONNECT  

if(TRUE){

	# WITH RPostgreSQL
	drv <- dbDriver("PostgreSQL") # load PostgreSQL driver
	con <- dbConnect(drv, host='linux352.itc.utwente.nl', port='5432', dbname='cssl', user='dimitris') # open a connection
	#dbListTables(con) # list existing tables

	
	# specific query with prepared statement
	query <- postgresqlExecStatement(con, "SELECT id AS id, ST_AsBinary(geom) AS geom from common.cs_fmu_studyarea_year($1,$2,2.0)", c(study_area_OID, year_OID))
	fmu <- fetch(query,n=-1)

	dbDisconnect(con)     # closes the connection
	dbUnloadDriver(drv)   # free all the resources on the driver

	EPSG = make_EPSG()
	p4s = EPSG[which(EPSG$code == 4326), "prj4"]

	row.names(fmu) = fmu$id

	ind <- which(!is.na(fmu$geom))
	fmu <- fmu[ind,]

	# fmu <- fmu[1:50,]
	
	for (i in seq(nrow(fmu))){
		if (i == 1) {
			#spTemp = readWKT(fmu$geom[i], fmu$id[i])
			# If the PROJ4 string has been set, use the following instead
			spTemp = readWKB(hex2raw(fmu$geom[i]), fmu$id[i], p4s)
		}
		else{
			spTemp = rbind(
			# spTemp, readWKT(fmu$geom[i], fmu$id[i])
			# If the PROJ4 string has been set, use the following instead
			spTemp, readWKB(hex2raw(fmu$geom[i]), fmu$id[i], p4s)
		)
	  }
	}

	# Create SpatialPolygonsDataFrame, drop WKT field from attributes
	vector_postgresql = SpatialPolygonsDataFrame(spTemp, fmu[-2])

	#writeOGR(vector_postgresql, ".", "test", driver="ESRI Shapefile")
}else{
	# Just read shapefile instead
	# Read FMU shapefile
	Path_fmu <- paste("/home/tolpekin/R","/Input/Nigeria",sep="")
	fn.polygons <- "nigeria_utm"
	Pol <- readOGR(dsn=Path_fmu, layer=fn.polygons)
	names(Pol)[1] <- "id"
	proj_vector <- Pol@proj4string@projargs
	n.pol <- nrow(Pol)

	vector_postgresql <- Pol
}

# writeOGR(vector_postgresql, "/home/stratouliasd/stars/derived/DG_v9/0_shapefiles_fmu", "Nigeria_fmu_extracted_from_CSSL", driver="ESRI Shapefile")

# ================================================================================
# raster
# ================================================================================
setwd(pathin)
image_list <- list.files(pattern = paste(input, ".*tif$", sep = ""), ignore.case = TRUE)

if(length(image_list)==0){
	stop("Image not found in the directory")
}

ms.imagefn <- image_list[1]
sensor <- metadata[2,2]
imagebase <- unlist(strsplit(ms.imagefn, "[.]"))
folder_name <- substr(ms.imagefn, 0, 12)

ms.imageinfo <- GDALinfo(paste(pathin,ms.imagefn,sep="/"), silent=TRUE)
ps.ms <- c(ms.imageinfo[["res.x"]],ms.imageinfo[["res.x"]])
N0.ms <- ms.imageinfo[["rows"]]
M0.ms <- ms.imageinfo[["columns"]]
Nb 	  <- ms.imageinfo[["bands"]]

proj_raster_ms <- attributes(ms.imageinfo)$projection
xyLL <- c(ms.imageinfo[["ll.x"]],ms.imageinfo[["ll.y"]])
ysign <- attributes(ms.imageinfo)$ysign #y sign

#pan.imageinfo <- GDALinfo(paste(Path_pan,pan.imagefn,sep="/"), silent=TRUE)
#ps.pan <- c(pan.imageinfo[["res.x"]],pan.imageinfo[["res.x"]])
# N0.pan <- pan.imageinfo[["rows"]]
# M0.pan <- pan.imageinfo[["columns"]]

# ================================================================================
# vector
# ================================================================================

# shapefile <- "Mali_field_2014_smoothed"
# vector <- readOGR(pathin_shp, layer = shapefile) # read vector as SpatialPolygonsDataFrame
# 
# # set vector's projection system the same as raster
# if (projection(raster) != projection(vector)) 
#   vector <- spTransform(vector, projection(raster))

#compareCRS(raster, vector)
#vector@proj4string
#raster@proj4string

#work.region.tif <- raster(work.region.tif)
#reprojected <- projectRaster(work.region.tif, newproj)
#test <- spTransform(vector, newproj)
#test <- spatial_sync_vector(vector, work.region.tif) # sync CRS of vector and raster datasets

# ================================================================================
# POSTGRESQL OBJECTS
# ================================================================================

# set spdfFinal's projection system the same as raster
 if (proj_raster_ms != projection(vector_postgresql)) 
  vector_postgresql <- spTransform(vector_postgresql, proj_raster_ms)

# ================================================================================
# Extract raster clip-polygons and statistics
# ================================================================================
setwd(pathout)

n_pol <- nrow(vector_postgresql)

Nc <- array(0,n_pol) #create null array
n_overhead <- 5
n_vi <- 9				# number of VIs
n_diag <- Nb*(Nb-1)/2

col_num <- n_overhead+ 2*n_vi + 3*Nb + 2*n_diag	#number of columns, calculate and update
# overhead + VIs(mean+st_dev) + {mean + variance + skewness} + {covariance + correlation}
fmu_stats <- data.frame(array(NA,c(n_pol,col_num)),stringsAsFactors = FALSE)

# Overhead
names_overhead <- c("image_identifier","satellite_sensor","target_of_interest_oid","number_of_all_pixels","number_of_contributing_pixels")
# VIs
names_vi <- c("ndvi","sd_ndvi","ndvi_green","sd_ndvi_green","evi","sd_evi","tcari","sd_tcari","nir_r","sd_nir_r","sarvi","sd_sarvi","savi","sd_savi","msavi","sd_msavi","wv2vi","sd_wv2vi")
names_mean <- array("",Nb)
names_var  <- array("",Nb)
names_skew <- array("",Nb)
for(k in 1:Nb){
	names_mean[k] <- paste("mean",k,sep="")
	names_var[k]  <- paste("variance",k,sep="")
	names_skew[k] <- paste("skewness",k,sep="")
}

cov_names <- array("", n_diag)
cor_names <- array("", n_diag)
counter <- 0
for(k in 1:(Nb-1)){
	for(l in (k+1):Nb){
		counter <- counter + 1
		cov_names[counter] <- paste("cov",k,l,sep="")
		cor_names[counter] <- paste("cor",k,l,sep="")
	}
}

names(fmu_stats) <- c(names_overhead,names_vi,names_mean,names_var,names_skew,cov_names,cor_names)

Pol <- vector_postgresql

#i.pol <- 3
for(i.pol in 1:n_pol){

	field_id <- vector_postgresql@data$id[[i.pol]]

	# read shapefiles marginal coordinates
	# In case of complex polygons (with holes)
	n_subpol <- length(vector_postgresql@polygons[[i.pol]])
	
	pol <- Pol@polygons[[i.pol]]@Polygons[[1]]
	xy_pol <- pol@coords
	xrl <- range(xy_pol[,1])
	yrl <- range(xy_pol[,2])

	# Add margins to ensure sufficiently large dataset
	margin <- 2*mean(ps.ms)
	xrl <- xrl + c(-1,1) * margin
	yrl <- yrl + c(-1,1) * margin
	
	ijr <- xy_to_rowcol(cbind(xrl,yrl),ms.imageinfo)

	#if((ijr[1]>M0.ms)&(ijr[2]<1)&(ijr[3]>N0.ms)&(ijr[4]<1)){ # Polygon is  completely outside the image
	if((ijr[1]>M0.ms)|(ijr[2]<=1)|(ijr[4]<=1)|(ijr[3]>N0.ms)){ # Polygon is  completely outside the image
		fmu_stats[i.pol,1:n_overhead] <- c(folder_name, sensor, field_id, 0, 0)
		next
	}

	if(ijr[1]<1) ijr[1] <- 1
	if(ijr[2]>M0.ms) ijr[2] <- M0.ms
	if(ijr[3]<1) ijr[3] <- 1
	if(ijr[4]>N0.ms) ijr[4] <- N0.ms
	
	Path_in <- pathin
	raster_fmu <- read_subset(ms.imagefn,ijr[1],ijr[2],ijr[3],ijr[4])

	# If all pixels are NA, skip
	sums <- rowSums(raster_fmu@data)
	ind <- which(!is.na(sums))
	if(length(ind)==0){
		fmu_stats[i.pol,1:n_overhead] <- c(folder_name, sensor, field_id, 0, 0)
		next
	}

	#check if polygon and subset raster intersect locally
	pol_extent <- extent(min(pol@coords[,1]), max(pol@coords[,1]), min(pol@coords[,2]), max(pol@coords[,2])) 
	intersect(pol_extent, raster_fmu)
	if(is.null(intersect(pol_extent, raster_fmu))){ # Polygon is  completely outside the subset image
		fmu_stats[i.pol,1:n_overhead] <- c(folder_name, sensor, field_id, 0, 0)
		next
	}
	
	# Read and apply tree & tree shadow mask
	setwd(path_tree_mask)
	Path_in <- path_tree_mask
	filename <- paste("Tree_mask_",input,".tif",sep="")
	if(!file.exists(filename)){
		stop("Tree mask not found")
	}
	# Careful! Assuming that tree mask raster has identical grid with the multispectral image 
	tree_mask <- read_subset(filename,ijr[1],ijr[2],ijr[3],ijr[4])
	# Test if the assumption is valid
	if(!all.equal(tree_mask@grid, raster_fmu@grid, tolerance = 1.0e-03)) stop("Something is wrong with tree mask")
	
	ind_mask <- which(tree_mask$band1>0)
	raster_fmu@data[ind_mask,] <- NA

	# Read and apply cloud mask
	setwd(path_cloud_mask)
	Path_in <- path_cloud_mask
	filename <- paste(input,"_cloud_mask",".tif",sep="")

	if(!file.exists(filename)){
		stop("Cloud mask not found")
	}
	# Careful! Assuming that tree mask raster has identical grid with the multispectral image 
	cloud_mask <- read_subset(filename,ijr[1],ijr[2],ijr[3],ijr[4])
	# Test if the assumption is valid
	if(!all.equal(cloud_mask@grid,raster_fmu@grid)) stop("Something is wrong with cloud mask")

	ind_mask <- which(cloud_mask$band1==0)
	raster_fmu@data[ind_mask,] <- NA

	# SUBSETTING BY POLYGON
	# Remove holes!
	xy <- coordinates(raster_fmu)
	contrib_pix <- array(0,nrow(raster_fmu))
	
	tmp.in <- point.in.polygon(xy[,1],xy[,2],xy_pol[,1],xy_pol[,2])
	contrib_pix[tmp.in==1] <- 1
	
	if(n_subpol > 1){
		for(i_sub in 2:n_subpol){
			if(Pol@polygons[[i.pol]]@Polygons[[i_sub]]@hole){
				h_pol <- Pol@polygons[[i.pol]]@Polygons[[i_sub]]
				xy_hole <- h_pol@coords
				tmp.in <- point.in.polygon(xy[,1],xy[,2],xy_hole[,1],xy_hole[,2])
				contrib_pix[tmp.in==1] <- 0
			}
		}
	}

	ind <- which(contrib_pix==1)
	n_pix <- length(ind)
	fmu_values <- raster_fmu@data[ind,]

	# Normalize to range [0,1]
	fmu_values <- fmu_values / 10000

	sums <- rowSums(fmu_values)
	ind <- which(is.na(sums))
	pix_NA <- length(ind) # number of NAs
	# Remove NAs
	fmu_values <- fmu_values[!is.na(sums),]
	sums <- rowSums(fmu_values)
	
	pix_zeros <- sum(sums==0L, na.rm = TRUE) # number of zero pixels
	pix_contr <- nrow(fmu_values) - pix_zeros # number of contributing pixels - not equal to zero

	# Remove zeros and negative values
	fmu_values <- fmu_values[sums>0,]
	
	# Compute VIs
	y <- fmu_values[,]
	# remove NA, 0
	
	# Arrange the bands for VI computation, depends on the sensor 
	if(Nb==4){
		nir <- 4
		red <- 3
		green <- 2
		blue <- 1
		#red_edge <- 0
		#nir2 <- 0
	}else{
	
		if(Nb==8){
			nir2 <- 8
			nir <- 7
			red_edge <- 6
			red <- 5
			green <- 3
			blue <- 2
		}
	}

	# NDVI
	i_col <- 6
	veg_index <- (y[,nir] - y[,red])/(y[,nir] + y[,red])
	fmu_stats[i.pol,i_col] <- mean(veg_index)
	fmu_stats[i.pol,i_col+1] <- sd(veg_index)
	i_col <- i_col + 2
	# NDVI green
	veg_index <- (y[,nir] - y[,green])/(y[,nir] + y[,green])
	fmu_stats[i.pol,i_col] <- mean(veg_index)
	fmu_stats[i.pol,i_col+1] <- sd(veg_index)
	i_col <- i_col + 2
	# EVI
	veg_index <- 2.5*(y[,nir]-y[,red])/(y[,nir]+6*y[,red]-7.5*y[,blue]+1)
	fmu_stats[i.pol,i_col] <- mean(veg_index)
	fmu_stats[i.pol,i_col+1] <- sd(veg_index)
	i_col <- i_col + 2
	# Red Edge: TCARI
	if(Nb==8){
		veg_index <- 3*((y[,red_edge]-y[,red])-0.2*(y[,red_edge]-y[,green])*(y[,red_edge]/y[,red]))
	} else veg_index <- NA
	fmu_stats[i.pol,i_col] <- mean(veg_index)
	fmu_stats[i.pol,i_col+1] <- sd(veg_index)
	i_col <- i_col + 2
	# NIR-R
	veg_index <- y[,nir] - y[,red]
	fmu_stats[i.pol,i_col] <- mean(veg_index)
	fmu_stats[i.pol,i_col+1] <- sd(veg_index)
	i_col <- i_col + 2
	# SARVI
	gama <- 1		# 
	L <- 0.5		# between -0.9 and +1.6
	r_rb <- y[,red] - gama*(y[,blue] - y[,red])
	veg_index <- (1+L)*(y[,nir]-r_rb)/(y[,nir] + r_rb+L)
	fmu_stats[i.pol,i_col] <- mean(veg_index)
	fmu_stats[i.pol,i_col+1] <- sd(veg_index)
	i_col <- i_col + 2
	#SAVI
	veg_index <- (1+L)*(y[,nir] - y[,red])/(y[,nir] + y[,red] + L)
	fmu_stats[i.pol,i_col] <- mean(veg_index)
	fmu_stats[i.pol,i_col+1] <- sd(veg_index)
	i_col <- i_col + 2
	# MSAVI
	veg_index <- 0.5*(2*y[,nir]+1-sqrt(((2*y[,nir]+1)^2)-8*(y[,nir]-y[,red])))
	fmu_stats[i.pol,i_col] <- mean(veg_index)
	fmu_stats[i.pol,i_col+1] <- sd(veg_index)
	i_col <- i_col + 2
	# WV2 VI
	if(Nb==8){
		veg_index <- (y[,nir2]-y[,red])/(y[,nir2]+y[,red])
	} else veg_index <- NA
	fmu_stats[i.pol,i_col] <- mean(veg_index)
	fmu_stats[i.pol,i_col+1] <- sd(veg_index)

	i_col <- n_overhead + 2*n_vi
	for(j in 1:Nb) fmu_stats[i.pol,i_col+j] <- mean(y[,j])
	i_col <- i_col + Nb
	for(j in 1:Nb) fmu_stats[i.pol,i_col+j] <- var(y[,j])
	i_col <- i_col + Nb
	for(j in 1:Nb) fmu_stats[i.pol,i_col+j] <- skewness(y[,j])
	i_col <- i_col + Nb

	covariance <- cov(fmu_values) # calculated on matrix bands~observations
	correlation <- cor(fmu_values) # calculated on matrix bands~observations
	
	cov_upper <- array(NA, n_diag)
	cor_upper <- array(NA, n_diag)
	
	counter <- 0
	if(Nb>1){
		for(k in 1:(Nb-1))
		for(l in (k+1):Nb){
			counter <- counter + 1
			cov_upper[counter] <- covariance[k,l]
			cor_upper[counter] <- correlation[k,l]
		}
	}

	fmu_stats[i.pol,(i_col+1):(i_col+n_diag)] <- cov_upper
	i_col <- i_col + n_diag
	fmu_stats[i.pol,(i_col+1):(i_col+n_diag)] <- cor_upper

	fmu_stats[i.pol,1] <- folder_name 
	fmu_stats[i.pol,2] <- sensor
	fmu_stats[i.pol,3] <- field_id
	fmu_stats[i.pol,4] <- n_pix
	fmu_stats[i.pol,5] <- pix_contr
	
	#fmu_stats[i.pol,1:n_overhead] <- c(folder_name, sensor, field_id, n_pix, pix_contr)
	
	# Round numeric values
	# for(n in 1:length(fmu_stats[1,]))
	# for(m in 1:length(fmu_stats[,1])){
		# if(class(fmu_stats[n,m])=="numeric"){
		# xround <- round(fmu_stats[n,m], digits = 2)
		# fmu_stats[n,m] <- xround
		# } else next 
	# }

}

fmu_stats[[3]] <- as.integer(fmu_stats[[3]])
fmu_stats[[4]] <- as.integer(fmu_stats[[4]])
fmu_stats[[5]] <- as.integer(fmu_stats[[5]])

setwd(pathout)
if(EXPORT_CSV==TRUE){
filename <- substr(imagebase[1], 1, nchar(imagebase[1])-21)
write.csv(fmu_stats, file = paste(filename, "_spectral", ".csv", sep = ""), row.names=FALSE)
}

#============================================================
# The END
#============================================================
