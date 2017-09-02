require(rgdal)
require(rgeos)
require("RPostgreSQL")
require("wkb")

args <- commandArgs(trailingOnly = TRUE)
input <- args[1]
nr_strips <- args[2] 

# input <- Sys.getenv("input")
# nr_strips <- as.numeric(Sys.getenv("nr_strips"))

# EPSG 4326
target_projection <- CRS("+init=epsg:4326")

root <- Sys.getenv("baseDG")
base <- Sys.getenv("base")
stars_dir <- Sys.getenv("stars")

path_area <- paste(base,"shapefiles", sep="/")
path_product <- paste(root,"0_categ",input,"GIS_FILES", sep="/")
path_out <- paste(root,"2_orthorectification_gdal",input, sep="/")
if(!file.exists(path_out))dir.create(path_out,showWarnings=FALSE, recursive=TRUE)
path_out_pan <- paste(root,"2_orthorectification_gdal_pan",input, sep="/")
if(!file.exists(path_out_pan))dir.create(path_out_pan,showWarnings=FALSE, recursive=TRUE)

# Determine number of strips and filenames
files <- list.files(path_product,pattern="P00[0-9]_PIXEL_SHAPE.SHP",ignore.case=TRUE)
ind <- grep("M2AS", files, ignore.case=TRUE)
ms_files <- files[ind]

ind <- grep("P2AS", files, ignore.case=TRUE)
pan_files <- files[ind]

# copy the TIL, IMD & RPB files
path_dest	<- paste(root, "1_atcor_6s", input, sep="/")

# Missing value flag
miss_val <- 2^16 - 1

for(strip_id in 1:nr_strips){
	strip_str <- paste("P00", strip_id, sep="")
	subprod_str <- paste(input,strip_str,"MUL",sep="_")

	path_src 	<- paste(root, "0_categ", input, subprod_str, sep="/")

	til_files <- list.files(path_src, pattern=".TIL", full.names = TRUE, ignore.case=TRUE)
	file.copy(til_files, path_dest, overwrite = TRUE, recursive = FALSE, copy.mode = TRUE, copy.date = FALSE)

	imd_files <- list.files(path_src, pattern=".IMD", full.names = TRUE, ignore.case=TRUE)
	file.copy(imd_files, path_dest, overwrite = TRUE, recursive = FALSE, copy.mode = TRUE, copy.date = FALSE)

	rpb_files <- list.files(path_src, pattern=".RPB", full.names = TRUE, ignore.case=TRUE)
	file.copy(rpb_files, path_dest, overwrite = TRUE, recursive = FALSE, copy.mode = TRUE, copy.date = FALSE)

	# Use relevant study area (determine from Metadata)
	path_meta <- paste(root, "0_categ", input, subprod_str, sep = "/")
	setwd(path_meta)
	metadata <- read.csv(paste("metadata_", input, ".csv", sep = ""), header=TRUE, stringsAsFactors=FALSE)

	nr_fixed <- ncol(metadata) - metadata$Nr_study_areas

	if(metadata$Nr_study_areas==0) stop("Error: Study area is unknown!")

	for(i_st_area in 1:metadata$Nr_study_areas){
		oid <- metadata[1,nr_fixed+i_st_area]

		# WITH RPostgreSQL
		drv <- dbDriver("PostgreSQL") # load PostgreSQL driver
		con <- dbConnect(drv, host='linux352.itc.utwente.nl', port='5432', dbname='cssl', user='dimitris')
		#dbListTables(con) # list existing tables

		query <- postgresqlExecStatement(con,"select s.name as name, s.dem_name as dem_name, ST_AsBinary(st_transform(s.geometry,4326)) as geom, common.cs_utmzone_from_polygon(s.geometry) as utmzonecode from common.study_area as s where oid=$1",oid)
		# query <- postgresqlExecStatement(con,"select * from common.study_area as s where oid=$1",oid)
		st_area_info_db <- fetch(query,n=-1)

		dbDisconnect(con)     # closes the connection
		dbUnloadDriver(drv)   # free all the resources on the driver

		study_area_shape <- SpatialPolygonsDataFrame(readWKB(hex2raw(st_area_info_db$geom), 1, target_projection), data.frame(oid=oid))

		product_type <- metadata$Processing_level

		tmp <- ms_files[strip_id]
		tmp <- unlist(strsplit(tmp,".",fixed=TRUE))
		filename <- tmp[-length(tmp)]

		strip_shape <- readOGR(dsn = path_product, layer = filename)

		# Check if this is required: ensure the shapefiles have the same projection
		if(!identical(CRS(study_area_shape@proj4string@projargs),CRS(strip_shape@proj4string@projargs))){
			strip_shape <- spTransform(strip_shape, study_area_shape@proj4string)
		}

		subset_pol <- gIntersection(strip_shape,study_area_shape)

		# convert from SpatialPolygons to SpatialPolygonsDataFrame
		datafr <- strip_shape@data
		row.names(datafr) <- row.names(subset_pol)
		subset_shape <- SpatialPolygonsDataFrame(subset_pol, datafr)

		tmp <- unlist(strsplit(filename,"_"))
		tmp <- tmp[length(tmp)-2]
		name_out <- paste("footprint_","ms_",tmp,"_st_area_",oid,sep="")
		writeOGR(subset_shape, path_out, name_out, driver="ESRI Shapefile", overwrite_layer=TRUE)

		# Panchromatic footprint
		tmp <- pan_files[strip_id]
		tmp <- unlist(strsplit(tmp,".",fixed=TRUE))
		filename <- tmp[-length(tmp)]

		strip_shape <- readOGR(dsn = path_product, layer = filename)

		# Check if this is required: ensure the shapefiles have the same projection
		if(!identical(CRS(study_area_shape@proj4string@projargs),CRS(strip_shape@proj4string@projargs))){
			strip_shape <- spTransform(strip_shape, study_area_shape@proj4string)
		}

		subset_pol <- gIntersection(strip_shape,study_area_shape)

		# convert from SpatialPolygons to SpatialPolygonsDataFrame
		datafr <- strip_shape@data
		row.names(datafr) <- row.names(subset_pol)
		subset_shape <- SpatialPolygonsDataFrame(subset_pol, datafr)

		tmp <- unlist(strsplit(filename,"_"))
		tmp <- tmp[length(tmp)-2]
		name_out <- paste("footprint_","pan_",tmp,"_st_area_",oid,sep="")
		writeOGR(subset_shape, path_out_pan, name_out, driver="ESRI Shapefile", overwrite_layer=TRUE)


		#=====================================================================
		# Apply subsetting and if necessary also orthorectification
		#=====================================================================
		path_dest	<- paste(root, "1_atcor_6s", input, sep="/")

		tmp <- ms_files[strip_id]
		tmp <- unlist(strsplit(tmp,".",fixed=TRUE))
		filename <- tmp[-length(tmp)]
		tmp <- unlist(strsplit(filename,"_"))
		tmp <- tmp[length(tmp)-2]
		name_out <- paste("footprint_","ms_",tmp,"_st_area_",oid,".shp",sep="")
		shape_filename_ms <- paste(path_out,name_out,sep="/")

		name_out <- paste("footprint_","pan_",tmp,"_st_area_",oid,".shp",sep="")
		shape_filename_pan <- paste(path_out_pan,name_out,sep="/")

		# Extract the projection info
		dem_name <- st_area_info_db$dem_name
		dem_file <- paste(stars_dir,"derived","DEM/DEM_SRTM",dem_name,sep="/")

		til_file <- list.files(path_dest, pattern=paste(strip_id,".TIL",sep=""),full.names = TRUE,ignore.case=TRUE)

		info <- GDALinfo(til_file, silent=TRUE)
		proj_code <- attr(info,"projection")
		proj_code <- paste("'",proj_code,"'",sep="")

		out_file <- paste(path_out,paste(input,"_P00",strip_id,"_st_area_",oid,".tif",sep=""),sep="/")

		# Decision is based on description from http://www.c-agg.org/cm_vault/files/docs/DigitalGlobe-Base-Product-FAQ.pdf
		# Only ORStandard2A products are orthorectified
		if(product_type=="OrthoRectified3"){
			# do not apply orthorectification. Only subsetting and mosaicing to single tif file (internally tiled).
			command <- paste("gdalwarp","-dstnodata",miss_val,"-co","TILED=YES","-cutline",shape_filename_ms,"-r","near","-t_srs",proj_code,"-crop_to_cutline","-overwrite",til_file,out_file,sep=" ")
			system(command, intern = FALSE,ignore.stdout = FALSE, ignore.stderr = FALSE, wait = TRUE, input = NULL)
		} else if(product_type=="Standard2A"){
			# do not apply orthorectification. Only subsetting and mosaicing to single tif file (internally tiled).
			command <- paste("gdalwarp","-dstnodata",miss_val,"-co","TILED=YES","-cutline",shape_filename_ms,"-r","near","-t_srs",proj_code,"-crop_to_cutline","-overwrite",til_file,out_file,sep=" ")
			system(command, intern = FALSE,ignore.stdout = FALSE, ignore.stderr = FALSE, wait = TRUE, input = NULL)
		} else if(product_type=="ORStandard2A"){
			# apply orthorectification based on RPC and DEM in addition to subsetting and mosaicing to single tif file (internally tiled).
			command <- paste("gdalwarp","-dstnodata",miss_val,"-co","TILED=YES","-cutline",shape_filename_ms,"-rpc","-to",paste("rpc_dem=",dem_file,sep=""),"-et",0.01,"-r","near","-t_srs",proj_code,"-crop_to_cutline","-overwrite",til_file,out_file,sep=" ")
			system(command, intern = FALSE,ignore.stdout = FALSE, ignore.stderr = FALSE, wait = TRUE, input = NULL)
		}

		# Panchromatic data
		subprod_str_pan <- paste(input,strip_str,"PAN",sep="_")

		path_src 	<- paste(stars_dir, "acquired", "DG", input, subprod_str_pan, sep="/")

		# files <- list.files(path_src, pattern="*", full.names = TRUE, ignore.case=TRUE)
		# file.copy(files, path_dest, overwrite = TRUE, recursive = FALSE, copy.mode = TRUE, copy.date = FALSE)

		til_file <- list.files(path_src, pattern=paste(strip_id,".TIL",sep=""),full.names = TRUE,ignore.case=TRUE)

		info <- GDALinfo(til_file, silent=TRUE)
		proj_code <- attr(info,"projection")
		proj_code <- paste("'",proj_code,"'",sep="")

		out_file <- paste(path_out_pan,paste(input,"_pan","_P00",strip_id,"_st_area_",oid,".tif",sep=""),sep="/")

		# Decision is based on description from http://www.c-agg.org/cm_vault/files/docs/DigitalGlobe-Base-Product-FAQ.pdf
		# Only ORStandard2A products are orthorectified
		if(product_type=="OrthoRectified3"){
			# do not apply orthorectification. Only subsetting and mosaicing to single tif file (internally tiled).
			command <- paste("gdalwarp","-dstnodata",miss_val,"-co","TILED=YES","-cutline",shape_filename_pan,"-r","near","-t_srs",proj_code,"-crop_to_cutline","-overwrite",til_file,out_file,sep=" ")
			system(command, intern = FALSE,ignore.stdout = FALSE, ignore.stderr = FALSE, wait = TRUE, input = NULL)
		} else if(product_type=="Standard2A"){
			# do not apply orthorectification. Only subsetting and mosaicing to single tif file (internally tiled).
			command <- paste("gdalwarp","-dstnodata",miss_val,"-co","TILED=YES","-cutline",shape_filename_pan,"-r","near","-t_srs",proj_code,"-crop_to_cutline","-overwrite",til_file,out_file,sep=" ")
			system(command, intern = FALSE,ignore.stdout = FALSE, ignore.stderr = FALSE, wait = TRUE, input = NULL)
		} else if(product_type=="ORStandard2A"){
			# apply orthorectification based on RPC and DEM in addition to subsetting and mosaicing to single tif file (internally tiled).
			command <- paste("gdalwarp","-dstnodata",miss_val,"-co","TILED=YES","-cutline",shape_filename_pan,"-rpc","-to",paste("rpc_dem=",dem_file,sep=""),"-et",0.01,"-r","near","-t_srs",proj_code,"-crop_to_cutline","-overwrite",til_file,out_file,sep=" ")
			system(command, intern = FALSE,ignore.stdout = FALSE, ignore.stderr = FALSE, wait = TRUE, input = NULL)
		}
	}

	# Remove redundant tiles
	path_dest	<- paste(root, "1_atcor_6s", input, sep="/")
	filename <- paste("redundant_tiles_","P00",strip_id,".RData",sep="")
	setwd(path_dest)
	if(file.exists(filename)){
		load(filename)
		ind_rem <- which(file.exists(tile_list))
		file.remove(tile_list[ind_rem])
	}
}

# Remove metadata
files <- list.files(pattern="P00[0-9].TIL", full.names = FALSE, ignore.case=TRUE)
all_files <- files

files <- list.files(pattern="P00[0-9].IMD", full.names = FALSE, ignore.case=TRUE)
all_files <- c(all_files,files)

files <- list.files(pattern="P00[0-9].RPB", full.names = FALSE, ignore.case=TRUE)
all_files <- c(all_files,files)

files <- list.files(pattern="P00[0-9].RData", full.names = FALSE, ignore.case=TRUE)
all_files <- c(all_files,files)

file.remove(all_files)

# The END
