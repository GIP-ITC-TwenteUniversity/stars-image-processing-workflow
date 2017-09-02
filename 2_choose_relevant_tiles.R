require(rgdal)
require(rgeos)
require("RPostgreSQL")
require("wkb")

# EPSG 4326
target_projection <- CRS("+init=epsg:4326")

args <- commandArgs(trailingOnly = TRUE)
input <- args[1]
nr_strips <- args[2]

# input <- Sys.getenv("input")
# nr_strips <- as.numeric(Sys.getenv("nr_strips"))

root <- Sys.getenv("baseDG")
base <- Sys.getenv("base")
stars_dir <- Sys.getenv("stars")

path_origin <- paste(stars_dir,"acquired/DG", input, sep="/")

for(strip_id in 1:nr_strips){
	strip_str <- paste("P00",strip_id,sep="")
	subprod_str <- paste(input,strip_str,"MUL",sep="_")

	path_area <- paste(base,"shapefiles", sep="/")
	path_product <- paste(stars_dir,"acquired","DG",input,subprod_str, sep="/")

	# Determine filenames
	files <- list.files(path_product,pattern="\\.TIF",ignore.case=TRUE)
	# Exclude the aux.xml files
	ind <- grep("aux.xml",files)
	if(length(ind)>0) files <- files[-ind]

	nr_tiles <- length(files)
	to_process <- array(FALSE,nr_tiles)

	# read the study area shapefile
	# Use relevant study areas (determine from Metadata)
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

		query <- postgresqlExecStatement(con,"select s.name as name, ST_AsBinary(st_transform(s.geometry,4326)) as geom, common.cs_utmzone_from_polygon(s.geometry) as utmzonecode from common.study_area as s where oid=$1",oid)

		st_area_info_db <- fetch(query,n=-1)

		dbDisconnect(con)     # closes the connection
		dbUnloadDriver(drv)   # free all the resources on the driver

		study_area_shape <- SpatialPolygonsDataFrame(readWKB(hex2raw(st_area_info_db$geom), 1, target_projection), data.frame(oid=oid))
		# writeOGR(study_area_shape, "/home/tolpekin/Documents", "test", driver="ESRI Shapefile")

		#reproject to match the image projection
		tile_id <- 1
		tile_nm <- files[tile_id]
		setwd(path_product)
		imageinfo <- GDALinfo(tile_nm, silent=TRUE)
		proj_raster <- attr(imageinfo,"projection")
		proj_vector <- study_area_shape@proj4string@projargs

		if(!identical(CRS(proj_raster), CRS(proj_vector))) study_area_shape <- spTransform(study_area_shape, proj_raster)

		# bb_area <- bbox(study_area_shape)

		near_tiles <- array(FALSE,nr_tiles)

		setwd(path_product)

		for(tile_id in 1:nr_tiles){
			# read tile info
			tile_nm <- files[tile_id]
			imageinfo <- GDALinfo(tile_nm, silent=TRUE)
			proj_raster <- attr(imageinfo,"projection")
			
			N0 <- imageinfo[["rows"]]
			M0 <- imageinfo[["columns"]]

			xyLL <- c(imageinfo[["ll.x"]],imageinfo[["ll.y"]])
			ps <- c(imageinfo[["res.x"]],imageinfo[["res.x"]])

			ysign <- attributes(imageinfo)$ysign

			# Entire coordinate range of the image
			xr <- c(xyLL[1],xyLL[1]+M0*ps[1])

			if(ysign<0){
				yr <- c(xyLL[2],xyLL[2]+N0*ps[2])
			}else{
				yr <- c(xyLL[2]-N0*ps[2],xyLL[2])
			}

			# bb_tile <- rbind(xr,yr)

			# xdist <- 0
			# if(max(xr)<min(bb_area[1,])){
				# # Case 1: tile is west of the study area
				# xdist <- min(bb_area[1,]) - max(xr)
			# } else if(min(xr)>max(bb_area[1,])){
				# # Case 2: tile is east of the study area
				# xdist <- min(xr) - max(bb_area[1,])
			# }

			# ydist <- 0
			# if(max(yr)<min(bb_area[2,])){
				# # Case 1: tile is south of the study area
				# ydist <- min(bb_area[2,]) - max(yr)
			# } else if(min(yr)>max(bb_area[2,])){
				# # Case 2: tile is north of the study area
				# ydist <- min(yr) - max(bb_area[2,])
			# }

			# max allowed distance, in meters
			margin <- 200.0

			# if((xdist<margin) & (ydist <margin)) near_tiles[tile_id] <- TRUE

			# Measure the distance from tile corners to the study area polygon
			xy <- data.frame(x=c(xr[1],xr[1],xr[2],xr[2]),y=c(yr[1],yr[2],yr[1],yr[2]))
			pts <- SpatialPoints(xy,proj4string=CRS(proj_raster))
			dst <- gDistance(pts,study_area_shape)
			# determine extent of the tile in m
			til_size <- min(c(M0,N0)*ps)

			if(dst < margin | dst < 0.5*til_size) near_tiles[tile_id] <- TRUE
		}
		
		# Identify the relevant tiles to be atmospherically processed
		ind_relevant <- which(near_tiles)
		to_process[ind_relevant] <- TRUE
	}

	# Place the relevant tiles into the ATCOR input folder 
	ind_relevant <- which(to_process)
	if(length(ind_relevant)>0){
		tiles <- paste(path_product,files[ind_relevant],sep="/")
		path_dest <- paste(root,"0_categ",input,subprod_str,sep="/")
		# check if the files are already in the destination folder
		setwd(path_dest)
		missing_files <- !file.exists(basename(tiles))
		tiles <- tiles[missing_files]
		if(length(tiles)>0)file.copy(tiles, path_dest, overwrite = TRUE, recursive = FALSE, copy.mode = TRUE, copy.date = FALSE)
	}
	
	# The irrelevant tiles are copied to the ATCOR output directory
	ind_irrel <- which(!to_process)
	if(length(ind_irrel)>0){
		tiles <- paste(path_product,files[ind_irrel],sep="/")
		path_dest <- paste(root,"1_atcor_6s",input,sep="/")
		# check if the files are already in the destination folder
		setwd(path_dest)
		missing_files <- !file.exists(basename(tiles))
		tiles <- tiles[missing_files]
		if(length(tiles)>0)	file.copy(tiles, path_dest, overwrite = TRUE, recursive = FALSE, copy.mode = TRUE, copy.date = FALSE)
		tile_list <- basename(tiles)
		save(tile_list,file=paste(path_dest,"/redundant_tiles_","P00",strip_id,".RData",sep=""),ascii=TRUE)
	}
}

# The END
