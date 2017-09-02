require(rgdal)
require(rgeos)
require(gdalUtils)

args <- commandArgs(trailingOnly = TRUE)
input <- args[1] 
nr_strips <- args[2]

# input <- Sys.getenv("input")
# nr_strips <- as.numeric(Sys.getenv("nr_strips"))

root <- Sys.getenv("baseDG")

path_in <- paste(root,"1_atcor_6s",input,sep="/")

for(strip_id in 1:nr_strips){

	# Convert the relevant tiles from hdf to tif
	setwd(path_in)
	list_hdf <- list.files(pattern = paste("\\P00",strip_id,".hdf$",sep=""))
	file.rename(list_hdf, gsub("[-]", "_", list_hdf))
	list_hdf <- list.files(pattern = paste("\\P00",strip_id,".hdf$",sep=""))

	strip_str <- paste("P00",strip_id,sep="")
	subprod_str <- paste(input,strip_str,"MUL",sep="_")

	path_categ <- paste(root,"0_categ",input,subprod_str,sep="/")
	tiles <- list.files(path_categ, pattern = ".tif$", ignore.case=TRUE, full.names=TRUE)
	rel_tiles <- basename(tiles)
	tmp_str <- gsub("[-]", "_", rel_tiles)

	for(tile_id in seq_along(rel_tiles)){
		simple_str <- tmp_str[tile_id]
		pat <- paste(unlist(strsplit(simple_str,"_",fixed=TRUE))[2],unlist(strsplit(simple_str,"_",fixed=TRUE))[3],sep="_")
		tilerc <- unlist(strsplit(simple_str,"_",fixed=TRUE))[3]
		tmp <- list.files(path_in,pattern = paste("P00",strip_id,".hdf$",sep=""))
		filename <- grep(tilerc,tmp,fixed=TRUE,value=TRUE)

		# if filename is empty str?
		if(nchar(filename)==0) next
		
		setwd(path_in)
		info <- gdalinfo(filename)
		data <- grep("HDF4_SDS", info, value = TRUE)
		datasub <- gsub("SUBDATASET_\\d+_NAME=|\\s+", "", data)

		# GET METADATA FROM RAW IMAGE
		setwd(path_categ)
		tif_file <- rel_tiles[tile_id]
		imageinfo <- GDALinfo(tif_file, silent=TRUE)

		N0 <- imageinfo[["rows"]]
		M0 <- imageinfo[["columns"]]
		nb <- imageinfo[["bands"]]
		npix <- M0*N0

		xyLL <- c(imageinfo[["ll.x"]],imageinfo[["ll.y"]])
		ps <- c(imageinfo[["res.x"]],imageinfo[["res.x"]])
		proj_ref <- attributes(imageinfo)$projection

		ysign <- attributes(imageinfo)$ysign

		# Entire coordinate range of the image
		xr <- c(xyLL[1],xyLL[1]+M0*ps[1])

		if(ysign<0)
		{
			yr <- c(xyLL[2],xyLL[2]+N0*ps[2])
		}else
		{
			yr <- c(xyLL[2]-N0*ps[2],xyLL[2])
		}
		
		bb_tile <- rbind(xr,yr)

		bounds_ullr <- c(bb_tile[1,1], bb_tile[2,2], bb_tile[1,2], bb_tile[2,1])

		my_data <- data.frame(array(0,c(npix,nb)))
		names(my_data) <- paste("band",1:nb,sep="")

		cellcentre <- bb_tile[,1] + 0.5*ps
		mygrid <- GridTopology(cellcentre, ps, c(M0,N0))
		A <- SpatialGridDataFrame(mygrid, my_data, proj_ref)
		# CONTINUE WITH CONVERSION
		setwd(path_in)

		# EXPORT INDIVIDUAL HDF4 FILES TO TIF
		temp <- paste(path_in, "/temp_", tile_id, "_P00", strip_id, "_", pat, sep = "")
		dir.create(temp, showWarnings = TRUE, recursive = FALSE, mode = "0777")
		for (j in 1:length(datasub)) {
			hdftotif <- gdal_translate(src_dataset = filename, dst_dataset = paste(temp, "/band0", j, "_surface_reflectance", ".tif", sep = ""), of = "GTiff", ot = "UInt16", output_Raster = TRUE, sd_index=j, a_srs = proj_ref, a_ullr = bounds_ullr)
		}

		# STACK INDIVIDUAL BANDS FROM HDF
		for (j in 1:nb) {
			y <- readGDAL(datasub[j], silent=TRUE)
			tmp_band <- y@data$band1	
			#tmp_band[tmp_band<0] <- 2^16-1	
			tmp_band[tmp_band<0] <- NA
			A@data[,j] <- tmp_band
		}

		writeGDAL(dataset = A, fname = tif_file, drivername = "GTiff", type = "UInt16", mvFlag = 2^16-1)

		# Delete content of the temporary directory
		unlink(temp, recursive = TRUE, force = FALSE)
	}
}
