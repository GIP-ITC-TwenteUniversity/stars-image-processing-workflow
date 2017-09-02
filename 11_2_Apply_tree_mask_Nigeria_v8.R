#======================================================================================
# Variable definitions, data import, preparation
#======================================================================================
rm(list=ls(all=TRUE))

graphics.off()

require(rgdal)
require(Rcpp)

WriteLog <- FALSE
Show_graph <- FALSE
Server <- TRUE

if(Server){
	Root <- "/home/tolpekin/stars/derived/RE_v8" 
	Path_lib <- "/home/tolpekin/stars/derived/scripts/RE_v8"
} else{
	Root <- "S:/derived/RE_v8"
	Path_lib <- "S:/derived/scripts/RE_v8"
}

setwd(Path_lib)
source("tree_measurement_lib_v10.R")
sourceCpp("tree_mask_lib.cpp")

# For Nigeria
#master_id <- "053734892380_01"
master_id <- "054399067010_01"

Path_images <- paste(Root,"4_categ","NG_Kofa",sep="/")
# Get image names from the input dir: filename contains tif but not aux
setwd(Path_images)
Files <- list.files(".",pattern=".tif",ignore.case = TRUE)
aux_files <- grep(".aux",Files,ignore.case = TRUE)
if(length(aux_files)>0) Files <- Files[-aux_files]
ind_master <- grep(master_id, Files,ignore.case = TRUE)
All_files <- Files

N_im <- length(All_files)
image_id_arr <- array("", N_im)
#MS_arr <- vector(mode="list", length = N_im)
#Pan_arr <- vector(mode="list", length = N_im)

for(im in 1:N_im){
	ms.imagefn  <- All_files[im]

	image_id <- paste(unlist(strsplit(ms.imagefn,"_"))[1],"01",sep="_")
	image_id_arr[im] <- image_id
}

# Set input & output directories
Path_in <- paste(Root, "5_tree_measurement",sep="/")
Path_out <- paste(Root, "8_tree_mask",sep="/")

# Open logfile
if(WriteLog){
	setwd(Path_out)
	sink(file="logfile_tree_mask_multitemporal.txt")
	cat(paste(Sys.time(),"starting","\n",sep=" "))
}
#===========================================================================
# Read tree mask. True position on the ground 
# (already corrected for topographic shift)
#===========================================================================
setwd(Path_in)
#datafile <- paste("corrected_trees_size_shape_top_100_perc_v6.RData",sep="")
datafile <- paste("Nigeria_updated_v2_054399067010_01",".RData",sep="")

if(!file.exists(datafile)){
	if(WriteLog) cat(paste(Sys.time(),"Error: blobs file for master image not found","\n", sep=" "))
	next
}

load(file=datafile)
Blobs_m <- Tree_mask

if(Show_graph){
	#windows(16,9,record=TRUE)
	windows(16,9,record=F)
	par(mfrow=c(1,2))
}

#im <- 2

#for(im in 1:N_im){
for(im in 12:N_im){
	image_id <- image_id_arr[im]

	if(WriteLog) cat(paste(Sys.time()," Processing image ",im," out of ",N_im," id=",image_id,"\n",sep=""))

	Path_ms_images <- paste(Root,"6_image_matching_v2",sep="/")
	Path_ms <- paste(Path_ms_images,image_id,sep="/")

	Path_pan_images <- paste(Root,"6_matched_pan_v2",sep="/")
	Path_pan <- paste(Path_pan_images,image_id,sep="/")

	Path_out_mask <- paste(Path_out,image_id,sep="/")
	if(!file.exists(Path_out_mask))dir.create(Path_out_mask,showWarnings=FALSE, recursive=TRUE)

	setwd(Path_ms)
	Files <- list.files(".",pattern=".tif",ignore.case = TRUE)
	aux_files <- grep(".aux",Files,ignore.case = TRUE)
	if(length(aux_files)>0) Files <- Files[-aux_files]
	if(length(Files)<1){
		if(WriteLog)cat(paste(Sys.time(),"Image ",image_id,"Error: cannot identify ms image files","\n",sep=" "))
		next
	}
	ms.imagefn  <- Files[1]

	setwd(Path_pan)
	Files <- list.files(".",pattern=".tif",ignore.case = TRUE)
	aux_files <- grep(".aux",Files,ignore.case = TRUE)
	if(length(aux_files)>0) Files <- Files[-aux_files]
	if(length(Files)<1){
		if(WriteLog)cat(paste(Sys.time(),"Image ",image_id,". Error: cannot identify pan image files","\n",sep=" "))
		next
	}
	pan.imagefn  <- Files[1]

	# Define tiles
	Mtile <- 200
	Ntile <- 200
	# overlap of tiles; prevents loosing points at the margins
	tile_over <- 10

	ms.imageinfo <- GDALinfo(paste(Path_ms,ms.imagefn,sep="/"), silent=TRUE)
	ps.ms <- c(ms.imageinfo[["res.x"]],ms.imageinfo[["res.x"]])

	pan.imageinfo <- GDALinfo(paste(Path_pan,pan.imagefn,sep="/"), silent=TRUE)
	ps.pan <- c(pan.imageinfo[["res.x"]],pan.imageinfo[["res.x"]])

	# Define tiles on the basis of master image
	N0.ms <- ms.imageinfo[["rows"]]
	M0.ms <- ms.imageinfo[["columns"]]

	N0.pan <- pan.imageinfo[["rows"]]
	M0.pan <- pan.imageinfo[["columns"]]

	xTL.ms <- ms.imageinfo[["ll.x"]] + 0.5*ps.ms[1]
	yTL.ms <- ms.imageinfo[["ll.y"]] + N0.ms*ps.ms[2] - 0.5*ps.ms[2]
	yBL.ms <- ms.imageinfo[["ll.y"]] + 0.5*ps.ms[2]
	
	proj.ms <- attributes(ms.imageinfo)$projection

	xTL.pan <- pan.imageinfo[["ll.x"]] + 0.5*ps.pan[1]
	yTL.pan <- pan.imageinfo[["ll.y"]] + N0.pan*ps.pan[2] - 0.5*ps.pan[2]
	yBL.pan <- pan.imageinfo[["ll.y"]] + 0.5*ps.pan[2]
	
	proj.pan <- attributes(pan.imageinfo)$projection
	
	ms_tree_mask <- array(0,N0.ms*M0.ms)
	pan_tree_mask <- array(0,N0.pan*M0.pan)
	ms_shad_mask <- array(0,N0.ms*M0.ms)
	pan_shad_mask <- array(0,N0.pan*M0.pan)

	Rbuffer <- 1.0*mean(ps.ms)
	shad_buffer <- 0.0*mean(ps.ms)
	#===========================================================================
	# Read image geometry: sun and satellite position (determined from metadata)
	#===========================================================================
	#input <- args[1] #input="053613698020_01"
	input <- image_id

	image <- paste(input, "_P001_MUL", ".tif", sep = "")
	imagebase <- unlist(strsplit(image, "[.]"))[1]

	pathmeta <- paste(Root, "0_categ", input, imagebase, sep = "/")
	# GET METADATA
	setwd(pathmeta)
	metadata <- read.table(paste("metadata_", input, ".txt", sep = ""), stringsAsFactors=FALSE)

	sun_sat <- as.numeric(unlist(metadata[4:8,2]))
	
	acq_date <- unlist(metadata[11,2])

	sun_az <- sun_sat[1] * pi/180
	sat_az <- sun_sat[3] * pi/180
	theta_sun <- (90-sun_sat[2]) * pi/180
	theta_sat <- (90-sun_sat[4]) * pi/180

	alpha_sun <- pi/2 - sun_az
	alpha_sat <- pi/2 - sat_az
	
	while(alpha_sun<0) alpha_sun <- alpha_sun + 2*pi
	#===========================================================================
	# Project master blobs onto master image geometry
	#===========================================================================
	if(WriteLog)cat(paste(Sys.time()," projecting the mask onto the image geometry","\n",sep=""))

	h_arr  <- Blobs_m$h
	hf_arr <- Blobs_m$hf
	n_arr  <- Blobs_m$n

	x_arr  <- Blobs_m$x
	y_arr  <- Blobs_m$y
	R_arr  <- 0.5 * Blobs_m$d

	h1_arr <- h_arr*hf_arr
	h2_arr <- h_arr-h1_arr

	#start_time <- Sys.time()
	topo_shift <- project_pollock_quantitative_matrixC(h1_arr, h2_arr, R_arr, theta_sat, n_arr)
	#Sys.time() - start_time

	x_proj <- x_arr - topo_shift * cos(alpha_sat)
	y_proj <- y_arr - topo_shift * sin(alpha_sat)

	Blobs_proj <- Blobs_m

	Blobs_proj$x <- x_proj
	Blobs_proj$y <- y_proj
	#===========================================================================
	# Convert blobs to mask: tree mask
	#===========================================================================
	if(WriteLog)cat(paste(Sys.time()," Obtain tree mask","\n",sep=""))

	buffer <- 0.0

	Blobs_matrix <- array(0,c(nrow(Blobs_proj),3))
	Blobs_matrix[,1] <- Blobs_proj$x
	Blobs_matrix[,2] <- Blobs_proj$y
	Blobs_matrix[,3] <- 0.5*Blobs_proj$d + Rbuffer

	#start_time <- Sys.time()
	ms_tree_mask  <- tree_mask_from_blobs(Blobs_matrix, M0.ms, N0.ms, xTL.ms, yTL.ms, ps.ms, buffer)
	#end_time <- Sys.time()
	#end_time - start_time

	#start_time <- Sys.time()
	pan_tree_mask  <- tree_mask_from_blobs(Blobs_matrix, M0.pan, N0.pan, xTL.pan, yTL.pan, ps.pan, buffer)
	#end_time <- Sys.time()
	#end_time - start_time

	#===========================================================================
	# Shadow mask
	#===========================================================================
	if(WriteLog)cat(paste(Sys.time()," Obtain tree shadow mask","\n",sep=""))

	for(id in 1:nrow(Blobs_proj)){
		R <- 0.5*Blobs_m$d[id] + Rbuffer

		xtrue <- Blobs_m$x[id]
		ytrue <- Blobs_m$y[id]
		h <- Blobs_m$h[id]
		hf <- Blobs_m$hf[id]
		n  <- Blobs_m$n[id]
		
		h1 <- h*hf
		h2 <- h - h1

		shad_pol <- shadow_contour(h1, h2, R, xtrue, ytrue, theta_sun, alpha_sun, n, 36)

		pix_tree <- pol_to_raster(shad_pol,M0.ms,N0.ms,xTL.ms,yTL.ms,ps.ms)
		#ms_shad_mask[pix_tree] <- 1
		ms_tree_mask[pix_tree] <- 1

		pix_tree <- pol_to_raster(shad_pol,M0.pan,N0.pan,xTL.pan,yTL.pan,ps.pan)
		#pan_shad_mask[pix_tree] <- 1
		pan_tree_mask[pix_tree] <- 1
	}	
	#===========================================================================
	# Save raster masks
	#===========================================================================
	if(WriteLog)cat(paste(Sys.time()," Save the masks","\n",sep=""))

	ms_grid <- GridTopology(c(xTL.ms,yBL.ms), ps.ms, c(M0.ms,N0.ms))
	Mask_ms <- SpatialGridDataFrame(ms_grid, data.frame(tree=ms_tree_mask), proj.ms)

	setwd(Path_out_mask)

	Mask_ms$tree <- Mask_ms$tree * 255

	OUT <- Mask_ms
	imagefn.out <- paste("Tree_mask_",image_id,".tif",sep="")

	OUT.tif<-create2GDAL(OUT,drivername="GTiff",type="Byte")
	saveDataset(OUT.tif,imagefn.out)
	GDAL.close(OUT.tif)

	pan_grid <- GridTopology(c(xTL.pan,yBL.pan), ps.pan, c(M0.pan,N0.pan))
	Mask_pan <- SpatialGridDataFrame(pan_grid, data.frame(tree=pan_tree_mask), proj.pan)

	Mask_pan$tree <- Mask_pan$tree * 255

	OUT <- Mask_pan
	imagefn.out <- paste("Tree_mask_pan_",image_id,".tif",sep="")

	OUT.tif<-create2GDAL(OUT,drivername="GTiff",type="Byte")
	saveDataset(OUT.tif,imagefn.out)
	GDAL.close(OUT.tif)

	if(WriteLog)cat(paste(Sys.time()," Image id=",image_id," is processed","\n",sep=""))

}


# Close logfile
if(WriteLog){
	cat(paste(Sys.time(),"Process ended","\n",sep=" "))
	sink()
}

#==================================================================================
# The End
#==================================================================================
