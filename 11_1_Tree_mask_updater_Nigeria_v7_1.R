#======================================================================================
# Variable definitions, data import, preparation
#======================================================================================
#rm(list=ls(all=TRUE))

graphics.off()

require(rgdal)
require(Rcpp)

WriteLog <- TRUE
Show_graph <- FALSE
Server <- FALSE

Show_ms <- TRUE

if(Server){
	Root <- "/home/tolpekin/stars/derived/DG_v8" 
	Path_lib <- "/home/tolpekin/stars/derived/scripts/v8"
} else{
	Root <- "S:/derived/DG_v8"
	Path_lib <- "S:/derived/scripts/v8"
}

setwd(Path_lib)
#source("scale_space_lib_v7.R")
source("tree_measurement_lib_v9.R")
sourceCpp("matcher_v5.cpp")
sourceCpp("blobs_lib_v4.cpp")
sourceCpp("tree_mask_from_blobs_v3.cpp")

sourceCpp("MRF_lib.cpp")


# For Nigeria
master_id <- "053734892380_01"

Path_images <- paste(Root,"4_categ","NG_Kofa",sep="/")
# Get image names from the input dir: filename contains tif but not aux
setwd(Path_images)
Files <- list.files(".",pattern=".tif",ignore.case = TRUE)
aux_files <- grep(".aux",Files,ignore.case = TRUE)
if(length(aux_files)>0) Files <- Files[-aux_files]
ind_master <- grep(master_id, Files,ignore.case = TRUE)
All_files <- Files

# Set input & output directories
Path_in <- paste(Root, "5_tree_measurement",sep="/")

# Two slave images used to update the mask. Run sequentially.
#image_id <- "054112895100_01"
image_id <- "054399067010_01"

Path_out <- paste(Root, "5_tree_measurement",image_id,sep="/")
if(!file.exists(Path_out))dir.create(Path_out,showWarnings=FALSE, recursive=TRUE)

Path_tmp <- paste(Path_out,"temp",sep="/")
if(!file.exists(Path_tmp))dir.create(Path_tmp,showWarnings=FALSE, recursive=TRUE)

# Open logfile
if(WriteLog){
	setwd(Path_out)
	sink(file=paste("update_tree_mask_with_im=",image_id,".txt",sep=""))
	cat(paste(Sys.time(),"starting","\n",sep=" "))
}

#===========================================================================
# Read tree mask. True position on the ground 
# (already corrected for topographic shift)
#===========================================================================
setwd(Path_in)
# The mask from master image. Resides in S:\derived\DG_v8\5_tree_measurement\053734892380_01
#datafile <- paste("corrected_trees_size_shape_top_100_perc_v6.RData",sep="")
# The mask updated with the first slave image
datafile <- paste("Nigeria_updated_","054112895100_01",".RData",sep="")

if(!file.exists(datafile)){
	if(WriteLog) cat(paste(Sys.time(),"Error: blobs file for master image not found","\n", sep=" "))
	next
}

load(file=datafile)
#Blobs_m <- Master_blobs
Blobs_m <- Tree_mask

Path_ms_images <- paste(Root,"6_image_matching_v2",sep="/")
Path_ms <- paste(Path_ms_images,image_id,sep="/")

Path_pan_images <- paste(Root,"6_matched_pan_v2",sep="/")
Path_pan <- paste(Path_pan_images,image_id,sep="/")

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
Mtile <- 400
Ntile <- 400
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
psi <- atan2(tan(theta_sat)*sin(alpha_sat)-(1+tan(theta_sun))*sin(alpha_sun),tan(theta_sat)*cos(alpha_sat)-(1+tan(theta_sun))*cos(alpha_sun))
corr_factor <- sqrt((tan(theta_sat))^2+(1+tan(theta_sun))^2-2*tan(theta_sat)*(1+tan(theta_sun))*cos(alpha_sat-alpha_sun))
#===========================================================================
# Set up band combinations
#===========================================================================
Nb <- ms.imageinfo[["bands"]]

if(Nb==8){
	# Set RGB composition
	nR <- 7
	nG <- 5
	nB <- 3
	# WV2 8 band NDVI
	nir <- 7
	red <- 5

}else{
	if(Nb==4){
		nR <- 4
		nG <- 3
		nB <- 2

		nir <- 4
		red <- 2
	}else{
		nR <- 1
		nG <- 1
		nB <- 1
	}
}

ntx <- ceiling(M0.ms/Mtile)
nty <- ceiling(N0.ms/Ntile)

ix_arr <- 1:ntx
iy_arr <- 1:nty

if(WriteLog) cat(paste(Sys.time()," Tile size ",Mtile," by ",Ntile,"; tile overlap ",tile_over,"\n",sep=""))
if(WriteLog) cat(paste(Sys.time()," Processing tiles ",min(ix_arr),":",max(ix_arr)," by ",min(iy_arr),":",max(iy_arr),"\n",sep=""))

Tree_mask <- data.frame()
Blobs_update <- data.frame()

# process a single tile
#ix <- 13
#iy <- 1

for(ix in ix_arr)
for(iy in iy_arr){
	i1 <- max((ix-1)*Mtile + 1 - tile_over,1)
	i2 <- min(ix*Mtile + tile_over,M0.ms) 
	j1 <- max((iy-1)*Ntile + 1 - tile_over,1)
	j2 <- min(iy*Ntile + tile_over,N0.ms)

	# read image ms subset 
	ijr <- c(i1,i2,j1,j2)

	Path_in <- Path_ms
	MS <- read_subset(ms.imagefn,ijr[1],ijr[2],ijr[3],ijr[4])

	if(Show_graph){
		MSdisp <- MS
		for(k in 1:Nb)MSdisp@data[,k] <- histstretch(MS@data[,k])
	}

	if(is.na(diff(range(MS@data))) || diff(range(MS@data))==0 || median(MS@data[,3])==0){
		if(WriteLog) cat(paste(Sys.time()," intersection of image",im," and tile x=",ix,"_y=",iy," is empty. Skipping","\n",sep=""))
		next
	}

	bb <- bbox(MS)
	xrl <- bb[1,]
	yrl <- bb[2,]

	# read pan image subset 
	ijr_pan <- xy_to_rowcol(cbind(xrl,yrl),pan.imageinfo)
	ijr_pan[ijr_pan<0] <- 0
	ijr_pan[1] <- ijr_pan[1] + 1
	ijr_pan[3] <- ijr_pan[3] + 1
	#ijr_pan <- (ijr[]-1)*4 + 1

	Path_in <- Path_pan
	Pan <- read_subset(pan.imagefn,ijr_pan[1],ijr_pan[2],ijr_pan[3],ijr_pan[4])

	MS_arr[[im]] <- MS
	Pan_arr[[im]] <- Pan

	xy.pan <- coordinates(Pan)

	if(Show_graph){
		Pandisp <- Pan
		k <- 1
		Pandisp@data[,k] <- histstretch(Pan@data[,k])
	}

	# Subset tree mask blobs
	ind <- which((Blobs_m$x>xrl[1])&(Blobs_m$x<xrl[2])&(Blobs_m$y>yrl[1])&(Blobs_m$y<yrl[2]))
	Blobs_tile <- Blobs_m[ind,]

	# Delete blobs with h=0
	ind <- which((Blobs_tile$h>0)&(Blobs_tile$hf>0))
	if(length(ind)>0){
		Blobs_tile <- Blobs_tile[ind,]

		# Project master blobs onto master image geometry
		h  <- Blobs_tile$h
		hf <- Blobs_tile$hf
		n  <- Blobs_tile$n

		x  <- Blobs_tile$x
		y  <- Blobs_tile$y
		R  <- 0.5 * Blobs_tile$d

		h1 <- h*hf
		h2 <- h-h1

		#start_time <- Sys.time()
		topo_shift <- project_pollock_quantitative_matrixC(h1, h2, R, theta_sat, n)
		#Sys.time() - start_time

		x_proj <- x - topo_shift * cos(alpha_sat)
		y_proj <- y - topo_shift * sin(alpha_sat)

		Blobs_tile_proj <- Blobs_tile

		Blobs_tile_proj$x <- x_proj
		Blobs_tile_proj$y <- y_proj
	} else{
		Blobs_tile <- data.frame()
		Blobs_tile_proj <- data.frame()
	}

	if(Show_graph){
		windows(record=TRUE)
		par(mfrow=c(1,2))

		if(Show_ms){
			image(MSdisp,red=nR,green=nG,blue=nB,axes=TRUE)
			title(main=paste("mask_v1_",acq_date," true positions",sep=""))
			display_all_blobs(Blobs_tile,"white")
			
			# Add shadow contour to master image
			if(nrow(Blobs_tile)>0) for(id in 1:nrow(Blobs_tile)){
				#xtrue <- Blobs_tile_proj_m$x[id]
				#ytrue <- Blobs_tile_proj_m$y[id]
				xtrue <- Blobs_tile$x[id]
				ytrue <- Blobs_tile$y[id]
				
				shad_pol <- project_pollockC(h1[id],h2[id],R[id],xtrue,ytrue,theta_sun,alpha_sun,n[id],0.25)
				#shad_pol <- project_pollock(h1[id],h2[id],R_m[id],xtrue,ytrue,theta_sun_m,n[id],alpha_sun_m)
				lines(shad_pol,col="yellow")
			}

			image(MSdisp,red=nR,green=nG,blue=nB,axes=TRUE)
			title(main=paste("tile x=",ix,"iy=",iy," and projected trees",sep=""))

			display_all_blobs(Blobs_tile_proj,"green")

			# Add shadow contour to master image
			if(nrow(Blobs_tile)>0) for(id in 1:nrow(Blobs_tile)){
				#xtrue <- Blobs_tile_proj_m$x[id]
				#ytrue <- Blobs_tile_proj_m$y[id]
				xtrue <- Blobs_tile$x[id]
				ytrue <- Blobs_tile$y[id]
				
				shad_pol <- project_pollockC(h1[id],h2[id],R[id],xtrue,ytrue,theta_sun,alpha_sun,n[id],0.25)
				#shad_pol <- project_pollock(h1[id],h2[id],R_m[id],xtrue,ytrue,theta_sun_m,n[id],alpha_sun_m)
				lines(shad_pol,col="yellow")
			}
		}

		if(!Show_ms){
		#if(TRUE){
			windows()
			par(mfrow=c(1,2))
			image(Pandisp,col=gray((0:255)/255),axes=TRUE)
			title(main=paste("image",im,"=",acq_date," true positions",sep=""))
			display_all_blobs(Blobs_tile,"white")

			# Add shadow contour to master image
			if(nrow(Blobs_tile)>0) for(id in 1:nrow(Blobs_tile)){
				#xtrue <- Blobs_tile_proj_m$x[id]
				#ytrue <- Blobs_tile_proj_m$y[id]
				xtrue <- Blobs_tile$x[id]
				ytrue <- Blobs_tile$y[id]
				
				shad_pol <- project_pollockC(h1[id],h2[id],R[id],xtrue,ytrue,theta_sun,alpha_sun,n[id],0.25)
				#shad_pol <- project_pollock(h1[id],h2[id],R_m[id],xtrue,ytrue,theta_sun_m,n[id],alpha_sun_m)
				lines(shad_pol,col="yellow")
			}

			image(Pandisp,col=gray((0:255)/255),axes=TRUE)
			title(main=paste("tile x=",ix," y=",iy," and projected trees",sep=""))

			display_all_blobs(Blobs_tile_proj,"green")

			# Add shadow contour to master image
			if(nrow(Blobs_tile)>0) for(id in 1:nrow(Blobs_tile)){
				#xtrue <- Blobs_tile_proj_m$x[id]
				#ytrue <- Blobs_tile_proj_m$y[id]
				xtrue <- Blobs_tile$x[id]
				ytrue <- Blobs_tile$y[id]
				
				shad_pol <- project_pollockC(h1[id],h2[id],R[id],xtrue,ytrue,theta_sun,alpha_sun,n[id],0.25)
				#shad_pol <- project_pollock(h1[id],h2[id],R_m[id],xtrue,ytrue,theta_sun_m,n[id],alpha_sun_m)
				lines(shad_pol,col="yellow")
			}

		}
	}

	# Detect new trees
	#================================================================================
	# Phase 1: detect larger trees
	#================================================================================
	# Run tree detection in the MS image
	# define range of scale values
	# t = sigma^2, in pixels
	Dmin <- 0.0
	#Dsmall <- 5.0
	Dmax <- 40.0
	
	Darr <- seq(from = Dmin, to = Dmax, by = 0.75*mean(ps.ms))
	#Darr1 <- seq(from = Dmin, to = Dsmall-0.125*mean(ps.ms), by = 0.125*mean(ps.ms))
	#Darr2 <- seq(from = Dsmall, to = Dmax, by = 0.25*mean(ps.ms))
	#Darr <- c(Darr1,Darr2)
	
	#tmin <- 0.5*(Dmin/sum(ps.ms))^2
	#tmax <- 0.5*(Dmax/sum(ps.ms))^2
	#tarr <- seq(from = tmin, to = tmax, length.out=100)

	tarr <- 0.5*(Darr/sum(ps.ms))^2

	Ns <- length(tarr)

	xy.ms <- coordinates(MS)
	xrl <- range(xy.ms[,1])
	yrl <- range(xy.ms[,2])

	y <- data.matrix(MS@data, rownames.force = NA)

	ndvi <- (y[,nir]-y[,red])/(y[,nir]+y[,red])
	if(Show_graph){
		MSdisp <- MS
		for(k in 1:Nb)MSdisp@data[,k] <- histstretch(MS@data[,k])
		MSdisp$ndvi <- histstretch(ndvi)
		
		if(FALSE){
			windows()
			par(mfrow=c(1,2))
			image(MSdisp,red=nR,green=nG,blue=nB,axes=TRUE)
			title(main=paste("tile x=",ix,"y=",iy,sep=""))
			
			#windows()
			image(MSdisp,attr="ndvi",col=gray((0:255)/255),axes=TRUE)
			title(main=paste("NDVI"," tile x=",ix,"y=",iy,sep=""))
		}
	}

	MS$ndvi <- ndvi

	#load(file=paste("allblobs_image_",im,"_tx_",ix,"_ty_",iy,".RData",sep=""))
	#if((ix==1)&(iy==1)) Blobs_all <- Blobs else Blobs_all <- rbind(Blobs_all,Blobs)

	M <- MS@grid@cells.dim[1]
	N <- MS@grid@cells.dim[2]

	P <- MS$ndvi

	Debug <- FALSE

	xTL <- xy.ms[1,1]
	yTL <- xy.ms[1,2]

	#Thresholds for magnitude and ndvi
	magn_thresh <- 1.0e-04
	ndvi_thresh <- 0.10

	# Interest point detection
	start_time <- Sys.time()
	Blobs <- detect_blobs_v3(P, M, N, ps.ms, xTL, yTL, tarr, magn_thresh)
	end_time <- Sys.time()
	end_time - start_time

	# Evaluate ndvi of blobs
	ndvi_blobs <- measure_blob_ndvi(Blobs, M, N, ndvi, ps.ms, xTL,yTL)

	Blobs$ndvi <- ndvi_blobs

	ind <- which(Blobs$ndvi >= ndvi_thresh)

	if(length(ind)>0){
		Blobs <- Blobs[ind,]
	} else Blobs <- data.frame(array(0,c(0,5)))

	if(nrow(Blobs)>1000){
		ind <- order(Blobs$magn, decreasing=TRUE)
		ind <- ind[1:1000]
		Blobs <- Blobs[ind,]
	}

	if(nrow(Blobs)>0){
		# Delete objects near image margins
		dx1 <- Blobs$x - xrl[1]
		dx2 <- -Blobs$x + xrl[2]
		
		dy1 <- Blobs$y - yrl[1]
		dy2 <- -Blobs$y + yrl[2]
		
		ind <- which(pmin(dx1,dx2,dy1,dy2) <= 0.5*Blobs$d + mean(ps.ms))
		if(length(ind)>0) Blobs <- Blobs[-ind,]
	}

	if(nrow(Blobs)>0){
		Blobs <- Blobs[Blobs$d>=5.0,]
		if(nrow(Blobs)>0) Blobs <- clean_cocentric_blobs(Blobs)
	}

	if(FALSE) if(Show_graph){
		windows()
		par(mfrow=c(1,2))
		image(MSdisp,red=nR,green=nG,blue=nB,axes=TRUE)
		title(main=paste("tile x=",ix,"y=",iy,sep=""))
		display_all_blobs(Blobs,"white")
		#if(nrow(Blobs)>0)text(Blobs$x,Blobs$y,labels=1:nrow(Blobs),pos=4,col="green")
		#if(nrow(Blobs)>0)text(Blobs$x,Blobs$y,labels=1:nrow(Blobs),pos=4)
		
		#windows()
		#image(MSdisp,attr="ndvi",col=gray((0:255)/255),axes=TRUE)
		image(MSdisp,red=nR,green=nG,blue=nB,axes=TRUE)
		title(main=paste("NDVI"," tile x=",ix,"y=",iy,sep=""))
		display_all_blobs(Blobs_tile_proj,"green")
		display_all_blobs(Blobs,"white")
		#if(nrow(Blobs)>0)text(Blobs$x,Blobs$y,labels=1:nrow(Blobs),pos=4,col="green")
	}

	if(nrow(Blobs)>0) Large_blobs <- clean_overlap_tree_mask(Blobs) else Large_blobs <- Blobs

	if(FALSE) if(Show_graph){
		windows()
		par(mfrow=c(1,2))
		image(MSdisp,red=nR,green=nG,blue=nB,axes=TRUE)
		title(main=paste("Larger_trees"," tile x=",ix,"y=",iy,sep=""))
		display_all_blobs(Large_blobs,"white")
		#if(nrow(Blobs)>0)text(Blobs$x,Blobs$y,labels=1:nrow(Blobs),pos=4,col="green")
		#if(nrow(Blobs)>0)text(Blobs$x,Blobs$y,labels=1:nrow(Blobs),pos=4)
		
		#windows()
		#image(MSdisp,attr="ndvi",col=gray((0:255)/255),axes=TRUE)
		image(MSdisp,red=nR,green=nG,blue=nB,axes=TRUE)
		title(main=paste("Larger_trees"," NDVI"," tile x=",ix,"y=",iy,sep=""))
		display_all_blobs(Blobs_tile_proj,"green")
		display_all_blobs(Large_blobs,"white")
		#if(nrow(Blobs)>0)text(Blobs$x,Blobs$y,labels=1:nrow(Blobs),pos=4,col="green")
	}

	#================================================================================
	# Phase 2: detect smaller trees
	#================================================================================
	if(TRUE){
		# Run tree detection in the MS image
		# define range of scale values
		# t = sigma^2, in pixels
		Dmin <- 0.0
		#Dsmall <- 5.0
		Dmax <- 8.0
		
		Darr <- seq(from = Dmin, to = Dmax, by = 0.1*mean(ps.ms))
		#Darr1 <- seq(from = Dmin, to = Dsmall-0.125*mean(ps.ms), by = 0.125*mean(ps.ms))
		#Darr2 <- seq(from = Dsmall, to = Dmax, by = 0.25*mean(ps.ms))
		#Darr <- c(Darr1,Darr2)
		
		#tmin <- 0.5*(Dmin/sum(ps.ms))^2
		#tmax <- 0.5*(Dmax/sum(ps.ms))^2
		#tarr <- seq(from = tmin, to = tmax, length.out=100)

		tarr <- 0.5*(Darr/sum(ps.ms))^2

		Ns <- length(tarr)

		#Thresholds for magnitude and ndvi
		magn_thresh <- 1.0e-03
		ndvi_thresh <- 0.15

		# Interest point detection
		start_time <- Sys.time()
		Blobs <- detect_blobs_v3(P, M,N,ps.ms, xTL, yTL, tarr, magn_thresh)
		end_time <- Sys.time()
		end_time - start_time

		# Evaluate ndvi of blobs
		ndvi_blobs <- measure_blob_ndvi(Blobs, M, N, ndvi, ps.ms, xTL,yTL)

		Blobs$ndvi <- ndvi_blobs

		ind <- which(Blobs$ndvi >= ndvi_thresh)
		if(length(ind)>0){
			Blobs <- Blobs[ind,]
		} else Blobs <- data.frame(array(0,c(0,5)))

		if(FALSE)if(nrow(Blobs)>0){
			# Delete objects near image margins
			dx1 <- Blobs$x - xrl[1]
			dx2 <- -Blobs$x + xrl[2]
			
			dy1 <- Blobs$y - yrl[1]
			dy2 <- -Blobs$y + yrl[2]
			
			ind <- which(pmin(dx1,dx2,dy1,dy2) <= 0.5*Blobs$d + mean(ps.ms))
			if(length(ind)>0) Blobs <- Blobs[-ind,]
		}

		if(nrow(Blobs)>0){
			Blobs <- clean_cocentric_blobs(Blobs)
			Small_blobs <- blobs_contained(Blobs, Large_blobs)
		}else Small_blobs <- Blobs

		#Blobs <- Small_blobs
	}
	Blobs <- rbind(Large_blobs, Small_blobs)
	#Blobs <- Large_blobs

	if(nrow(Blobs)>1000){
		ind <- order(Blobs$magn, decreasing=TRUE)
		ind <- ind[1:1000]
		Blobs <- Blobs[ind,]
	}

	# Check which of the newly detected trees are not yet present in the mask
	# Identify non-redundant trees
	x1 <- Blobs_tile_proj$x
	y1 <- Blobs_tile_proj$y
	R1 <- 0.5*Blobs_tile_proj$d

	x2 <- Blobs$x
	y2 <- Blobs$y
	R2 <- 0.5*Blobs$d

	ind_red <- array(1,nrow(Blobs))

	if(nrow(Blobs)>0) for(k in 1:nrow(Blobs)){
		d_arr <- sqrt(((x1-x2[k])^2) + ((y1-y2[k])^2))
		#if(all(d_arr > 0.5*(R1+R2[k]))) ind_red[k] <- 0
		if(all(d_arr > pmax(R1,R2[k]))) ind_red[k] <- 0
	}

	if(Show_graph){
		windows()
		par(mfrow=c(1,2))
		image(MSdisp,red=nR,green=nG,blue=nB,axes=TRUE)
		title(main=paste("new candidates tile x=",ix,"y=",iy,sep=""))
		display_all_blobs(Blobs,"white")
		if(nrow(Blobs)>0)text(Blobs$x,Blobs$y,labels=1:nrow(Blobs),pos=4,col="white")
		#if(nrow(Blobs)>0)text(Blobs$x,Blobs$y,labels=1:nrow(Blobs),pos=4)
		
		#windows()
		#image(MSdisp,attr="ndvi",col=gray((0:255)/255),axes=TRUE)
		image(MSdisp,red=nR,green=nG,blue=nB,axes=TRUE)
		title(main=paste("redundancy check",sep=""))
		display_all_blobs(Blobs_tile_proj,"white")
		display_all_blobs(Blobs[ind_red==0,],"green")
		display_all_blobs(Blobs[ind_red==1,],"blue")
		#if(nrow(Blobs)>0)text(Blobs$x,Blobs$y,labels=1:nrow(Blobs),pos=4,col="green")
	}

	if(WriteLog)cat(paste(Sys.time()," Tile ix=",ix," iy=",iy," candidates found: ",nrow(Blobs),"\n",sep=""))
	# Measure tree size for the non-redundant trees
	Blobs <- Blobs[ind_red==0,]
	if(WriteLog)cat(paste(Sys.time()," non-redundant candidates: ",nrow(Blobs),"\n",sep=""))

	if(nrow(Blobs)>0){
		Blobs$h <- NA
		Blobs$hf <- NA
		Blobs$n <- NA

		#if(Show_graph) image(Pan,col=gray((0:255)/255),axes=TRUE)

		xyTL <- c(min(xy.pan[,1]),max(xy.pan[,2]))
		M <- Pan@grid@cells.dim[1]
		N <- Pan@grid@cells.dim[2]

		#if(Show_graph) display_all_blobs(Blobs[id,],"white")

		#Show_graph <- FALSE
		if(Show_graph) windows()
		#id <- 11

		for(id in 1:nrow(Blobs)){
			#if(Show_graph) windows()
			A <- Pan
			pan <- A$band1
			xy <- xy.pan
			ps <- ps.pan

			xy0 <- c(Blobs$x[id],Blobs$y[id])	# Centroid of observed blob
			R <- as.numeric(0.5*Blobs$d[id])	# and its radius

			i <- 1 + round((xy0[1]-xyTL[1])/ps.pan[1])
			j <- 1 + round((xyTL[2]-xy0[2])/ps.pan[2])
			pn <- i + (j-1)*M
			#if(Show_graph) points(xy.pan[pn,1],xy.pan[pn,2],col="blue",pch=16)

			# Identify relevant subset of the pan image
			S_max <- R * corr_factor * 2.0

			i1 <- i - round(2*R/ps.pan[1])
			i2 <- i + round(2*R/ps.pan[1])

			j1 <- j - round(2*R/ps.pan[2])
			j2 <- j + round(2*R/ps.pan[2])

			if(cos(psi)>0){
				i2 <- i2 + round(S_max*cos(psi)/ps.pan[1])
			}else i1 <- i1 + round(S_max*cos(psi)/ps.pan[1])

			if(sin(psi)>0){
				j1 <- j1 - round(S_max*sin(psi)/ps.pan[2])
			}else j2 <- j2 - round(S_max*sin(psi)/ps.pan[2])

			j1 <- max(c(j1,1))
			j2 <- min(c(j2),N)

			i1 <- max(c(i1,1))
			i2 <- min(c(i2,M))

			Pan_sub <- Pan[j1:j2,i1:i2]

			# Analyse shadow region		
			# Randomize segmentation
			pan_sub <- Pan_sub@data$band1
			dsub <- Pan_sub@grid@cells.dim
			Msub <- dsub[1]
			Nsub <- dsub[2]

			xy <- coordinates(Pan_sub)

			if(Show_graph){
				#windows()
				image(Pan_sub,col=gray((0:255)/255),axes=TRUE)
				title(main=id)
				display_all_blobs(Blobs[id,],"white")
			}

			xyc <- xy
			xyc[,1] <- xyc[,1] - xy0[1]
			xyc[,2] <- xyc[,2] - xy0[2]
			rc <- sqrt(rowSums(xyc^2))
			phic <- atan2(xyc[,2],xyc[,1])

			ind <- which((abs(phic-psi)<=pi/4)&(rc<=2*R)&(rc>0.75*R))
			ind_seed <- ind[which.min(pan_sub[ind])]

			if(length(ind_seed)==0) next

			min_val <- pan_sub[ind_seed]+1
			max_val <- max(pan_sub,na.rm=TRUE)-1
			if(min_val>=max_val-5){
				next
			}
			shad_thr_arr <- seq(min_val,max_val,5)
			#shad_thr_arr <- seq(pan_sub[ind_seed]+1,270,5)

			Nobs <- length(shad_thr_arr)

			cover_fun <- array(0,Msub*Nsub)
			area_arr <- array(0,0)

			for(i in 1:Nobs){
				#f <- pan
				f <- pan_sub
				f[] <- 0
				f[pan_sub<shad_thr_arr[i]] <- 1

				#if(f[ind_seed]==0){
				#	# no shadow found
				#	next
				#}

				shad <- Grow_region_seedC(f,Msub,Nsub,ind_seed)

				area_arr <- c(area_arr,length(shad))

				cover_fun[shad] <- cover_fun[shad]+1

				#if(Show_graph) points(xy[shad,,drop=F],col="green",cex=0.2,pch=16)

				# Too large area, neglect
				if(length(shad)*prod(ps.pan)> 0.5*pi*(R^2)){
					if(i>5){
						area_rate <- diff(area_arr)
						if(area_rate[i-1] > 5*area_rate[i-2]) break
					}

					if(length(shad)*prod(ps.pan)> 3.75*pi*(R^2)/corr_factor) break
				}

			}

			#windows()
			#plot(shad_thr_arr[1:length(area_arr)],area_arr)
			#windows()
			#plot(shad_thr_arr[1:length(area_rate)],area_rate)

			Nobs_act <- i-3
			if(Nobs_act<2) next

			cover_fun <- cover_fun/Nobs_act
			#summary(cover_fun)

			# p-level set, median
			ind_median <- which(cover_fun>=0.5)
			#if(Show_graph) points(xy[ind_median,],col="green",cex=0.2,pch=16)
			#points(xy[,],col="green",cex=cover_fun,pch=16)

			shad <- ind_median

			if(length(shad)<2) next

			# delete pixels that are inside the apparent crown
			xys <- xy[shad,,drop=FALSE]

			xys[,1] <- xys[,1] - xy0[1]
			xys[,2] <- xys[,2] - xy0[2]

			r <- sqrt(rowSums(xys^2))
			ind <- which(r>=R)
			shad <- shad[ind]
			if(Show_graph) points(xy[shad,],col="red",cex=0.2,pch=16)

			if(length(shad)==0) next

			xy_shad_pix <- xy[shad,,drop=FALSE]

			# Initial estimate of h
			h0 <- 30.0

			hf <- 0.75
			n <- 2.0
			#h <- 0.5
			h <- 2*R/h0

			h_min <- 0
			h_max <- 2.0
			hf_min <- 0
			hf_max <- 0.9
			n_min <- 2.0
			n_max <- 2.0

			dh <- 0.1
			dhf <- 0.1
			dn <- 0.0

			obj_fun <- eval_shad_v6

			if(Debug){
				params <- c(h,hf,n)
				dpar <- c(dh,dhf,dn)
				lower <- c(h_min,hf_min,n_min)
				upper <-c(h_max,hf_max,n_max)
			}

			grid_fit <- grid_optim_v4(c(h,hf,n),obj_fun,dpar=c(dh,dhf,dn),lower=c(h_min,hf_min,n_min),upper=c(h_max,hf_max,n_max))

			not_optim <- TRUE
			iter <- 0
			eps <- 0.001
			err_arr <- array(0,0)

			err_arr <- c(err_arr,grid_fit[4])
			
			h  <- grid_fit[1]
			h <- max(c(h,0))
			h <- min(c(h,2))
			
			#h_min <- max(c(h-0.2,0))
			#h_max <- min(c(h+0.2,2))
			h_min <- grid_fit[5]
			h_max <- grid_fit[6]

			dh <- max(c((h_max-h_min)/10,0.01))
			
			hf <- grid_fit[2]
			hf <- max(c(hf,0))
			hf <- min(c(hf,0.9))

			dhf <- max(c((hf_max-hf_min)/5,0.01))

			#hf_min <- max(c(hf-0.1,0))
			#hf_max <- min(c(hf+0.1,0.9))
			hf_min <- grid_fit[7]
			hf_max <- grid_fit[8]
			
			n  <- grid_fit[3]
			n <- max(c(n,1))
			n <- min(c(n,3))

			#n_min <- max(c(n-0.2,1))
			#n_max <- min(c(n+0.2,3))
			n_min <- grid_fit[9]
			n_max <- grid_fit[10]

			#dn <- max(c((n_max-n_min)/10,0.1))
			dn <- 0

			if(Debug) if(Show_graph) draw_shad_fit(paste("ix=",ix," iy=",iy," id=",id,sep=""))

			while(not_optim){
				if(Debug){
					params <- c(h,hf,n)
					dpar <- c(dh,dhf,dn)
					lower <- c(h_min,hf_min,n_min)
					upper <-c(h_max,hf_max,n_max)
				}
			
				grid_fit2 <- grid_optim_v4(c(h,hf,n),obj_fun,dpar=c(dh,dhf,dn),lower=c(h_min,hf_min,n_min),upper=c(h_max,hf_max,n_max))

				iter <- iter+1
				
				err_arr <- c(err_arr,grid_fit[4])
			
				h  <- grid_fit[1]
				h <- max(h,h_min)
				h <- min(h,h_max)
				
				hf <- grid_fit[2]
				hf <- max(hf,hf_min)
				hf <- min(hf,hf_max)
				
				n  <- grid_fit[3]
				n <- max(n,n_min)
				n <- min(n,n_max)

				if(Debug) draw_shad_fit(paste("ix=",ix," iy=",iy," id=",id,sep=""))
				
				#converg <- sqrt(mean((grid_fit2-grid_fit)^2))
				converg <- abs(grid_fit2[4]-grid_fit[4])
				grid_fit <- grid_fit2
				if((converg<=eps)|(iter>10)) not_optim <- FALSE
			}

			#if(Show_graph){
			#	#windows()
			#	image(Pan_sub,col=gray((0:255)/255),axes=TRUE)
			#	display_all_blobs(Blobs[id,],"green")
			#}

			if(Show_graph){
				draw_shad_fit(paste("ix=",ix," iy=",iy," id=",id,sep=""))
				#points(xy[ind_median,],col="green",cex=0.2,pch=16)
				#points(xy[shad,],col="red",cex=0.1,pch=16)
			}

			# Add the min margin
			#draw_shadow(xy0,R,h*h0,hf_min,n,theta_sat,alpha_sat,theta_sun,alpha_sun)
			#draw_shadow(xy0,R,h*h0,hf_max,n,theta_sat,alpha_sat,theta_sun,alpha_sun)

			#draw_shadow(xy0,R,h*h0,hf,3,theta_sat,alpha_sat,theta_sun,alpha_sun)
			#draw_shadow(xy0,R,h*h0,0,3,theta_sat,alpha_sat,theta_sun,alpha_sun)

			Blobs$h[id] <- h*h0
			Blobs$hf[id] <- hf
			Blobs$n[id] <- n
		}

		#ind <- which(!is.na(Blobs$h))
		#ind <- which(!is.na(Blobs$h) & Blobs$h>0 & Blobs$hf>0)
		ind <- which(!is.na(Blobs$h) & Blobs$h>0)

		if(length(ind)>0) Blobs <- Blobs[ind,] else Blobs <- data.frame()
	}

	if(WriteLog)cat(paste(Sys.time()," add ",nrow(Blobs)," candidates to the mask","\n",sep=""))

	#Show_graph <- TRUE
	
	if(FALSE)if(Show_graph){
		windows()
		par(mfrow=c(1,2))
		image(MSdisp,red=nR,green=nG,blue=nB,axes=TRUE)
		title(main=paste("tile x=",ix,"y=",iy,sep=""))
		display_all_blobs(Blobs,"white")
		#if(nrow(Blobs)>0)text(Blobs$x,Blobs$y,labels=1:nrow(Blobs),pos=4,col="green")
		#if(nrow(Blobs)>0)text(Blobs$x,Blobs$y,labels=1:nrow(Blobs),pos=4)
		
		#windows()
		#image(MSdisp,attr="ndvi",col=gray((0:255)/255),axes=TRUE)
		image(MSdisp,red=nR,green=nG,blue=nB,axes=TRUE)
		title(main=paste("NDVI"," tile x=",ix,"y=",iy,sep=""))
		display_all_blobs(Blobs_tile_proj,"white")
		display_all_blobs(Blobs,"green")
		#if(nrow(Blobs)>0)text(Blobs$x,Blobs$y,labels=1:nrow(Blobs),pos=4,col="green")
	}

	# Shift new trees to the true position
	Blobs_add <- Blobs
	
	if(nrow(Blobs)>0){
		
		h  <- Blobs$h
		hf <- Blobs$hf
		n  <- Blobs$n

		x  <- Blobs$x
		y  <- Blobs$y
		R  <- 0.5 * Blobs$d

		h1 <- h*hf
		h2 <- h-h1

		#start_time <- Sys.time()
		topo_shift <- project_pollock_quantitative_matrixC(h1, h2, R, theta_sat, n)
		#Sys.time() - start_time

		x <- x + topo_shift * cos(alpha_sat)
		y <- y + topo_shift * sin(alpha_sat)

		Blobs_add$x <- x
		Blobs_add$y <- y
	}	

	# Add new trees to the mask
	#Blobs_m <- rbind(Blobs_m,Blobs_add)
	Blobs_update <- rbind(Blobs_update,Blobs_add)
	
	Blobs_tile <- rbind(Blobs_tile,Blobs_add)
	
	Tree_mask <- rbind(Tree_mask,Blobs_tile)
	setwd(Path_tmp)
	save(Blobs_tile, file=paste("updated_trees_tile_","ix=",ix,"iy=",iy,".RData",sep=""))

	if(WriteLog)cat(paste(Sys.time()," tile ix=",ix," iy=",iy," is updated and contains ",nrow(Blobs_tile)," trees","\n",sep=""))

	# Display updated tree mask
	# Project master blobs onto master image geometry
	Blobs_tile_proj <- Blobs_tile
	
	if(nrow(Blobs_tile)>0){
		h  <- Blobs_tile$h
		hf <- Blobs_tile$hf
		n  <- Blobs_tile$n

		x  <- Blobs_tile$x
		y  <- Blobs_tile$y
		R  <- 0.5 * Blobs_tile$d

		h1 <- h*hf
		h2 <- h-h1

		#start_time <- Sys.time()
		topo_shift <- project_pollock_quantitative_matrixC(h1, h2, R, theta_sat, n)
		#Sys.time() - start_time

		x_proj <- x - topo_shift * cos(alpha_sat)
		y_proj <- y - topo_shift * sin(alpha_sat)

		Blobs_tile_proj$x <- x_proj
		Blobs_tile_proj$y <- y_proj
	}

	#windows(record=TRUE)
	#par(mfrow=c(1,2))
	if(Show_graph){
		if(Show_ms){
			image(MSdisp,red=nR,green=nG,blue=nB,axes=TRUE)
			title(main=paste("updated mask"," image",im,"=",acq_date," true positions",sep=""))
			display_all_blobs(Blobs_tile,"white")
			
			# Add shadow contour to master image
			if(nrow(Blobs_tile)>0) for(id in 1:nrow(Blobs_tile)){
				#xtrue <- Blobs_tile_proj_m$x[id]
				#ytrue <- Blobs_tile_proj_m$y[id]
				xtrue <- Blobs_tile$x[id]
				ytrue <- Blobs_tile$y[id]
				
				shad_pol <- project_pollockC(h1[id],h2[id],R[id],xtrue,ytrue,theta_sun,alpha_sun,n[id],0.25)
				#shad_pol <- project_pollock(h1[id],h2[id],R_m[id],xtrue,ytrue,theta_sun_m,n[id],alpha_sun_m)
				lines(shad_pol,col="yellow")
			}

			image(MSdisp,red=nR,green=nG,blue=nB,axes=TRUE)
			title(main=paste("tile x=",ix,"iy=",iy," and projected trees",sep=""))

			display_all_blobs(Blobs_tile_proj,"green")

			# Add shadow contour to master image
			if(nrow(Blobs_tile)>0) for(id in 1:nrow(Blobs_tile)){
				#xtrue <- Blobs_tile_proj_m$x[id]
				#ytrue <- Blobs_tile_proj_m$y[id]
				xtrue <- Blobs_tile$x[id]
				ytrue <- Blobs_tile$y[id]
				
				shad_pol <- project_pollockC(h1[id],h2[id],R[id],xtrue,ytrue,theta_sun,alpha_sun,n[id],0.25)
				#shad_pol <- project_pollock(h1[id],h2[id],R_m[id],xtrue,ytrue,theta_sun_m,n[id],alpha_sun_m)
				lines(shad_pol,col="yellow")
			}
		}

		if(!Show_ms){
		#if(TRUE){
			windows()
			par(mfrow=c(1,2))
			image(Pandisp,col=gray((0:255)/255),axes=TRUE)
			title(main=paste("image",im,"=",acq_date," true positions",sep=""))
			display_all_blobs(Blobs_tile,"white")

			# Add shadow contour to master image
			if(nrow(Blobs_tile)>0) for(id in 1:nrow(Blobs_tile)){
				#xtrue <- Blobs_tile_proj_m$x[id]
				#ytrue <- Blobs_tile_proj_m$y[id]
				xtrue <- Blobs_tile$x[id]
				ytrue <- Blobs_tile$y[id]
				
				shad_pol <- project_pollockC(h1[id],h2[id],R[id],xtrue,ytrue,theta_sun,alpha_sun,n[id],0.25)
				#shad_pol <- project_pollock(h1[id],h2[id],R_m[id],xtrue,ytrue,theta_sun_m,n[id],alpha_sun_m)
				lines(shad_pol,col="yellow")
			}

			image(Pandisp,col=gray((0:255)/255),axes=TRUE)
			title(main=paste("tile x=",ix," y=",iy," and projected trees",sep=""))

			display_all_blobs(Blobs_tile_proj,"green")

			# Add shadow contour to master image
			if(nrow(Blobs_tile)>0) for(id in 1:nrow(Blobs_tile)){
				#xtrue <- Blobs_tile_proj_m$x[id]
				#ytrue <- Blobs_tile_proj_m$y[id]
				xtrue <- Blobs_tile$x[id]
				ytrue <- Blobs_tile$y[id]
				
				shad_pol <- project_pollockC(h1[id],h2[id],R[id],xtrue,ytrue,theta_sun,alpha_sun,n[id],0.25)
				#shad_pol <- project_pollock(h1[id],h2[id],R_m[id],xtrue,ytrue,theta_sun_m,n[id],alpha_sun_m)
				lines(shad_pol,col="yellow")
			}

		}
	}

}

setwd(Path_out)
save(Tree_mask, file=paste("updated_trees",".RData",sep=""))

# Close logfile
if(WriteLog){
	cat(paste(Sys.time(),"Updated tree mask contains",nrow(Tree_mask),"trees","\n",sep=" "))
	cat(paste(Sys.time(),"Process ended","\n",sep=" "))
	sink()
}

#==================================================================================
# The End
#==================================================================================
