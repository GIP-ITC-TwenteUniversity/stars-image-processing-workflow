#======================================================================================
# Variable definitions, data import, preparation
#======================================================================================
rm(list=ls(all=TRUE))
graphics.off()

require(rgdal)
require(Rcpp)

args <- commandArgs(trailingOnly = TRUE)
oid <- args[1]

# oid <- Sys.getenv("oid")

WriteLog <- TRUE
Show_graph <- FALSE
Debug <- FALSE

root <- Sys.getenv("baseDG")
path_lib <- Sys.getenv("basescript")

setwd(root)
masters_df <- read.csv("masters.csv", header = TRUE, sep = ",", quote = "\"", dec = ".", fill = TRUE, comment.char = "", stringsAsFactors =FALSE)

ind <- which(masters_df$oid==oid)
master_id <- masters_df$master_id[ind]
strip_str <- masters_df$master_strip[ind]
strip_id <- as.numeric(substr(strip_str,2,4))

strip_str <- paste("P00", strip_id, sep="")
subprod_str <- paste(master_id,strip_str,"MUL",sep="_")

# Set input & output directories
Path_in_ms <- paste(root,"2_orthorectification_gdal",master_id,sep="/")
Path_blobs <- paste(root,"3_blob_detection",master_id,sep="/")
Path_in_pan <- paste(root,"2_orthorectification_gdal_pan",master_id,sep="/")
Path_out <- paste(root, "4_tree_measurement",oid,sep="/")

if(!file.exists(Path_out))dir.create(Path_out,showWarnings=FALSE, recursive=TRUE)

# Open logfile
if(WriteLog){
	setwd(Path_out)
	sink(file=paste("tree_height_from_shadow_image_",master_id,"_",strip_str,".txt",sep=""))
	cat(paste(Sys.time(),"starting","\n",sep=" "))
}

if(is.na(master_id) | identical(master_id, character(0))){
	cat(paste("No master image identified for study area ",oid,"\n",sep=""))
	sink()
	q(save="no")
}

setwd(path_lib)
source("tree_measurement_lib_v10.R")
sourceCpp("MRF_lib.cpp")
sourceCpp("blobs_lib_v4.cpp")
sourceCpp("matcher_v5.cpp")

# Set RGB composition
nR <- 7
nG <- 5
nB <- 3

#Thresholds for magnitude and ndvi
magn_thresh <- 3.0e-03
ndvi_thresh <- 0.3
shad_thr <- 300

dpsi <- 45*pi/180

# Get image names from the input dir: filename contains tif but not aux
Path_in <- Path_in_ms

setwd(Path_in)
Files <- list.files(".",pattern=paste(strip_str,"_st_area_",oid,".tif",sep=""), ignore.case=TRUE)

aux_files <- grep(".aux",Files,ignore.case = TRUE)
if(length(aux_files)>0) Files <- Files[-aux_files]

# Choose image
if(length(Files)!=1){
	cat("Error: file not found or not unique (multispectral master)","\n")
}
ms.imagefn <- Files[1]

setwd(Path_in_pan)
Files <- list.files(".",pattern=paste(strip_str,"_st_area_",oid,".tif",sep=""), ignore.case=TRUE)

aux_files <- grep(".aux",Files,ignore.case = TRUE)
if(length(aux_files)>0) Files <- Files[-aux_files]

if(length(Files)!=1){
	cat("Error: file not found or not unique (panchromatic master)","\n")
}
pan.imagefn <- Files[1]

#===========================================================================
# Read image geometry: sun and satellite position (determined from metadata)
#===========================================================================
path_meta <- paste(root, "0_categ", master_id, subprod_str, sep = "/")

setwd(path_meta)
metadata <- read.csv(paste("metadata_", master_id, ".csv", sep = ""), header=TRUE, stringsAsFactors=FALSE)

sun_az <- metadata$Sun_azimuth * pi/180
sat_az <- metadata$Satellite_azimuth * pi/180
theta_sun <- (90-metadata$Sun_elevation) * pi/180
theta_sat <- (90-metadata$Satellite_elevation) * pi/180

alpha_sun <- pi/2 - sun_az
alpha_sat <- pi/2 - sat_az

while(alpha_sun<0) alpha_sun <- alpha_sun + 2*pi
while(alpha_sat<0) alpha_sat <- alpha_sat + 2*pi

alpha_shad <- alpha_sun + pi
while(alpha_shad>=2*pi) alpha_shad <- alpha_shad - 2*pi
while(alpha_shad<=-2*pi) alpha_shad <- alpha_shad + 2*pi

#psi <- atan2(tan(theta_sat)*sin(alpha_sat)-tan(theta_sun)*sin(alpha_sun),tan(theta_sat)*cos(alpha_sat)-tan(theta_sun)*cos(alpha_sun))
#corr_factor <- cos(psi) / (tan(theta_sat)*cos(alpha_sat)-tan(theta_sun)*cos(alpha_sun))

#psi <- atan2(0.5*tan(theta_sat)*sin(alpha_sat)-tan(theta_sun)*sin(alpha_sun),0.5*tan(theta_sat)*cos(alpha_sat)-tan(theta_sun)*cos(alpha_sun))
#corr_factor <- cos(psi) / (0.5*tan(theta_sat)*cos(alpha_sat)-tan(theta_sun)*cos(alpha_sun))

psi <- atan2(tan(theta_sat)*sin(alpha_sat)-(1+tan(theta_sun))*sin(alpha_sun),tan(theta_sat)*cos(alpha_sat)-(1+tan(theta_sun))*cos(alpha_sun))
corr_factor <- sqrt((tan(theta_sat))^2+(1+tan(theta_sun))^2-2*tan(theta_sat)*(1+tan(theta_sun))*cos(alpha_sat-alpha_sun))

ms.imageinfo <- GDALinfo(paste(Path_in_ms,"/",ms.imagefn,sep=""),silent=TRUE)
pan.imageinfo <- GDALinfo(paste(Path_in_pan,"/",pan.imagefn,sep=""),silent=TRUE)

Nb <- ms.imageinfo[["bands"]]
proj_raster_ms <- attributes(ms.imageinfo)$projection

ps.ms <- c(ms.imageinfo[["res.x"]],ms.imageinfo[["res.x"]])

N0.ms <- ms.imageinfo[["rows"]]
M0.ms <- ms.imageinfo[["columns"]]

proj_raster_pan <- attributes(pan.imageinfo)$projection
ps.pan <- c(pan.imageinfo[["res.x"]],pan.imageinfo[["res.x"]])
N0.pan <- pan.imageinfo[["rows"]]
M0.pan <- pan.imageinfo[["columns"]]

ysign_pan <- attributes(pan.imageinfo)$ysign

# Read all blobs
setwd(Path_blobs)
load(file=paste("all_blobs_image_",master_id,"_",strip_str,"_st_area_",oid,".RData",sep=""))

Blobs_all <- Blobs

# Define tiles
Mtile <- 400
Ntile <- 400
# overlap of tiles; prevents loosing points at the margins
tile_over <- 10

ntx <- ceiling(M0.ms/Mtile)
nty <- ceiling(N0.ms/Ntile)

ix_arr <- 1:ntx
iy_arr <- 1:nty

if(WriteLog) cat(paste(Sys.time()," Tile size ",Mtile," by ",Ntile,"; tile overlap ",tile_over,"\n",sep=""))
if(WriteLog) cat(paste(Sys.time()," Processing tiles ",min(ix_arr),":",max(ix_arr)," by ",min(iy_arr),":",max(iy_arr),"\n",sep=""))

# process a single tile
# ix <- 1
# iy <- 1

All_blobs_3d <- data.frame()

for(ix in ix_arr)
for(iy in iy_arr){

	i1 <- max((ix-1)*Mtile + 1 - tile_over,1)
	i2 <- min(ix*Mtile + tile_over,M0.ms) 
	j1 <- max((iy-1)*Ntile + 1 - tile_over,1)
	j2 <- min(iy*Ntile + tile_over,N0.ms)
	
	# read multspectral image subset 
	ijr <- c(i1,i2,j1,j2)

	Path_in <- Path_in_ms
	MS <- read_subset(ms.imagefn,ijr[1],ijr[2],ijr[3],ijr[4])

	if(Show_graph){
		MSdisp <- MS
		for(k in 1:Nb)MSdisp@data[,k] <- histstretch(MS@data[,k])
		#windows()
		#par(mfrow=c(1,1))
		#image(MSdisp,red=nR,green=nG,blue=nB,axes=TRUE)
		#title(main=paste("Tile x=",ix," y=",iy,sep=""))
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

	if(ijr_pan[2]>=M0.pan) ijr_pan[2] <- M0.pan
	if(ijr_pan[4]>=N0.pan) ijr_pan[4] <- N0.pan

	Path_in <- Path_in_pan
	Pan <- read_subset(pan.imagefn,ijr_pan[1],ijr_pan[2],ijr_pan[3],ijr_pan[4])

	xy.pan <- coordinates(Pan)

	if(Show_graph){
		Pandisp <- Pan
		k<-1
		Pandisp@data[,k] <- histstretch(Pan@data[,k])
	}

	xy.ms <- coordinates(MS)
	y <- data.matrix(MS@data, rownames.force = NA)

	# IKONOS NDVI
	nir <- 4
	red <- 2

	# WV2 8 band NDVI
	#nir <- 8
	#red <- 5

	ndvi <- (y[,nir]-y[,red])/(y[,nir]+y[,red])
	if(Show_graph) MSdisp$ndvi <- histstretch(ndvi)
	MS$ndvi <- ndvi

	# Subset master blobs
	ind <- which((Blobs_all$x>xrl[1])&(Blobs_all$x<xrl[2])&(Blobs_all$y>yrl[1])&(Blobs_all$y<yrl[2]))
	Blobs <- Blobs_all[ind,]
	
	# Thresholding
	ind <- which((Blobs$magn > magn_thresh) & (Blobs$ndvi > ndvi_thresh))
	if(length(ind)==0) next
	Blobs <- Blobs[ind,]
	
	# Sort on magnitude 
	ind <- order(Blobs$magn,decreasing = TRUE)

	# Choose the top 100%
	ntop <- round(nrow(Blobs) * 1.00,0)
	ind <- ind[1:ntop]
	Blobs <- Blobs[ind,]

	rownames(Blobs) <- 1:nrow(Blobs)

	if(nrow(Blobs)==0) next

	#======================================================================
	#Detect shadow by threshodling (varies among the images)
	#======================================================================
	Blobs$h <- NA
	Blobs$hf <- NA
	Blobs$n <- NA

	#if(Show_graph) image(Pan,col=gray((0:255)/255),axes=TRUE)

	xyTL <- c(min(xy.pan[,1]),max(xy.pan[,2]))
	M <- Pan@grid@cells.dim[1]
	N <- Pan@grid@cells.dim[2]

	#if(Show_graph) display_all_blobs(Blobs[id,],"white")

	#id <- 1

	for(id in 1:nrow(Blobs)){
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
		pan_sub[is.na(pan_sub)] <- 0
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

		shad_thr_arr <- seq(pan_sub[ind_seed]+1,max(pan_sub)-1,5)
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

			#if(Show_graph) points(xy[shad,],col="green",cex=0.2,pch=16)

			# Too large area, neglect
			if(length(shad)*prod(ps.pan)> 0.5*pi*(R^2)){
				if(i>3){
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
		if(Nobs_act<=3) next

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


	ind <- which(!is.na(Blobs$h))
	Blobs_3d <- Blobs[ind,]
	All_blobs_3d <- rbind(All_blobs_3d,Blobs_3d)

	if(Show_graph){
		image(Pandisp,col=gray((0:255)/255),axes=TRUE)
		title(main=paste("Tile x=",ix," y=",iy,sep=""))
		display_all_blobs(Blobs,"green")
		if(nrow(Blobs)>0)text(Blobs$x,Blobs$y,labels=1:nrow(Blobs),pos=4,col="white")
	}

	setwd(Path_out)
	save(Blobs_3d,file=paste("trees_size_shape_",master_id,"_",strip_str,"_ix_",ix,"_iy_",iy,".RData",sep=""))
	if(WriteLog)cat(paste(Sys.time(),"tile ix",ix,"iy",iy,"processed with",nrow(Blobs_3d)," out of ",nrow(Blobs),"trees","\n",sep=" "))

}

All_blobs_3d <- unique(All_blobs_3d)

setwd(Path_out)
save(All_blobs_3d,file=paste("all_trees_size_shape_",master_id,"_",strip_str,"_st_area_",oid,".RData",sep=""))

# Shift the blobs to fit the reference master image
# dx_m <- -6.999
# dy_m <- 4.013
dx_m <- 0.0
dy_m <- 0.0

h  <- All_blobs_3d$h
hf <- All_blobs_3d$hf
n  <- All_blobs_3d$n

x_m  <- All_blobs_3d$x
y_m  <- All_blobs_3d$y
R_m  <- 0.5 * All_blobs_3d$d

h1 <- h*hf
h2 <- h-h1

start_time <- Sys.time()
topo_shift_m <- project_pollock_quantitative_matrixC(h1, h2, R_m, theta_sat, n)
Sys.time() - start_time

x_m <- x_m + dx_m + topo_shift_m * cos(alpha_sat)
y_m <- y_m + dy_m + topo_shift_m * sin(alpha_sat)

Master_blobs <- All_blobs_3d
Master_blobs$x <- x_m
Master_blobs$y <- y_m

# Save corrected blobs from master image
setwd(Path_out)
datafile <- paste("corrected_trees_size_shape_",master_id,"_P00",strip_id,".RData",sep="")

save(Master_blobs,file=datafile)

# Close the logfile
if(WriteLog){
	cat(paste(Sys.time(),"Process ended",sep=" "))
	sink()
}

#==================================================================================
# The End
#==================================================================================
