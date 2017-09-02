#================================================================
# Variable definitions, data import, preparation
#================================================================
rm(list=ls(all=TRUE))

require(rgdal)
require(mmand)
require(Rcpp)
require(MASS)
#================================================================
# IMAGE INPUT
#================================================================
args <- commandArgs(trailingOnly = TRUE)
input <- args[1]
nr_strips <- args[2]

# input <- Sys.getenv("input")
# nr_strips <- as.numeric(Sys.getenv("nr_strips"))

WriteLog <- TRUE
Show_graph <- FALSE

root <- Sys.getenv("baseDG")
Path_lib <- Sys.getenv("basescript")

setwd(Path_lib)

sourceCpp("blobs_lib_v4.cpp")
source("scale_space_lib_v7_1.R")
#sourceCpp("tree_mask_from_blobs_v3.cpp")

# Set input & output directories
Path_in <- paste(root, "2_orthorectification_gdal",input, sep="/")
Path_out <- paste(root, "3_blob_detection", input, sep="/")
if(!file.exists(Path_out))dir.create(Path_out,showWarnings=FALSE, recursive=TRUE)

for(strip_id in 1:nr_strips){
	strip_str <- paste("P00", strip_id, sep="")
	subprod_str <- paste(input,strip_str,"MUL",sep="_")

	pathmeta <- paste(root, "0_categ", input, subprod_str, sep = "/")

	# GET METADATA
	setwd(pathmeta)
	metadata <- read.csv(paste("metadata_", input, ".csv", sep = ""), header=TRUE, stringsAsFactors=FALSE)

	nr_fixed <- ncol(metadata) - metadata$Nr_study_areas

	for(i_st_area in 1:metadata$Nr_study_areas){
		oid <- metadata[1,nr_fixed+i_st_area]

		Satellite <- metadata$Satellite
		Band_nr <- metadata$Number_of_bands

		# Set NIR and RED bands according to band number
		if(Band_nr==4){
			nir <- 4
			red <- 3
		} else if(Band_nr==8){
			nir <- 7
			red <- 5
		} else{
			nir <- "Unknown"
			red <- "Unknown"
		}

		# Set RGB composition
		if(Show_graph){
			if(Band_nr==4){
				nR <- 4
				nG <- 3
				nB <- 2
			} else if(Band_nr==8){
				nR <- 7
				nG <- 5
				nB <- 3
			} else if(Band_nr==5){
				nR <- 5
				nG <- 3
				nB <- 2
			} else {
				nR <- "Unknown"
				nG <- "Unknown"
				nB <- "Unknown"
			}
		} 

		# Define tiles
		Mtile <- 400
		Ntile <- 400
		# overlap of tiles; prevents loosing points at the margins
		tile_over <- 10

		# Get image names from the input dir: filename contains tif but not aux
		setwd(Path_in)
		Files <- list.files(".",pattern = paste("P00",strip_id,"_st_area_",oid,".tif", sep = ""), ignore.case = TRUE)

		# Choose image
		if(length(Files)!=1){
			stop("Error: file not found or not unique (multispectral master)","\n")
		}

		ms.imagefn <- Files[1]

		# Open logfile
		if(WriteLog){
			setwd(Path_out)
			sink(file=paste("blob_detection_",ms.imagefn,".txt",sep=""))
			cat(paste(Sys.time(),"starting processing of image ",ms.imagefn,"\n",sep=" "))
		}

		ms.imageinfo <- GDALinfo(paste(Path_in,"/",ms.imagefn,sep=""))
		#pan.imageinfo <- GDALinfo(paste(Path_in,"/",pan.imagefn,sep=""))

		Nb <- ms.imageinfo[["bands"]]
		proj_raster_ms <- attributes(ms.imageinfo)$projection
		ps.ms <- c(ms.imageinfo[["res.x"]],ms.imageinfo[["res.x"]])
		N0.ms <- ms.imageinfo[["rows"]]
		M0.ms <- ms.imageinfo[["columns"]]

		ntx <- ceiling(M0.ms/Mtile)
		nty <- ceiling(N0.ms/Ntile)

		ix_arr <- 1:ntx
		iy_arr <- 1:nty

		if(WriteLog) cat(paste(Sys.time()," Tile size ",Mtile," by ",Ntile,"; tile overlap ",tile_over,"\n",sep=""))
		if(WriteLog) cat(paste(Sys.time()," Processing tiles ",min(ix_arr),":",max(ix_arr)," by ",min(iy_arr),":",max(iy_arr),"\n",sep=""))

		setwd(Path_out)

		#ix <- 2
		#iy <- 4

		start_time0 <- Sys.time()

		All_blobs <- data.frame()

		for(ix in ix_arr)
		for(iy in iy_arr){

			i1 <- max((ix-1)*Mtile + 1 - tile_over,1)
			i2 <- min(ix*Mtile + tile_over,M0.ms) 
			j1 <- max((iy-1)*Ntile + 1 - tile_over,1)
			j2 <- min(iy*Ntile + tile_over,N0.ms)
			
			# read multspectral image subset 
			ijr <- c(i1,i2,j1,j2)

			MS <- read_subset(ms.imagefn,ijr[1],ijr[2],ijr[3],ijr[4])

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
				
				if(TRUE){
					windows()
					par(mfrow=c(1,2))
					image(MSdisp,red=nR,green=nG,blue=nB,axes=TRUE)
					title(main=paste(Satellite," tile x=",ix,"y=",iy,sep=""))
					
					#windows()
					image(MSdisp,attr="ndvi",col=gray((0:255)/255),axes=TRUE)
					title(main=paste(Satellite," NDVI"," tile x=",ix,"y=",iy,sep=""))
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
			magn_thresh <- 1.0e-03
			ndvi_thresh <- 0.30

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

			if(nrow(Blobs)>2000){
				ind <- order(Blobs$magn, decreasing=TRUE)
				ind <- ind[1:2000]
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

			if(nrow(Blobs)>0) Blobs <- Blobs[Blobs$d>=5.0,]
			if(nrow(Blobs)>0) Blobs <- clean_cocentric_blobs(Blobs)

			if(FALSE)if(Show_graph){
				windows()
				par(mfrow=c(1,2))
				image(MSdisp,red=nR,green=nG,blue=nB,axes=TRUE)
				title(main=paste(Satellite," tile x=",ix,"y=",iy,sep=""))
				display_all_blobs(Blobs,"white")
				#if(nrow(Blobs)>0)text(Blobs$x,Blobs$y,labels=1:nrow(Blobs),pos=4,col="green")
				#if(nrow(Blobs)>0)text(Blobs$x,Blobs$y,labels=1:nrow(Blobs),pos=4)
				
				#windows()
				image(MSdisp,attr="ndvi",col=gray((0:255)/255),axes=TRUE)
				title(main=paste(Satellite," NDVI"," tile x=",ix,"y=",iy,sep=""))
				display_all_blobs(Blobs,"red")
				#if(nrow(Blobs)>0)text(Blobs$x,Blobs$y,labels=1:nrow(Blobs),pos=4,col="green")
			}

			if(nrow(Blobs)>0) Large_blobs <- clean_overlap_tree_mask(Blobs) else Large_blobs <- Blobs

			Blobs <- Large_blobs

			if(Show_graph){
				windows()
				par(mfrow=c(1,2))
				image(MSdisp,red=nR,green=nG,blue=nB,axes=TRUE)
				title(main=paste("larger & smaller trees"," tile x=",ix,"y=",iy,sep=""))
				display_all_blobs(Blobs,"white")
				#if(nrow(Blobs)>0)text(Blobs$x,Blobs$y,labels=1:nrow(Blobs),pos=4,col="green")
				#if(nrow(Blobs)>0)text(Blobs$x,Blobs$y,labels=1:nrow(Blobs),pos=4)

				#windows()
				image(MSdisp,attr="ndvi",col=gray((0:255)/255),axes=TRUE)
				title(main=paste(Satellite," NDVI"," tile x=",ix,"y=",iy,sep=""))
				display_all_blobs(Blobs,"red")
				#if(nrow(Blobs)>0)text(Blobs$x,Blobs$y,labels=1:nrow(Blobs),pos=4,col="green")
			}

			setwd(Path_out)
			save(Blobs, file=paste("nonovl_blobs_image_",input,"_P00",strip_id,"_st_area_",oid,"_tx_",ix,"_ty_",iy,".RData",sep=""))

			if(WriteLog) cat(as.character(Sys.time())," tile ",ix," ",iy," successfully processed",nrow(Blobs)," blobs detected","\n")

			All_blobs <- rbind(All_blobs,Blobs)
		}

		All_blobs <- unique(All_blobs)

		Blobs <- All_blobs

		save(Blobs, file=paste("all_blobs_image_",input,"_P00",strip_id,"_st_area_",oid,".RData",sep=""))

		# end_time0 <- Sys.time()
		# end_time0 - start_time

		# Close logfile
		if(WriteLog){
			cat(paste(Sys.time(),"Process ended",sep=" "))
			sink()
		}
	}
}

#==================================================================================
# The End
#==================================================================================
