#==============================================================
# Variable definitions, data import, preparation
#==============================================================
rm(list=ls(all=TRUE))

require(rgdal)
require(Rcpp)

WriteLog <- TRUE

args <- commandArgs(trailingOnly = TRUE)
input <- args[1]
nr_strips <- args[2]

# Slave image: id and strip
# input     <- Sys.getenv("input")
# nr_strips <- as.numeric(Sys.getenv("nr_strips"))

slave_id  <- input

root 	  <- Sys.getenv("baseDG")
Path_lib  <- Sys.getenv("basescript")
Path_in_s <- paste(root, "3_blob_detection", slave_id,sep="/")
Path_out  <- paste(root, "4_image_matching",slave_id,sep="/")
if(!file.exists(Path_out))dir.create(Path_out,showWarnings=FALSE, recursive=TRUE)

setwd(Path_lib)
source("scale_space_lib_v7_1.R")
sourceCpp("matcher_v5.cpp")

setwd(root)
masters_df <- read.csv("masters.csv", header = TRUE, sep = ",", quote = "\"", dec = ".", fill = TRUE, comment.char = "", stringsAsFactors =FALSE)

for(strip_id_s in 1:nr_strips){
	image <- paste(slave_id, "_P00",strip_id_s,"_MUL", ".tif", sep = "")
	imagebase <- unlist(strsplit(image, "[.]"))[1]

	pathmeta  <- paste(root, "0_categ", slave_id, imagebase, sep = "/")

	# Define area
	setwd(pathmeta)
	metadata <- read.csv(paste("metadata_", input, ".csv", sep = ""), header=TRUE, stringsAsFactors=FALSE)

	nr_fixed <- ncol(metadata) - metadata$Nr_study_areas

	for(i_st_area in 1:metadata$Nr_study_areas){
		oid <- metadata[1,nr_fixed+i_st_area]

		ind <- which(masters_df$oid==oid)
		master_id <- masters_df$master_id[ind]
		strip_str_m <- masters_df$master_strip[ind]
		strip_id_m <- as.numeric(substr(strip_str_m,2,4))

		# Babban gona is a separate case where single master is not cvoering all study area
		if(oid=="1060"){
			if(slave_id=="056027177010_01"){
				master_id <- "056027177010_01"
				strip_id_m <- 1
			} else if(slave_id=="056027177020_01"){
				master_id <- "056027177020_01"
				strip_id_m <- 1
			} else if(slave_id=="056027177030_01"){
				master_id <- "056027177020_01"
				strip_id_m <- 1
			} else if(slave_id=="056027177040_01"){
				master_id <- "056027177020_01"
				strip_id_m <- 2
			}
		}
		
		imagebase <- unlist(strsplit(paste(master_id, "_P00",strip_id_m,"_MUL", ".tif", sep = ""), "[.]"))[1]
		pathmeta_master <- paste(root, "0_categ", master_id, imagebase, sep = "/") 
		Path_images <- paste(root,"2_orthorectification_gdal",sep="/")
		Path_in_m <- paste(root, "4_tree_measurement",oid,sep="/")

		setwd(Path_images)

		# Read master and slave raster datasets info
		master_file <- list.files(paste(Path_images,master_id,sep="/"),pattern=paste("P00",strip_id_m,"_st_area_",oid,".tif",sep=""), ignore.case=TRUE)
		ms.imagefn_m <- master_file
		ms.imageinfo_m <- GDALinfo(paste(Path_images,master_id,ms.imagefn_m,sep="/"),silent=TRUE)
		ps.ms_m <- c(ms.imageinfo_m[["res.x"]],ms.imageinfo_m[["res.x"]])

		slave_file <- list.files(paste(Path_images,slave_id,sep="/"),pattern=paste("P00",strip_id_s,"_st_area_",oid,".tif",sep=""), ignore.case=TRUE)
		ms.imagefn_s <- slave_file
		ms.imageinfo_s <- GDALinfo(paste(Path_images,slave_id,ms.imagefn_s,sep="/"),silent=TRUE)
		ps.ms_s <- c(ms.imageinfo_s[["res.x"]],ms.imageinfo_s[["res.x"]])

		#===========================================================================
		# Read image geometry: sun and satellite position (determined from metadata)
		#===========================================================================

		# GET METADATA
		setwd(pathmeta_master)
		metadata <- read.csv(paste("metadata_", master_id, ".csv", sep = ""), stringsAsFactors=FALSE)

		sun_sat <- as.numeric(unlist(metadata[4:8,2]))

		sun_az_m <- sun_sat[1] * pi/180
		sat_az_m <- sun_sat[3] * pi/180
		theta_sun_m <- (90-sun_sat[2]) * pi/180
		theta_sat_m <- (90-sun_sat[4]) * pi/180

		alpha_sun_m <- pi/2 - sun_az_m
		alpha_sat_m <- pi/2 - sat_az_m

		psi_m <- atan2(tan(theta_sat_m)*sin(alpha_sat_m)-tan(theta_sun_m)*sin(alpha_sun_m),tan(theta_sat_m)*cos(alpha_sat_m)-tan(theta_sun_m)*cos(alpha_sun_m))
		corr_factor_m <- cos(psi_m) / (tan(theta_sat_m)*cos(alpha_sat_m)-tan(theta_sun_m)*cos(alpha_sun_m))
		#===========================================================================

		# GET METADATA
		setwd(pathmeta)
		metadata <- read.csv(paste("metadata_", slave_id, ".csv", sep = ""), stringsAsFactors=FALSE, header=TRUE)

		sun_az_s <- metadata$Sun_azimuth * pi/180
		sat_az_s <- metadata$Satellite_azimuth * pi/180
		theta_sun_s <- (90-metadata$Sun_elevation) * pi/180
		theta_sat_s <- (90-metadata$Satellite_elevation) * pi/180

		alpha_sun_s <- pi/2 - sun_az_s
		alpha_sat_s <- pi/2 - sat_az_s

		while(alpha_sun_s<0) alpha_sun_s <- alpha_sun_s + 2*pi
		while(alpha_sat_s<0) alpha_sat_s <- alpha_sat_s + 2*pi

		alpha_shad_s <- alpha_sun_s + pi
		while(alpha_shad_s>=2*pi) alpha_shad_s <- alpha_shad_s - 2*pi
		while(alpha_shad_s<=-2*pi) alpha_shad_s <- alpha_shad_s + 2*pi
		#===========================================================================
		# Open logfile
		if(WriteLog){
			setwd(Path_out)
			sink(file=paste("matching_log_master=",master_id,"_P00",strip_id_m,"_slave=",slave_id,"_P00",strip_id_s,"_st_area_",oid,".txt",sep=""))
			cat(paste(Sys.time(),"starting","\n",sep=" "))
		}

		# Master image info
		N0.ms <- ms.imageinfo_m[["rows"]]
		M0.ms <- ms.imageinfo_m[["columns"]]

		# Read blobs from master image (already corrected for topographic shift)
		setwd(Path_in_m)
		datafile <- paste("corrected_trees_size_shape__P00",strip_id_m,"_st_area_",oid,".RData",sep="")

		if(!file.exists(datafile)){
			if(WriteLog) cat(as.character(Sys.time()),"Error: blobs file for master image not found","\n")
			next
		}

		load(file=datafile)
		#Master_blobs <- Blobs_all
		Blobs_A <- Master_blobs

		# Read blobs from slave image
		setwd(Path_in_s)
		datafile <- paste("all_blobs_image_",slave_id,"_P00",strip_id_s,"_st_area_",oid,".RData",sep="")

		if(!file.exists(datafile)){
			if(WriteLog) cat(as.character(Sys.time()),"Error: blobs file for slave image not found","\n")
			next
		}

		load(file=datafile)
		Blobs_B <- Blobs

		Nb <- ms.imageinfo_m[["bands"]]

		if(Nb==8){
			# Set RGB composition
			nR <- 7
			nG <- 5
			nB <- 3
		}else{
			if(Nb==4){
				nR <- 4
				nG <- 3
				nB <- 2
			}else{
				nR <- 1
				nG <- 1
				nB <- 1
			}
		}

		N0.ms_m <- ms.imageinfo_m[["rows"]]
		M0.ms_m <- ms.imageinfo_m[["columns"]]

		N0.ms_s <- ms.imageinfo_s[["rows"]]
		M0.ms_s <- ms.imageinfo_s[["columns"]]

		Nb_s <- ms.imageinfo_s[["bands"]]
		if(Nb_s==8){
			# Set RGB composition
			nR_s <- 7
			nG_s <- 5
			nB_s <- 3
		}else{
			if(Nb_s==4){
				nR_s <- 4
				nG_s <- 3
				nB_s <- 2
			}else{
				nR_s <- 1
				nG_s <- 1
				nB_s <- 1
			}
		}

		# Shift the master blobs according to slave image geometry
		h  <- Master_blobs$h
		hf <- Master_blobs$hf
		n  <- Master_blobs$n

		R <- 0.5 * Master_blobs$d

		h1 <- h*hf
		h2 <- h-h1

		start_time <- Sys.time()
		topo_shift_s <- project_pollock_quantitative_matrixC(h1, h2, R, theta_sat_s, n)
		Sys.time() - start_time

		Blobs_A$x <- Master_blobs$x - topo_shift_s * cos(alpha_sat_s)
		Blobs_A$y <- Master_blobs$y - topo_shift_s * sin(alpha_sat_s)

		# Define ranges of blobs in master and in slave
		xr_m <- range(Blobs_A$x)
		yr_m <- range(Blobs_A$y)

		xr_s <- range(Blobs_B$x)
		yr_s <- range(Blobs_B$y)

		# Intersection, buffer with the expected search area
		xmin <- max(c(xr_m[1],xr_s[1]))
		xmax <- min(c(xr_m[2],xr_s[2]))

		ymin <- max(c(yr_m[1],yr_s[1]))
		ymax <- min(c(yr_m[2],yr_s[2]))

		r_search <- 40.0

		xmin <- xmin - r_search
		xmax <- xmax + r_search
		ymin <- ymin - r_search
		ymax <- ymax + r_search

		# Discard blobs outside the intersection area
		ind <- which((Blobs_A$x>=xmin)&(Blobs_A$x<=xmax)&(Blobs_A$y>=ymin)&(Blobs_A$y<=ymax))
		Blobs_tile_A <- Blobs_A[ind,]

		ind <- which((Blobs_B$x>=xmin)&(Blobs_B$x<=xmax)&(Blobs_B$y>=ymin)&(Blobs_B$y<=ymax))
		Blobs_tile_B <- Blobs_B[ind,]

		if(nrow(Blobs_tile_A)==0 | nrow(Blobs_tile_B)==0){
			cat(as.character(Sys.time()), paste(" No blobs to match ","\n",sep=""))
			stop("Exiting")
		}

		Blobs_tile_A <- unique(Blobs_tile_A)
		Blobs_tile_B <- unique(Blobs_tile_B)

		# Reduce the number of elements in master and slave.
		# Activate when Rcpp function complains about vector that is too long (negative length)
		# Inputs are fine, but output of mathing candidates generates this error

		# Master
		# max_n_rows <- 32000
		# if(nrow(Blobs_tile_A)>max_n_rows){
			# ind <- order(Blobs_tile_A$magn,decreasing=TRUE)
			# ind <- ind[1:max_n_rows]
			# Blobs_tile_A <- Blobs_tile_A[ind,]
		# }

		# Slave
		max_n_rows <- 32000
		if(nrow(Blobs_tile_B)>max_n_rows){
			ind <- order(Blobs_tile_B$magn,decreasing=TRUE)
			ind <- ind[1:max_n_rows]
			Blobs_tile_B <- Blobs_tile_B[ind,]
		}

		prec <- max(c(mean(ps.ms_s),mean(ps.ms_m)))
		n_samp <- 0.5*min(c(nrow(Blobs_tile_A),nrow(Blobs_tile_B)))

		# Handle errors due to too large input/output matrices

		# consensus_set <- blob_match_v2(Blobs_tile_A,Blobs_tile_B,r_search,prec,n_samp)
		consensus_set <- tryCatch(blob_match_v2(Blobs_tile_A,Blobs_tile_B,r_search,prec,n_samp),
              warning = function(w) {cat(paste("Some warning", "\n"))},
              error = function(e) {cat(paste("Error: matrix is too large", "A", nrow(Blobs_tile_A), "B", nrow(Blobs_tile_A)),"\n");stop()})
 
		n_cons <- nrow(consensus_set)
		# n_cons

		if(n_cons==1 | n_cons==-1) stop("No match found")

		Blobs_A1 <- Blobs_tile_A[consensus_set[,1],]
		Blobs_B1 <- Blobs_tile_B[consensus_set[,2],]

		xydiff <- Blobs_A1[,1:2]-Blobs_B1[,1:2]
		dxy <- colMeans(xydiff)

		dx <- dxy[1]
		dy <- dxy[2]
		# c(n_cons,dx,dy)

		rmse <- sqrt(mean((xydiff[,1]-dx)^2 + (xydiff[,2]-dy)^2))

		newdf <- Blobs_A1
		for(i in 1:5)names(newdf)[i] <- paste(names(newdf)[i],"_m",sep="") 

		newdf <- cbind(newdf,Blobs_B1)
		for(i in (ncol(Blobs_tile_A)+1):(ncol(Blobs_tile_A)+ncol(Blobs_B1)))names(newdf)[i] <- paste(names(newdf)[i],"_s",sep="") 

		Matched_Blobs <- newdf

		if(WriteLog) cat(as.character(Sys.time()),"success", dx,dy,n_cons,rmse,"\n")
		summary_matcher <- c(dx,dy,n_cons,rmse)

		names(summary_matcher) <- c("dx","dy","n_cons","rmse")

		setwd(Path_out)
		save(Matched_Blobs,file=paste("matched_blobs_master_",master_id,"_P00",strip_id_m,"_slave_",slave_id,"_P00",strip_id_s,"_st_area_",oid,".RData",sep=""))
		# save(summary_matcher,file=paste("summary_of_matching",".RData",sep=""))

		# str(Matched_Blobs)
		# summary_matcher

		A <- array(0,c(nrow(Blobs_tile_A),3))
		for(k in 1:3) A[,k] <- Blobs_tile_A[,k]

		B <- array(0,c(nrow(Blobs_tile_B),3))
		for(k in 1:3) B[,k] <- Blobs_tile_B[,k]

		# Subset by area
		# A <- A[(A[,1]>=xmin)&(A[,1]<=xmax)&(A[,2]>=ymin)&(A[,2]<=ymax),]
		# B <- B[(B[,1]+dx>=xmin)&(B[,1]+dx<=xmax)&(B[,2]+dy>=ymin)&(B[,2]+dy<=ymax),]

		# Find possible candidates based on the distance and diameter
		matched_pairs <- match_candC(A, B, dx, dy, 3.0*max(c(mean(ps.ms_m),mean(ps.ms_s))))

		# Discard the outliers
		A2 <- A[matched_pairs[,1],]
		B2 <- B[matched_pairs[,2],]

		d <- sqrt((A2[,1]-B2[,1]-dx)^2 + (A2[,2]-B2[,2]-dy)^2)
		dratio <- pmin(A2[,3],B2[,3]) / pmax(A2[,3],B2[,3])

		dmin <- 3*max(c(mean(ps.ms_m),mean(ps.ms_s)))
		ind <- which(d<dmin & dratio >=0.0)

		A3 <- A2[ind,]
		B3 <- B2[ind,]

		dx_all <- A3[,1] - B3[,1]
		dy_all <- A3[,2] - B3[,2]
		# cut tails:
		dr <- quantile(dx_all,c(0.05,0.95))
		ind <- which(dx_all>=dr[1] & dx_all<=dr[2])
		dx_ref <- dx_all[ind]
		dx_final <- mean(dx_ref,na.rm=TRUE)

		dr <- quantile(dy_all,c(0.05,0.95))
		ind <- which(dy_all>=dr[1] & dy_all<=dr[2])
		dy_ref <- dy_all[ind]
		dy_final <- mean(dy_ref,na.rm=TRUE)

		if(FALSE){
			dx_all <- Matched_Blobs$x_m - Matched_Blobs$x_s
			dy_all <- Matched_Blobs$y_m - Matched_Blobs$y_s

			# remove outliers
			dr <- quantile(dx_all,c(0.05,0.95))
			ind <- which(dx_all>dr[1] & dx_all<dr[2])
			dx_ref <- dx_all[ind]
			dx_final <- mean(dx_ref)

			dr <- quantile(dy_all,c(0.05,0.95))
			ind <- which(dy_all>dr[1] & dy_all<dr[2])
			dy_ref <- dy_all[ind]
			dy_final <- mean(dy_ref)
		}

		# Apply topographic shift of master GCPs
		dx_m <- 0.0
		dy_m <- 0.0

		# Shift master to fit the FMU polygons
		# if(study_area=="ML_Sukumba"){
		if(oid==1000){
			dx_m <- -3.5*ps.ms_m[1]
			dy_m <- +2.0*ps.ms_m[2]
		}

		# if(study_area=="NG_Babban_Gona" & master_id=="056027177020_01" & strip_id_m==2){
		if(oid==1060 & master_id=="056027177020_01" & strip_id_m==2){
			# from the log file: dx=0.50875538912703 dy=1.01785207394032
			dx_m <- 0.50875538912703
			dy_m <- 1.01785207394032
		}

		dx <- dx_final
		dy <- dy_final

		dx_s <- dx + dx_m
		dy_s <- dy + dy_m

		if(WriteLog) cat(paste(Sys.time()," Determined overall transformation parameters: dx=",dx_s," dy=",dy_s,"\n",sep=""))

		# Clean memory
		tmp <- sort(sapply(ls(),function(x){object.size(get(x))}))
		# rm(Master_blobs,Blobs_A,Blobs_B,Blobs_A1,Blobs_B1,Blobs_tile_A,Blobs_tile_B,Blobs,newdf,Matched_Blobs,xydiff,A,B,A2,B2,A3,B3)
		# rm(matched_pairs,consensus_set,topo_shift_s,R,n,hf,h1,h2,h,d,dx_all,dy_all,dx_ref,dy_ref)

		# rm(blob_match,height_shadow,ind,dratio,grid_optim_v4)
		# rm(shadow_tip,interpolate_point_v2,shadow_profile,clean_overlap_v3,clean_overlap,grid_optim_v3)
		# rm(project_pollock,project_pollock_draw, clean_overlap_tree_mask,global_glcm_fast)
		
		ind <- which(tmp>10e3)
		rm(list=names(tmp[ind]))
		
		# Wait: let OS free the memory for R
		# Sys.sleep(30)
		cleanMem <- function(n=10) { for (i in 1:n) gc() }
		# gc()
		cleanMem()

		#====================================================================
		# Transform entire scenes
		#====================================================================
		if(WriteLog) cat(paste(Sys.time()," Applying transformation to slave image","\n",sep=""))

		# This should be done earlier! After ATCOR?
		# for(i_b in 1:Nb_s){
			# ind <- which(Ssub@data[,i_b]==65535)
			# Ssub@data[ind,i_b] <- NA
		# }

		# Using gdal_translate:
		src_file <- paste(Path_images,slave_id,ms.imagefn_s,sep="/")

		text_str <- unlist(strsplit(slave_file, split = "[.]"))[[1]]
		# dst_file <- paste(Path_out,"/",text_str,"_st_area_",oid,"_transformed_topocorr",".tif",sep="")
		dst_file <- paste(Path_out,"/",text_str,"_transformed_topocorr",".tif",sep="")

		command <- paste("gdalinfo",src_file,sep=" ")
		tmp <- system(command, intern = TRUE, ignore.stdout = FALSE, ignore.stderr = FALSE, wait = TRUE, input = NULL)
		str1 <- grep("Upper Left",tmp,value=TRUE)
		str2 <- grep("Lower Right",tmp,value=TRUE)

		str1 <- gsub(","," ",str1)
		str1 <- gsub("\\)"," ",str1)
		tmp <- as.numeric(unlist(strsplit(str1," ")))
		tmp <- tmp[!is.na(tmp)]
		ulx <- tmp[1]
		uly <- tmp[2]

		str2 <- gsub(","," ",str2)
		str2 <- gsub("\\)"," ",str2)
		tmp <- as.numeric(unlist(strsplit(str2," ")))
		tmp <- tmp[!is.na(tmp)]
		lrx <- tmp[1]
		lry <- tmp[2]

		# Apply transformation of coordinates
		ulx <- ulx + dx_s
		uly <- uly + dy_s
		lrx <- lrx + dx_s
		lry <- lry + dy_s

		#-a_ullr ulx uly lrx lry:
		new_corners <- paste(ulx,uly,lrx,lry,sep=" ")
		command <- paste("gdal_translate","-a_ullr",new_corners,src_file,dst_file,sep=" ")
		system(command, intern = FALSE, ignore.stdout = FALSE, ignore.stderr = FALSE, wait = TRUE, input = NULL)

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
