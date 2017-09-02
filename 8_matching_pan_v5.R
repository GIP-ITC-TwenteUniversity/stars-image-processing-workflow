#======================================================================================
# Variable definitions, data import, preparation
#======================================================================================
rm(list=ls(all=TRUE))

require(rgdal)
# require(Rcpp)

WriteLog <- TRUE

# Slave image: id and strip
args <- commandArgs(trailingOnly = TRUE)
input <- args[1]
nr_strips <- args[2]

# input     <- Sys.getenv("input")
# nr_strips <- as.numeric(Sys.getenv("nr_strips"))

slave_id  <- input

root      <- Sys.getenv("baseDG")
Path_lib  <- Sys.getenv("basescript")
Path_out  <- paste(root,"4_image_matching_pan", slave_id, sep="/")
if(!file.exists(Path_out))dir.create(Path_out,showWarnings=FALSE, recursive=TRUE)

setwd(Path_lib)
source("scale_space_lib_v7_1.R")
#sourceCpp("matcher_v5.cpp")

setwd(root)
masters_df <- read.csv("masters.csv", header = TRUE, sep = ",", quote = "\"", dec = ".", fill = TRUE, comment.char = "", stringsAsFactors =FALSE)

for(strip_id_s in 1:nr_strips){
	image <- paste(slave_id, "_P00",strip_id_s,"_MUL", ".tif", sep = "")
	imagebase <- unlist(strsplit(image, "[.]"))[1]

	pathmeta <- paste(root, "0_categ", slave_id, imagebase, sep = "/")

	# Define area
	setwd(pathmeta)
	metadata <- read.csv(paste("metadata_", slave_id, ".csv", sep = ""), header=TRUE, stringsAsFactors=FALSE)

	nr_fixed <- ncol(metadata) - metadata$Nr_study_areas

	for(i_st_area in 1:metadata$Nr_study_areas){
		oid <- metadata[1,nr_fixed+i_st_area]

		ind <- which(masters_df$oid==oid)
		master_id <- masters_df$master_id[ind]
		strip_str_m <- masters_df$master_strip[ind]
		strip_id_m <- as.numeric(substr(strip_str_m,2,4))


		Path_in <- paste(root,"2_orthorectification_gdal_pan", slave_id ,sep="/")
		setwd(Path_in)

		# Read slave raster info
		slave_file <- list.files(Path_in,pattern=paste("P00",strip_id_s,"_st_area_",oid,".tif",sep=""), ignore.case=TRUE)
		ms.imageinfo_s <- GDALinfo(paste(Path_in, slave_file, sep="/"),silent=TRUE)
		ps.ms_s <- c(ms.imageinfo_s[["res.x"]],ms.imageinfo_s[["res.x"]])

		# Open logfile
		if(WriteLog){
			setwd(Path_out)
			sink(file=paste("matching_log_pan_master=",master_id,"_P00",strip_id_m,"_slave=",slave_id,"_P00",strip_id_s,"_st_area_",oid,".txt",sep=""))
			cat(paste(Sys.time(),"starting","\n",sep=" "))
		}

		if(length(slave_file)!=1){
			cat(paste("Error: pan image is not found or is not unique"),"\n",sep="")
			if(WriteLog){
				cat(paste(Sys.time(),"Process failed","\n",sep=" "))
				sink()
			}
			next
		}

		# Set input & output directories
		Path_in_ms_match <- paste(root, "4_image_matching", slave_id, sep="/")

		if(!file.exists(Path_in_ms_match)){
			cat(paste("Error: this file was not matched"),"\n",sep="")
			if(WriteLog){
				cat(paste(Sys.time(),"Process failed","\n",sep=" "))
				sink()
			}
			next
		}

		setwd(Path_in_ms_match)

		# Read log file
		all_text <- readLines(paste("matching_log_master=",master_id,"_P00",strip_id_m,"_slave=",slave_id,"_P00",strip_id_s,"_st_area_",oid,".txt",sep=""),-1,warn=FALSE)
		n_lines <- length(all_text)

		# Check if the process completed
		if(length(grep("Process ended",all_text))!=1){
			cat(paste("Error: this file was not matched. Quitting now."),"\n",sep="")
			if(WriteLog){
				cat(paste(Sys.time(),"Process failed","\n",sep=" "))
				sink()
			}
			next
		}

		if(grep("Process ended",all_text)!=n_lines){
			cat(paste("Error: this file was not matched. Quitting now."),"\n",sep="")
			if(WriteLog){
				cat(paste(Sys.time(),"Process failed","\n",sep=" "))
				sink()
			}
			next
		}

		# Read transformation parameters
		text_str <- all_text[grep("Determined overall transformation parameters:",all_text)]
		str_arr <- unlist(strsplit(text_str," "))
		n_str <- length(str_arr)
		str_x <- str_arr[n_str-1]
		str_y <- str_arr[n_str]
		str_x2 <- unlist(strsplit(str_x,"="))
		dx <- as.numeric(str_x2[2])
		str_y2 <- unlist(strsplit(str_y,"="))
		dy <- as.numeric(str_y2[2])

		dx_m <- 0.0
		dy_m <- 0.0

		# Apply topographic shift of master GCPs
		# Now the master has been alterred
		#dx=0.581699130823836 dy=1.27540983469225
		# dx_m <- 0.5817
		# dy_m <- 1.2754

		#====================================================================
		# Transform entire scene
		#====================================================================
		if(WriteLog) cat(paste(Sys.time()," Applying transformation to slave image","\n",sep=""))
		if(WriteLog) cat(paste(Sys.time()," Transformation parameters: dx=",dx,", dy=",dy,"\n",sep=""))
		dx_s <- dx + dx_m
		dy_s <- dy + dy_m

		# Using gdal_translate:
		src_file <- paste(Path_in,"/",slave_file,sep="")

		text_str <- unlist(strsplit(slave_file, split = "[.]"))[[1]]
		dst_file <- paste(Path_out,"/",text_str,"_gdal_transformed_topocorr",".tif",sep="")

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
			cat(paste(Sys.time(),"Process ended","\n",sep=" "))
			sink()
		}
	}
}
#==================================================================================
# The End
#==================================================================================
