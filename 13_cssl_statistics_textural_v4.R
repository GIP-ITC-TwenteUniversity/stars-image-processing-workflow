# Extract texture and rows from image using FMU polygons
rm(list=ls(all=TRUE))

# INSTALL AND LOAD R PACKAGES
ipak <- function(pkg){ # check to see if packages are installed. Install them if they are not, then load them into the R session.
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, repos="http://cran.rstudio.com/", dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}
packages <- c("rgdal", "rgeos", "RPostgreSQL", "wkb", "Rcpp")
ipak(packages)

Server <- TRUE # run on server
EXPORT_CSV <- TRUE # export csv files
EXPORT_POSTGRESQL_table <- FALSE # export to POSTGRESQL as table

Show_graph <- FALSE
#============================================================
# Texture computation parameters
#============================================================
# How many grey levels should be considered in GLCM analysis?
ql <- 256 # 64 or 256
# Dependence on lag (h)
h_arr <- 1:12
# Dependence on angle (\phi)
nphi <- 4
phi_arr <- seq(0,pi - pi/nphi,length.out=nphi)
# A very small number
eps <- 1e-32


##### DEFINE WORKING DIRECTORIES
if(Server){
	root <- Sys.getenv('baseDG')
	Path_lib <- Sys.getenv('basescript')
} else{
	root <- "S:/derived/DG_v9"
	Path_lib <- "S:/derived/scripts/v9"
}

setwd(Path_lib)
source("scale_space_lib_v7_1.R")
sourceCpp("glcm_v3.cpp")

# source("tree_measurement_lib_v9.R")
# sourceCpp("matcher_v5.cpp")
# sourceCpp("blobs_lib_v4.cpp")
# sourceCpp("tree_mask_from_blobs_v3.cpp")
#========================================================
# IMAGE INPUT
#========================================================
args <- commandArgs(trailingOnly = TRUE)
input <- Sys.getenv('input') # input="054112895070_01"

pan.imagefn <- paste(input, "_P001_MUL", ".tif", sep = "")
imagebase <- unlist(strsplit(pan.imagefn, "[.]"))[1]

pathin_ms <- paste(root, "6_image_matching", input, sep = "/")
pathin_pan <- paste(root, "6_image_matching_pan", input, sep = "/")
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
year <- as.integer(substr(date, 1, 4))
erosion <- as.numeric(2.0) # it does not work at the moment, erosion is set in the query line 
sensor_name <- metadata[2,2]
study_area_OID <- as.integer(unlist(metadata[14,2]))

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

tex_names <- c("asm","contrast","correlation","variance","homogeneity","sum_avg","sum_var","sum_entropy","entropy","diff_var","diff_entropy","infor_corr1","infor_corr2","max_cor_coef")

haralick_textures_v3 <- function(A){
	tex_names <- c("asm","contrast","correlation","variance","homogeneity","sum_avg","sum_var","sum_entropy","entropy","diff_var","diff_entropy","infor_corr1","infor_corr2","max_cor_coef")

	eps <- 1.0e-32

	Ntex <- length(tex_names)

	out <- array(NA,Ntex)
	names(out) <- tex_names

	ql <- dim(A)[1]
	if(dim(A)[2]!=ql) stop("Not square GCLM")

	rownames(A) <- colnames(A) <- 0:(ql-1)

	ind <- which(rowSums(A)>0)

	S <- A[ind,ind,drop=F]

	levels_x <- levels_y <- as.numeric(rownames(S))

	rows_s <- row(A)-1
	rows_s <- rows_s[ind,ind]
	cols_s <- col(A)-1
	cols_s <- cols_s[ind,ind]

	px <- colSums(A)
	py <- rowSums(A)

	px_s <- colSums(S)
	py_s <- rowSums(S)

	#p_(x+y)(k) in Haralick's paper
	p_plus_arr <- as.numeric(p_plus_arr_c_v2(A))

	iarr <- 2:(2*ql)
	ind <- which(p_plus_arr>0)
	iarr_s <- iarr[ind]
	p_plus_arr_s <- p_plus_arr[ind]

	#p_(x-y)(k) in Haralick's paper
	#p_minus_arr <- as.numeric(p_minus_arr_c_v2(A))
	# Remark: computed discarding zero rows and columns!
	p_minus_arr <- as.numeric(p_minus_arr_c_v2(S))
	p_minus_arr_s <- p_minus_arr[p_minus_arr>0]

	# f1 in Haralick et al. (1973)
	# Angular Second Moment
	# Note: some sources apply the square root to get energy from ASM
	#out[1] <- sum(A^2)
	out[1] <- sum(S^2)

	# f2 in Haralick et al. (1973)
	# formula is from elsewhere, simplified computationally
	out[2] <- sum(S*((rows_s-cols_s)^2))

	# f3 in Haralick et al. (1973)
	# Correlation
	mx <- sum(px_s*levels_x)
	my <- sum(py_s*levels_y)

	sx <- sqrt(sum(px_s*((levels_x-mx)^2)))
	sy <- sqrt(sum(py_s*((levels_y-my)^2)))

	if(sx*sy>eps){
		val <- sum((rows_s-my)*(cols_s-mx)*S)
		out[3] <- val / (sx*sy)
	}

	# f4 in Haralick et al. (1973)
	# Sum of Squares: variance
	# Note: m is not defined in the source. Use mx
	# Note: the summation is from 0 to ql-1, not from 1 to ql
	out[4] <- sum(S*((cols_s-mx)^2))

	# f5 in Haralick et al. (1973)
	# Inverse Difference Moment
	out[5] <- sum(S/(1+((rows_s-cols_s)^2)))

	# f6 in Haralick et al. (1973)
	# Sum Average
	#out[6] <- sum(p_plus_arr*iarr)
	out[6] <- sum(p_plus_arr_s*iarr_s)

	# f7 in Haralick et al. (1973)
	# Sum Variance
	mu <- out[6]
	#out[7] <- sum(((iarr-mu)^2) * p_plus_arr)
	out[7] <- sum(((iarr_s-mu)^2) * p_plus_arr_s)
	
	# f8 in Haralick et al. (1973)
	# Sum Entropy
	# Remark: different log bases are found in literature
	# Follow base 2 here
	#val <- p_plus_arr * log2(p_plus_arr)
	#out[8] <- -sum(val[p_plus_arr>0])
	out[8] <- -sum(p_plus_arr_s * log2(p_plus_arr_s))

	# f9 in Haralick et al. (1973)
	# Remark: different log bases are found in literature
	# Follow base 2 here
	val <- as.vector(S)
	val <- val[val>0]
	val <- - val*log2(val)
	out[9] <- sum(val)

	# f10 in Haralick et al. (1973)
	# Difference variance
	out[10] <- var(p_minus_arr_s)

	# f11 in Haralick et al. (1973)
	# Difference entropy
	val <- p_minus_arr_s
	out[11] <- -sum(val*log2(val))

	# f12 in Haralick et al. (1973)
	# Information measures of correlation (1)
	HXY <- out[9]
	#tmp <- f12(A,px,py,HXY)
	tmp <- f12(S,px_s,py_s,HXY)
	out[12] <- tmp$f12
	# f13 in Haralick et al. (1973)
	# Information measures of correlation (2)
	out[13] <- tmp$f13

	# f14 in Haralick et al. (1973)
	# Maximal correlation coefficient
	Q <- Qfun(S,px_s,py_s)
	out[14] <- Re(eigen(Q,only.values=TRUE)$values[2])

	return(out)
}

# RULES TO CONNECT  

if(TRUE){
	# ATTEMPT 1 WITH RGDAL - FAIL NO POSTGIS DRIVER INSTALLED
	# connect to postgres
	#dsn="PG:dbname='cssl' host=linux352.itc.utwente.nl port=5432 user=dimitris password=random_pwd"
	#ogrDrivers() # list available drivers
	#ogrListLayers(dsn) #list tables in the database
	#commonfmu = readOGR(dsn = dsn, "common.fmu") #read a table

	# ATTEMPT 2 WITH RPostgreSQL
	drv <- dbDriver("PostgreSQL") # load PostgreSQL driver
	con <- dbConnect(drv, host='linux352.itc.utwente.nl', port='5432', dbname='cssl', user='dimitris', password='random_pwd') # open a connection
	#dbListTables(con) # list existing tables

	# 1. generic query
	# query <- "SELECT * from common.fmu"
	# query <- "SELECT target_of_interest_oid AS id, ST_AsBinary(geometry) AS geom from common.fmu" # ST_AsText or ST_AsBinary
	# rs <- dbSendQuery(con, query) # submit a query
	# fmu<- fetch(rs,n=-1)  # fetch all elements from the result set

	# 2. specific query
	#query_s <- "SELECT common.cs_fmu_studyarea_year(1005, 2015, 0.0)" # Mali 1000, Nigeria 1005
	# query <- "SELECT id AS id, ST_AsBinary(geom) AS geom from common.cs_fmu_studyarea_year(1005,2015,2.0)"
	# rs <- dbSendQuery(con, query) # submit a query
	# fmu <- fetch(rs,n=-1)  # fetch all elements from the result set
	
	# 3. specific query with prepared statement
	query <- postgresqlExecStatement(con, "SELECT id AS id, ST_AsBinary(geom) AS geom from common.cs_fmu_studyarea_year($1,$2,2.0)", c(study_area_OID, year_OID))
	fmu <- fetch(query,n=-1)

	# 3. test
	#query_s <- "SELECT common.cs_fmu_studyarea_year(1005, 2015, 0.0)" # Mali 1000, Nigeria 1005
	# query <- "SELECT target_of_interest_oid AS id, ST_AsBinary(geometry) AS geom from common.fmu" # ST_AsText or ST_AsBinary
	# rs <- dbSendQuery(con, query) # submit a query
	# fmu<- fetch(rs,n=-1)  # fetch all elements from the result set
	# str(fmu)

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
}else{
	# Just read shapefile instead
	# Read FMU shapefile
	if(Server){
		#Path_fmu <- paste("/home/tolpekin/R","/Input/Nigeria",sep="")
		Path_fmu <- paste("/home/tolpekin/R","/Input/Mali/mali_field",sep="")
	}else{
		#Path_fmu <- paste("D:/Programming/STARS","/Input/Nigeria",sep="")
		#Path_fmu <- paste("D:/Programming/STARS","/Input/Mali/mali_field",sep="")
		Path_fmu <- paste("D:/Programming/STARS","/texture",sep="")
	}
	#fn.polygons <- "nigeria_utm"
	fn.polygons <- "Mali_fmu_utm"
	#fn.polygons <- "mali2014full"
	Pol <- readOGR(dsn=Path_fmu, layer=fn.polygons)
	proj_vector <- Pol@proj4string@projargs
	n.pol <- nrow(Pol)

	names(Pol)[1] <- "id"

	vector_postgresql <- Pol
	rm(Pol)
}
#================================================================================
# raster
#================================================================================
setwd(pathin_pan)
image_list <- list.files(pattern = paste(input, ".*TIF$", sep = ""), ignore.case = TRUE)

if(length(image_list)==0){
	stop("Pan image not found in the directory")
}

pan.imagefn <- image_list[1]
sensor <- metadata[2,2]
imagebase <- unlist(strsplit(pan.imagefn, "[.]"))
folder_name <- substr(pan.imagefn, 0, 12)

pan.imageinfo <- GDALinfo(paste(pathin_pan,pan.imagefn,sep="/"), silent=TRUE)
ps.pan <- c(pan.imageinfo[["res.x"]],pan.imageinfo[["res.x"]])
N0.pan <- pan.imageinfo[["rows"]]
M0.pan <- pan.imageinfo[["columns"]]

# Determine data range
# Not working (-32768,32767). Fix this.
if(FALSE){
	#pan_min <- attributes(pan.imageinfo)$df$Bmin
	#pan_max <- attributes(pan.imageinfo)$df$Bmax
	
	#pan_min <- 0
	#pan_max <- 650

} else{
	#pan_min <- 0
	#pan_max <- 2047
	A <- readGDAL(paste(pathin_pan,pan.imagefn,sep="/"), silent=TRUE)
	y <- A@data[,1]
	tmp <-quantile(y,c(0.01,0.99))
	tmp <- unname(tmp, force = TRUE)
	pan_min <- tmp[1]
	pan_max <- tmp[2]
	rm(y,A, tmp)
}

proj_raster <- attributes(pan.imageinfo)$projection
xyLL.pan <- c(pan.imageinfo[["ll.x"]],pan.imageinfo[["ll.y"]])
ysign <- attributes(pan.imageinfo)$ysign #y sign

setwd(pathin_ms)
image_list <- list.files(pattern = paste(input, ".*TIF$", sep = ""), ignore.case = TRUE)
if(length(image_list)==0){
	stop("MS image not found in the directory")
}

ms.imagefn <- image_list[1]
sensor <- metadata[2,2]
imagebase <- unlist(strsplit(ms.imagefn, "[.]"))
folder_name <- substr(ms.imagefn, 0, 12)

ms.imageinfo <- GDALinfo(paste(pathin_ms,ms.imagefn,sep="/"), silent=TRUE)
ps.ms <- c(ms.imageinfo[["res.x"]],ms.imageinfo[["res.x"]])
N0.ms <- ms.imageinfo[["rows"]]
M0.ms <- ms.imageinfo[["columns"]]
Nb 	  <- ms.imageinfo[["bands"]]

xyLL.ms <- c(ms.imageinfo[["ll.x"]],ms.imageinfo[["ll.y"]])
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
if(proj_raster != proj4string(vector_postgresql)){
	vector_postgresql <- spTransform(vector_postgresql, proj_raster)
}
# ================================================================================
# Extract raster clip-polygons and statistics
# ================================================================================
setwd(pathout)

n_pol <- nrow(vector_postgresql)

Nc <- array(0,n_pol) #create null array

Ntex <- length(tex_names)

n_overhead <- 5
n_tex <- Ntex*length(h_arr)*nphi

col_num <- n_overhead + n_tex	#number of columns
pan_stats <- data.frame(array(NA,c(n_pol,col_num)),stringsAsFactors = FALSE)

# Overhead
names_overhead <- c("image_identifier","satellite_sensor","target_of_interest_oid","number_of_all_pixels","number_of_contributing_pixels")
# Textures
names_tex <- array("",n_tex)
i_col <- 1
for(i in h_arr){
	for(j in phi_arr){
		for(k in 1:Ntex){
			names_tex[i_col] <- paste(tex_names[k],"_h",i,"_phi",round(j*180/pi,digits=1),sep="")
			i_col <- i_col + 1
		}
	}
}

names(pan_stats) <- c(names_overhead,names_tex)

for(i.pol in 1:n_pol){
	field_id <- vector_postgresql@data$id[[i.pol]]

	# read shapefiles marginal coordinates
	# In case of complex polygons (with holes)
	n_subpol <- length(vector_postgresql@polygons[[i.pol]])

	pol <- vector_postgresql@polygons[[i.pol]]@Polygons[[1]]
	xy_pol <- pol@coords
	
	xy_pol[,1] <- xy_pol[,1] + (-6.999)
	xy_pol[,2] <- xy_pol[,2] + (4.013)
	
	xrl <- range(xy_pol[,1])
	yrl <- range(xy_pol[,2])

	# Add margins to ensure sufficiently large dataset
	margin <- 1.0 * mean(ps.ms)
	#margin <- 20.0 * mean(ps.ms)
	xrl <- xrl + c(-1,1) * margin
	yrl <- yrl + c(-1,1) * margin

	# read multspectral image subset 
	ijr_ms <- xy_to_rowcol(cbind(xrl,yrl),ms.imageinfo)

	if((ijr_ms[1]>M0.ms)||(ijr_ms[2]<=1)||(ijr_ms[3]>N0.ms)||(ijr_ms[4]<=1)){ # Polygon is  completely outside the image
		pan_stats[i.pol,1:n_overhead] <- c(folder_name, sensor, field_id, 0, 0)
		next
	}

	if(ijr_ms[1]<1) ijr_ms[1] <- 1
	if(ijr_ms[2]>M0.ms) ijr_ms[2] <- M0.ms
	if(ijr_ms[3]<1) ijr_ms[3] <- 1
	if(ijr_ms[4]>N0.ms) ijr_ms[4] <- N0.ms

	Path_in <- pathin_ms
	MS <- read_subset(ms.imagefn,ijr_ms[1],ijr_ms[2],ijr_ms[3],ijr_ms[4])

	# If all pixels are NA, skip
	sums <- rowSums(MS@data)
	ind <- which(!is.na(sums))
	if(length(ind)==0){
		pan_stats[i.pol,1:n_overhead] <- c(folder_name, sensor, field_id, 0, 0)
		next
	}
	
	# If all pixels are 0
	if(max(MS@data) == 0 & min(MS@data) == 0){
		pan_stats[i.pol,1:n_overhead] <- c(folder_name, sensor, field_id, 0, 0)
		next
	}
	
	M <- MS@grid@cells.dim[1]
	N <- MS@grid@cells.dim[2]

	xrl <- bbox(MS)[1,]+0.5*c(1,-1)*ps.pan
	yrl <- bbox(MS)[2,]+0.5*c(1,-1)*ps.pan

	# read pan image subset 
	ijr_pan <- xy_to_rowcol(cbind(xrl,yrl),pan.imageinfo)

	if(ijr_pan[1]<1) ijr_pan[1] <- 1
	if(ijr_pan[2]>M0.pan) ijr_pan[2] <- M0.pan
	if(ijr_pan[3]<1) ijr_pan[3] <- 1
	if(ijr_pan[4]>N0.pan) ijr_pan[4] <- N0.pan

	Path_in <- pathin_pan
	Pan <- read_subset(pan.imagefn,ijr_pan[1],ijr_pan[2],ijr_pan[3],ijr_pan[4])

	if(Show_graph){
		nR <- 7
		nG <- 5
		nB <- 3
		MSdisp <- MS
		for(k in 1:Nb)MSdisp@data[,k] <- histstretch(MS@data[,k])
	}
	if(Show_graph){
		Pandisp <- Pan
		k <- 1
		Pandisp@data[,k] <- histstretch(Pan@data[,k])
	}
	if(Show_graph){
		pol <- vector_postgresql@polygons[[i.pol]]@Polygons[[1]]

		windows()
		par(mfrow=c(1,2))
		image(MSdisp,red=nR,green=nG,blue=nB,axes=TRUE)
		title(main=paste("polygon ",i.pol,sep=""))
		polygon(xy_pol,border="blue",lwd=2)

		image(Pandisp,col=gray((0:255)/255),axes=TRUE)
		polygon(xy_pol,border="blue",lwd=2)
	}

	xy.ms <- coordinates(MS)
	xy.pan <- coordinates(Pan)

	xyTL.ms  <- xy.ms[1,]
	xyTL.pan <- xy.pan[1,]

	# Read and apply tree & tree shadow mask
	setwd(path_tree_mask)
	Path_in <- path_tree_mask
	filename <- paste("Tree_mask_pan_",input,".tif",sep="")
	if(!file.exists(filename)){
		stop("Tree mask not found")
	}

	# Careful! Assuming that tree mask raster has identical grid with the multispectral image 
	tree_mask <- read_subset(filename,ijr_pan[1],ijr_pan[2],ijr_pan[3],ijr_pan[4])
	ind_mask <- which(tree_mask$band1>0)

	if(length(ind_mask)>0) Pan@data[ind_mask,] <- NA

	# Read and apply cloud mask
	setwd(path_cloud_mask)
	Path_in <- path_cloud_mask
	filename <- paste(input,"_cloud_mask",".tif",sep="")

	if(!file.exists(filename)){
		stop("Cloud mask not found")
	}
	
	# Careful! Assuming that cloud mask raster has identical grid with the multispectral image 
	cloud_mask <- read_subset(filename,ijr_ms[1],ijr_ms[2],ijr_ms[3],ijr_ms[4])

	ind_mask <- which(cloud_mask$band1==0)
	if(length(ind_mask)>0){
		# Convert tree mask to pan grid
		ind_pan <- get_pan_pixels(ind_mask-1, M, N) + 1	# note +1 and -1 to convert sequence 1,...,n to 0,...,n-1
		Pan@data[ind_pan,] <- NA
	}

	# SUBSETTING BY POLYGON
	# Remove holes!
	xy <- coordinates(Pan)
	contrib_pix <- array(0,nrow(Pan))
	
	tmp.in <- point.in.polygon(xy[,1],xy[,2],xy_pol[,1],xy_pol[,2])
	contrib_pix[tmp.in==1] <- 1
	
	if(n_subpol > 1){
		for(i_sub in 2:n_subpol){
			if(vector_postgresql@polygons[[i.pol]]@Polygons[[i_sub]]@hole){
				h_pol <- vector_postgresql@polygons[[i.pol]]@Polygons[[i_sub]]
				xy_hole <- h_pol@coords
				tmp.in <- point.in.polygon(xy[,1],xy[,2],xy_hole[,1],xy_hole[,2])
				contrib_pix[tmp.in==1] <- 0
			}
		}
	}

	ind <- which(contrib_pix==1)
	n_pix <- length(ind)
	fmu_values <- Pan@data[ind,]

	ind <- which(contrib_pix==0)
	Pan@data[ind,] <- NA

	ind <- which(is.na(fmu_values))
	pix_NA <- length(ind) # number of NAs
	# Remove NAs
	fmu_values <- fmu_values[!is.na(fmu_values)]
	
	pix_zeros <- sum(fmu_values==0L, na.rm = TRUE) # number of zero pixels
	pix_contr <- length(fmu_values) - pix_zeros # number of contributing pixels - not equal to zero

	if(Show_graph){
		pol <- vector_postgresql@polygons[[i.pol]]@Polygons[[1]]

		windows()
		par(mfrow=c(1,2))
		#image(MSdisp,red=nR,green=nG,blue=nB,axes=TRUE)
		image(Pan,col=gray((0:255)/255),axes=TRUE)
		polygon(xy_pol,border="blue",lwd=2)

		image(Pandisp,col=gray((0:255)/255),axes=TRUE)
		polygon(xy_pol,border="blue",lwd=2)
	}

	# Texture analysis
	#Quantize the image into ql grey levels
	i_col <- n_overhead + 1

	M <- Pan@grid@cells.dim[1]
	N <- Pan@grid@cells.dim[2]

	z <- Pan@data[,1]
	zq <- z

	z[z<pan_min] <- pan_min
	z[z>pan_max] <- pan_max

	breaks <- seq(from = pan_min, to = pan_max, length.out = ql+1)

	zq <- cut(z,breaks,include.lowest=TRUE,labels=FALSE)
	zq <- zq-1

	qdata <- zq

	bad_val <- -1
	qdata[is.na(qdata)] <- bad_val

	P <- array(qdata,c(M,N))

	for(i in 1:length(h_arr))
	for(j in 1:length(phi_arr)){
		lag_dist <- h_arr[i]
		phi <- phi_arr[j]

		dx <- round(lag_dist*cos(phi))
		dy <- round(lag_dist*sin(phi))

		GLCM <- glcmC(P, ql, dx, dy, bad_val)

		GLCM <- GLCM/sum(GLCM)

		#GLCM <- 0.5*(GLCM+t(GLCM))

		if(all(is.na(GLCM))){
			i_col <- i_col + Ntex
			next
		}

		tex <- haralick_textures_v3(GLCM)
		pan_stats[i.pol,i_col:(i_col+Ntex-1)] <- tex
		i_col <- i_col + Ntex
	}

	pan_stats[i.pol,1:n_overhead] <- c(folder_name, sensor_name, field_id, n_pix, pix_contr)
}

setwd(pathout)
if(EXPORT_CSV==TRUE){
filename <- substr(imagebase[1], 1, nchar(imagebase[1])-21)
write.csv(pan_stats, file = paste(filename, "_textural_ql", ql, ".csv", sep = ""), row.names=FALSE)}

#============================================================
# The END
#============================================================
