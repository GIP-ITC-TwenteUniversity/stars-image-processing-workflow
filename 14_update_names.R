# THIS FILE PRODUCES STATISTICS FROM THE IMAGES BASED ON THE FMUs

# REMOVE ALL VARIABLES FROM THE WORKSPACE
rm(list=ls(all=TRUE))

##### INSTALL AND LOAD R PACKAGES
ipak <- function(pkg){ # check to see if packages are installed. Install them if they are not, then load them into the R session.
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, repos="http://cran.rstudio.com/", dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}
packages <- c("wkb")
ipak(packages)

Server <- TRUE
EXPORT_CSV <- TRUE

# DEFINE WORKING DIRECTORIES
if(Server==TRUE){
	root <- Sys.getenv('baseDG')
	#root <- "/home/stratouliasd/stars/derived/DG_v8"9# FROM linux352.itc.utwente.nl
}else{
	root <- "S:/derived/DG_v9"
}

pathin <- paste(root, "9_cssl_statistics_backup", sep = "/")
pathout <- paste(root, "9_cssl_statistics", sep = "/")

#========================================================
# READ SPECTRAL TABLES
#========================================================
setwd(pathin)
files <- list.files(pathin, pattern = '*_spectral.csv$', recursive = FALSE)
  
for (i in 1:length(files)) {
  setwd(pathin)

  #========================================================
  # MANIPULATE
  #========================================================
  
  file <- files[i]
  x <- read.csv(file, stringsAsFactors= FALSE)
  
  x[[1]] <- as.character(x[[1]])
  x[[1]] <- paste("0", x[[1]], sep = "") 

  colnames(x)[2] <- "satellite_sensor"
  colnames(x)[6] <- "ndvi"
  colnames(x)[7] <- "sd_ndvi"
  colnames(x)[8] <- "ndvi_green"
  colnames(x)[9] <- "sd_ndvi_green"
  colnames(x)[10] <- "evi"
  colnames(x)[11] <- "sd_evi"
  colnames(x)[12] <- "tcari"
  colnames(x)[13] <- "sd_tcari"
  colnames(x)[14] <- "nir_r"
  colnames(x)[15] <- "sd_nir_r"
  colnames(x)[16] <- "sarvi"
  colnames(x)[17] <- "sd_sarvi"
  colnames(x)[18] <- "savi"
  colnames(x)[19] <- "sd_savi"
  colnames(x)[20] <- "msavi"
  colnames(x)[21] <- "sd_msavi"
  colnames(x)[22] <- "wv2vi"
  colnames(x)[23] <- "sd_wv2vi"

  #========================================================
  # EXPORT UPDATED SPECTRAL TABLES
  #========================================================
  
  setwd(pathout)
  if(EXPORT_CSV==TRUE){
  write.csv(x, file = file, row.names=FALSE)
  }
}

#============================================================
# The END
#============================================================
