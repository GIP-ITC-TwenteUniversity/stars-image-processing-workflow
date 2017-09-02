#================================================================
args <- commandArgs(trailingOnly = TRUE)
input <- args[1] 
# input <- Sys.getenv("input")

#================================================================
# DEFINE WORKING DIRECTORIES
#================================================================
root <- Sys.getenv("baseDG")
stars_dir <- Sys.getenv("stars")
path_in <- paste(stars_dir,"acquired/DG", input, sep="/")
setwd(path_in)
#================================================================
path_out <- paste(root, "0_categ", input, sep = "/")
if(!file.exists(path_out))dir.create(path_out,showWarnings=FALSE, recursive=TRUE, mode = "0777")

filenames <- list.files(pattern=".XML", ignore.case=TRUE, full.names = FALSE, recursive = FALSE, include.dirs = FALSE)

setwd(path_out)
missing_files <- !file.exists(filenames)
filenames <- filenames[missing_files]
setwd(path_in)
file.copy(filenames, to = path_out, recursive = TRUE, overwrite = TRUE)

sub_dirs <- list.dirs(path = path_in, full.names = TRUE, recursive = FALSE)

for(sub_id in seq_along(sub_dirs)){

	sub_dir <- sub_dirs[sub_id]
	setwd(sub_dir)

	sub_dir_out <- paste(path_out, basename(sub_dir), sep = "/")
	if(!file.exists(sub_dir_out))dir.create(sub_dir_out,showWarnings=FALSE, recursive=TRUE, mode = "0777")

	filenames <- list.files(full.names = FALSE, recursive = TRUE, include.dirs = FALSE) 

	# exclude tif files
	ind <- grep("*.tif", filenames, ignore.case=TRUE, fixed=FALSE, value=FALSE)

	if(length(ind)>0) sel_files <- filenames[-ind] else sel_files <- filenames

	setwd(sub_dir_out)
	missing_files <- !file.exists(sel_files)
	sel_files <- sel_files[missing_files]
	setwd(sub_dir)
	file.copy(sel_files,  to = sub_dir_out, recursive = TRUE, overwrite = TRUE)
}
