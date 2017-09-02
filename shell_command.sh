#!/bin/bash
#bash
#nohup ./shell_command.sh &>/dev/null &
#cd $base/scripts/v10
# nohup Rscript not_in_WF_execute_workflow.R &>/dev/null &

# SET PATHNAME 
#This is non-permanent change, LD_LIBRARY_PATH reset to default value at next session)
# export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib # extends the session's library-path
# export LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu:/usr/local/lib:/usr/lib
export stars=~/stars
export base=~/stars/derived
export baseDG=~/stars/derived/DG_v10
export basescript=~/stars/derived/scripts/v10
export R_LIBS=~/stars/rlibs

# SET IMAGE FILENAME
export input=$1
# export input="054330851010_01" # BG ortho shift 
# export input="054393058010_01" # BG ortho shift
# export input="054294362010_01" # Nigeria cssl issue
# export input="053828840150_01" # Nigeria (Kofa) OR test
# export input="053734892380_01" # Nigeria master, standard
# export input="053613698020_01" # Mali
export input="056027177010_01" # NG Babban Gona image 1
# export input="056027177020_01" # NG Babban Gona image 2
# export input="056027177030_01" # NG Babban Gona image 3
# export input="056027177040_01" # NG Babban Gona image 4
# export input="053734892020_01" 
# export input="054112895030_01" # Mali Sukumba master
# export input="054123995110_01" # BG master
# export input="054551817010_01" # TZ Kilosa master
# export input="054448707010_01" # TZ Njmoby master
# export input="054157627140_01" # TZ Same master
# export input="053734892310_01" # UG Moroto master

# export input="053828840120_01"
# export input="054112895030_01"
echo $input

# Determine number of strips and filenames
path_origin=$stars"/acquired/DG/"$input
cd $path_origin
shopt -s nullglob
strip_array=(*_P00*_MUL)
shopt -u nullglob # Turn off nullglob to make sure it doesn't interfere with anything later
export nr_strips=${#strip_array[@]}
# echo $nr_strips



# IMAGE MATCHING
Rscript $basescript/8_matching_v5.R $input $nr_strips
Rscript $basescript/8_matching_pan_v5.R $input $nr_strips















# COPY DATA FROM ARCHIVE TO PROCESSING NODE
Rscript $basescript/0_copy_data.R $input

# DERIVE METADATA
Rscript $basescript/1_metadata_v2.R $input $nr_strips

# CREATE IMAGE-SPECIFIC DIRECTORIES
mkdir -p $baseDG/1_atcor_6s/"$input"
mkdir -p $baseDG/2_orthorectification_gdal/"$input"
mkdir -p $baseDG/3_blob_detection/"$input"
mkdir -p $baseDG/4_image_matching/"$input"
mkdir -p $baseDG/5_cloud_mask/"$input"
mkdir -p $baseDG/6_tree_mask/"$input"

# Split the process into the relevant and irrelevant tiles. Relevant tiles are within 200 m from the area of interest
Rscript $basescript/2_choose_relevant_tiles.R $input $nr_strips

# ATMOSPHERIC CORRECTION
cd $base/atcor_6s/PythonWrapper
	
for (( strip_id=1; strip_id<=$nr_strips; strip_id++ ));
do
	python $base/atcor_6s/PythonWrapper/pythonWrV8.py --band4n8alldir=''$baseDG'/0_categ/'$input'/'$input'_P00'$strip_id'_MUL' --outputdir=''$baseDG'/1_atcor_6s/'$input'/'
done

Rscript $basescript/4_hdf2tif.R $input $nr_strips

# Copy files and
# intersect the image footprint and the study area
# Apply orthorectification if necessary 
# use GDAL
Rscript $basescript/5_subset_and_ortho.R $input $nr_strips

# BLOB DETECTION
Rscript $basescript/6_blob_detector_v8.R $input $nr_strips

##### THE MASTER IMAGE FOR EACH STUDY AREA IS PROCESSSED HERE
# !!! Move out of this script !!!
# export oid=1040
# Rscript $basescript/7_tree_size_shape_tiles_v5.R $oid

# nohup Rscript $basescript/7_tree_size_shape_tiles_v5.R $oid &>/dev/null &

# Manual: process the four large Babangona products
# nohup Rscript $basescript/7_tree_size_shape_tiles_v5_Babangona_manual.R &>/dev/null &


# IMAGE MATCHING
Rscript $basescript/8_matching_v5.R $input $nr_strips
Rscript $basescript/8_matching_pan_v5.R $input $nr_strips

# CLOUD MASKING
# Rscript $basescript/11_cloud_masking_v3.r $input

# TREE MASKING
# Rscript $basescript/12_tree_masking.R $input

# EXTRACT POLYGON-BASED SPECTRAL STATISTICS FOR THE CSSL
# Rscript $basescript/13_cssl_statistics_spectral_v7.R $input

# EXTRACT POLYGON-BASED TEXTURE STATISTICS FOR THE CSSL
# Rscript $basescript/13_cssl_statistics_textural_v4.R $input

# # # #################################
# # # #report  processing time
# # # # proc.time()

# # # #ls -LR | grep .tif$ >> list_tif.txt # list all tif files
# # # #find . -maxdepth 1 -type d >> subdirectories.txt # list subdirectories
# # # # find . -maxdepth 3 -iname "metadata*" > 0_metadata.txt
# # # # find . -iname "*.tif" | cut -d'/' -f 3-4 > list_tif.txt
# # # #ls -l --block-size=MB -LR | grep .tif$ >> tif_size.log # size of .tif files
# # # #du -h --max-depth=1 | sort -hr >> dir_size.txt # size of directories
