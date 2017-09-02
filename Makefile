SHELL = /bin/sh

export baseDG := $(addprefix $(HOME),/stars/derived/DG_v8)
export basescript :=$(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))
export R_LIBS := $(addprefix $(HOME),/stars/rlibs)
# extend the session's library-path.
export LD_LIBRARY_PATH:=$(addsuffix $(LD_LIBRARY_PATH),:/usr/local/lib)

folders := $(wildcard ~/stars/derived/DG_v8/0_categ/*)
items := $(notdir $(folders))
MULfolders := $(addprefix /,$(items))
MULfolders := $(addsuffix _P001_MUL,$(MULfolders))
MULfolders := $(join $(folders),$(MULfolders))
txtitems := $(addprefix /metadata_,$(items))
txtitems := $(addsuffix .txt,$(txtitems))
txtitems := $(join $(MULfolders),$(txtitems))

# only TAB as a delimiter in IFS, to influence the parsing of metadata.txt files by the /bin/sh read function
export IFS := '	'

all: $(txtitems)

-include generated.mk

generated.mk: Makefile
	# Rscript $basescript/1_metadata.R $input # derive metadata
	@for f in $(items); do \
		echo ${baseDG}/0_categ/$${f}/$${f}_P001_MUL/metadata_$${f}.txt: ${baseDG}/0_categ/$${f}; \
		echo '\t'"@echo Rscript ${basescript}/0_copy_data.R $${f}"; \
		echo ;\
	done > $@
	# cd $baseDG/1_atcor_6s_source/PythonWrapper
	# python $baseDG/1_atcor_6s_source/PythonWrapper/pythonWrV8.py --band4n8alldir='/home/stratouliasd/stars/derived/DG_v8/0_categ/'$input'/'$input'_P001_MUL' --outputdir='/home/stratouliasd/stars/derived/DG_v8/1_atcor_6s/'$input'/'
	@for f in $(items); do \
		echo ${baseDG}/1_atcor_6s/$${f}: ${baseDG}/0_categ/$${f}/$${f}_P001_MUL/metadata_$${f}.txt; \
		echo '\t'"@echo python ${baseDG}/1_atcor_6s_source/PythonWrapper/pythonWrV8.py --band4n8alldir='${baseDG}/0_categ/$${f}/$${f}_P001_MUL' --outputdir='${baseDG}/1_atcor_6s/$${f}/'"; \
		echo ;\
	done >> $@
	# Rscript $basescript/3_export_tif.R $input # convert to .tif
	@for f in $(items); do \
		echo ${baseDG}/1_atcor_6s/$${f}/$${f}_P001_M2AS_R1C1.tif: ${baseDG}/1_atcor_6s/$${f}; \
		echo '\t'"@echo Rscript ${basescript}/3_export_tif.R $${f}"; \
		echo ;\
	done >> $@
	# Rscript $basescript/4_mosaic.R $input # stitch the tiles of the MUL delivery together
	@for f in $(items); do \
		echo ${baseDG}/2_mosaic_r/$${f}/$${f}_P001_M2AS.tif: ${baseDG}/1_atcor_6s/$${f}/$${f}_P001_M2AS_R1C1.tif; \
		echo '\t'"@echo Rscript ${basescript}/4_mosaic.R $${f}"; \
		echo ;\
	done >> $@
	# Rscript $basescript/5_rename.R $input # prepare and rename files
	@for f in $(items); do \
		echo ${baseDG}/3_orthorectification_otb/$${f}/$${f}_P001_MUL.bat: ${baseDG}/2_mosaic_r/$${f}/$${f}_P001_M2AS.tif; \
		echo '\t'"@echo Rscript ${basescript}/5_rename.R $${f}"; \
		echo ;\
	done >> $@
	# /usr/local/bin/otbcli_OrthoRectification -io.in ""$baseDG"/2_mosaic_r/"$input"/"$input"_P001_M2AS.tif?&skipcarto=true" -io.out ""$baseDG"/3_orthorectification_otb/"$input"/"$input"_otb_skipcarto_DEM_nn.TIF" -elev.dem /home/stratouliasd/stars/derived/DEM/DEM_SRTM -interpolator nn
	@for f in $(items); do \
		echo ${baseDG}/3_orthorectification_otb/$${f}/$${f}_otb_skipcarto_DEM_nn.TIF: ${baseDG}/3_orthorectification_otb/$${f}/$${f}_P001_MUL.bat; \
		echo '\t'"@echo /usr/local/bin/otbcli_OrthoRectification -io.in "${baseDG}/2_mosaic_r/$${f}/$${f}_P001_M2AS.tif\?\&skipcarto=true" -io.out "${baseDG}/3_orthorectification_otb/$${f}/$${f}_otb_skipcarto_DEM_nn.TIF" -elev.dem $(HOME)/stars/derived/DEM/DEM_SRTM -interpolator nn"; \
		echo ;\
	done >> $@
	# Rscript $basescript/7_categorize_area.R $input # categorize according to area of study, this step can be skipped if we integrate the study area abbrev. in the filename
	@for f in $(items); do \
		file=${baseDG}/0_categ/$${f}/$${f}_P001_MUL/metadata_$${f}.txt; \
		sed "s/\"//g" $${file} | sed "s/ /_/g" | sed "s/off-nadir/off_nadir/g" | { while read name value; do export $${name}="$${value}"; done; \
		echo ${baseDG}/4_categ/$${Study_area}/$${f}_$${Study_area}_$${Acquisition_date}_$${Satellite}_$${Processing_level}.TIF: ${baseDG}/3_orthorectification_otb/$${f}/$${f}_otb_skipcarto_DEM_nn.TIF; \
		echo '\t'"@echo Rscript ${basescript}/7_categorize_area.R $${f}"; \
		echo ;\
		} \
	done >> $@
	# Rscript $basescript/8_blob_detector_server_clean_v7.R $input # blob detection
	@for f in $(items); do \
		file=${baseDG}/0_categ/$${f}/$${f}_P001_MUL/metadata_$${f}.txt; \
		sed "s/\"//g" $${file} | sed "s/ /_/g" | sed "s/off-nadir/off_nadir/g" | { while read name value; do export $${name}="$${value}"; done; \
		echo ${baseDG}/5_blob_detection/$${f}/all_blobs_image_$${f}.RData: ${baseDG}/4_categ/$${Study_area}/$${f}_$${Study_area}_$${Acquisition_date}_$${Satellite}_$${Processing_level}.TIF; \
		echo '\t'"@echo Rscript ${basescript}/8_blob_detector_server_clean_v7.R $${f}"; \
		echo ;\
		} \
	done >> $@
	# Rscript $basescript/9_blob_combine_server_v3.R $input # stitch tiles (integrated with step_8 at version 8)
	# at this stage the master image for each site is processed at 5_tree_measurement
	# Rscript $basescript/10_matching.R $input # image matching
	@for f in $(items); do \
		file=${baseDG}/0_categ/$${f}/$${f}_P001_MUL/metadata_$${f}.txt; \
		sed "s/\"//g" $${file} | sed "s/ /_/g" | sed "s/off-nadir/off_nadir/g" | { while read name value; do export $${name}="$${value}"; done; \
		if [ $${Study_area} != "NG_Kofa" -a $${Study_area} != "ML_Sukumba" ]; then continue; fi; \
		echo ${baseDG}/6_image_matching/$${f}/$${f}_$${Study_area}_$${Acquisition_date}_$${Satellite}_$${Processing_level}_transformed_topocorr.tif: ${baseDG}/5_blob_detection/$${f}/all_blobs_image_$${f}.RData; \
		echo '\t'"@echo Rscript ${basescript}/10_matching.R $${f}"; \
		echo ;\
		} \
	done >> $@
	# Rscript $basescript/11_cloud_masking.R $input # prepare and rename files
	@for f in $(items); do \
		file=${baseDG}/0_categ/$${f}/$${f}_P001_MUL/metadata_$${f}.txt; \
		sed "s/\"//g" $${file} | sed "s/ /_/g" | sed "s/off-nadir/off_nadir/g" | { while read name value; do export $${name}="$${value}"; done; \
		if [ $${Study_area} != "NG_Kofa" -a $${Study_area} != "ML_Sukumba" ]; then continue; fi; \
		echo ${baseDG}/7_cloud_mask/$${f}/$${f}_cloud_mask.tif: ${baseDG}/6_image_matching/$${f}/$${f}_$${Study_area}_$${Acquisition_date}_$${Satellite}_$${Processing_level}_transformed_topocorr.tif; \
		echo '\t'"@echo Rscript ${basescript}/11_cloud_masking.R $${f}"; \
		echo ;\
		} \
	done >> $@
	# Rscript $basescript/12_tree_masking.R $input # prepare and rename files
	@for f in $(items); do \
		file=${baseDG}/0_categ/$${f}/$${f}_P001_MUL/metadata_$${f}.txt; \
		sed "s/\"//g" $${file} | sed "s/ /_/g" | sed "s/off-nadir/off_nadir/g" | { while read name value; do export $${name}="$${value}"; done; \
		if [ $${Study_area} != "NG_Kofa" -a $${Study_area} != "ML_Sukumba" ]; then continue; fi; \
		echo ${baseDG}/8_tree_mask/$${f}/Tree_mask_$${f}.tif: ${baseDG}/7_cloud_mask/$${f}/$${f}_cloud_mask.tif; \
		echo '\t'"@echo Rscript ${basescript}/12_tree_masking.R $${f}"; \
		echo ;\
		} \
	done >> $@
	# Rscript $basescript/13_cssl_statistics.R $input # extract polygon-based statistics for CSSL
	@for f in $(items); do \
		file=${baseDG}/0_categ/$${f}/$${f}_P001_MUL/metadata_$${f}.txt; \
		sed "s/\"//g" $${file} | sed "s/ /_/g" | sed "s/off-nadir/off_nadir/g" | { while read name value; do export $${name}="$${value}"; done; \
		if [ $${Study_area} != "NG_Kofa" -a $${Study_area} != "ML_Sukumba" ]; then continue; fi; \
		export itemid=$$(echo $${f} | cut -c1-12); \
		echo ${baseDG}/9_cssl_statistics/$${f}/$${itemid}_basicstats.csv: ${baseDG}/6_image_matching/$${f}/$${f}_$${Study_area}_$${Acquisition_date}_$${Satellite}_$${Processing_level}_transformed_topocorr.tif; \
		echo '\t'"@echo Rscript ${basescript}/13_cssl_statistics.R $${f}"; \
		echo ;\
		} \
	done >> $@
	
#IFS := $(OLD_IFS)

