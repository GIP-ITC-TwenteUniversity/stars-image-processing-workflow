#!/bin/bash
# cd /home/tolpekin/stars/derived/scripts/v10
# nohup ./looper_sequential.sh &>/dev/null &
export stars=~/stars
export base=~/stars/derived
export baseDG=~/stars/derived/DG_v10
export basescript=~/stars/derived/scripts/v10
export R_LIBS=~/stars/rlibs

cd $stars/acquired/DG
#Rscript $basescript/0_multi_specs.R # export .csv with multi-image specifications 

# IN BASH

#DIRS=(`find . -maxdepth 1 -type d -name '*_01' | cut -c 3- `)
DIRS=( $(find . -maxdepth 1 -type d -name '*_01' | cut -c 3-) )
echo ${#DIRS[@]}

# get length of an array
arraylength=${#DIRS[@]}
echo $arraylength

for (( j=0; j<${arraylength}; j++ ));
do
  echo "Processing ${DIRS[j]} file...";
  $basescript/shell_command.sh ${DIRS[j]};
done
