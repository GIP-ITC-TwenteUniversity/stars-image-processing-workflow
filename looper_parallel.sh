#!/bin/bash
# nohup ./looper_parallel.sh &>/dev/null &
#cd /home/stratouliasd/stars/derived/scripts/v9
export base=~/stars/derived
export baseDG=~/stars/derived/DG_v9
export basescript=~/stars/derived/scripts/v9
export R_LIBS=~/stars/rlibs

cd $baseDG/0_categ
#Rscript $basescript/0_multi_specs.R

# DIRS=$(find . -maxdepth 1 -type d -name '*_01' | cut -c 3-)
# echo $DIRS
#DIRS=( `find . -maxdepth 1 -type d -name '*_01' | cut -c 3- `)
DIRS=( $(find . -maxdepth 1 -type d -name '*_01' | cut -c 3-) )
echo ${DIRS[142]}

for j in $(seq 0 144) # process 145 images
do
	#( baseDG="/home/stratouliasd/stars/derived/DG_v9"; basescript="/home/stratouliasd/stars/derived/scripts/v9"; cd $baseDG/0_categ; DIRS=( `find . -maxdepth 1 -type d -name '*_01' | cut -c 3- `); $basescript/shell_command.sh ${DIRS[j]}; ) &
	( cd $baseDG/0_categ; DIRS=( `find . -maxdepth 1 -type d -name '*_01' | cut -c 3- `); $basescript/shell_command.sh ${DIRS[j]}; ) &
	if (( $j % 5 == 0 )); then wait; fi
done
wait

