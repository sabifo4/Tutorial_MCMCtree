#!/bin/bash

# Get args
aln=$1         # 1, 2, etc.
pipedir=$2     # Path to pipeline dir
name_wd=$3     # Name of the working directory, e.g., `euk110`
baseml_exec=$4 # Name of the BASEML binary
bool_paml=$5   # Boolean, TRUE if PAML is in PATH, FALSE otherwise
# Replace vars in template bash script for job array
cp pipeline_Hessian_BASEML_template_PC.sh $pipedir/pipeline_Hessian.sh
if [[ $aln -eq 1 ]]
then 
sed -i 's/.seq\ 1\ NUM..*/1/' $pipedir/pipeline_Hessian.sh
else 
sed -i 's/NUM/'${aln}'/' $pipedir/pipeline_Hessian.sh
fi
# Replace name of working directory
upd_wd=$( echo $name_wd | sed 's/\//\\\//g' | sed 's/_/\\_/g' )
sed -i 's/WDNAME/'${upd_wd}'/g' $pipedir/pipeline_Hessian.sh
# Add name for BASEML binary
sed -i 's/BASEMLEXEC/'${baseml_exec}'/g' $pipedir/pipeline_Hessian.sh
# Replace relative path to main dir by just the executable name
if [[ $bool_paml =~ "Y" ]]
then
sed -i 's/\$main_dir\/'${baseml_exec}'/'${baseml_exec}'/' $pipedir/pipeline_Hessian.sh
fi