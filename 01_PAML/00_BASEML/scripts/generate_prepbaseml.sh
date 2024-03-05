#!/bin/bash

# Get args
dir=$1 # Alignment #1, #2, #3... The alignment being parsed at the moment

# Find global dirs for paths
scripts_dir=$( pwd )
cd .. 
home_dir=$( pwd )
cd $home_dir/trees/uncalibrated/$dir
tree_dir=$( pwd )
tree_name=`ls *tree | sed 's/..*\///'`
cd $home_dir/Hessian/$dir/prepare_baseml
baseml_dir=$( pwd ) 
cd $home_dir/control_files
ctl_dir=$( pwd )
cd $home_dir/alignments/$dir
aln_dir=$( pwd )
aln_name=`ls *phy | sed 's/..*\///'`

# Move to dir where BASEML will run
cd $baseml_dir 
curr_dir=$( pwd )

# Get path to alignments
printf "You are generating the BASEML files in the following directory:\n"
printf "Dir = "$baseml_dir" ... ...\n"
# Prepare paths to replace in ctl file
name_aln=$( echo $aln_dir/$aln_name | sed 's/\//\\\//g' | sed 's/\_/\\\_/g' | sed 's/\./\\\./g' )
name_tree=$( echo $tree_dir/$tree_name | sed 's/\//\\\//g' | sed 's/\_/\\\_/g' | sed 's/\./\\\./g' )
cp $ctl_dir/$dir/prepbaseml*ctl $baseml_dir
printf "\n" 
sed -i 's/ALN/'${name_aln}'/' $baseml_dir/*ctl
sed -i 's/TREE/'${name_tree}'/' $baseml_dir/*ctl
