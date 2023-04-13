#!/bin/bash

# Get args
dir=$1 # Alignment #1, #2, #3... The alignment being parsed at the moment

# Find global dirs for paths
scripts_dir=$( pwd )
cd .. 
home_dir=$( pwd )
cd $home_dir/trees/uncalibrated
tree_dir=$( pwd )
tree_name=`ls *tree | sed 's/..*\///'`
cd $home_dir/Hessian/$dir/prepare_codeml
codeml_dir=$( pwd ) 
cd $home_dir/control_files
ctl_dir=$( pwd )
cd $home_dir/alignments/$dir
aln_dir=$( pwd )
aln_name=`ls *phy | sed 's/..*\///'`

# Move to dir where CODEML will run
cd $codeml_dir 
curr_dir=$( pwd )

# Get path to alignments
printf "You are generating the CODEML files in the following directory:\n"
printf "Dir = "$codeml_dir" ... ...\n"
# Prepare paths to replace in ctl file
name_aln=$( echo $aln_dir/$aln_name | sed 's/\//\\\//g' | sed 's/\_/\\\_/g' | sed 's/\./\\\./g' )
name_tree=$( echo $tree_dir/$tree_name | sed 's/\//\\\//g' | sed 's/\_/\\\_/g' | sed 's/\./\\\./g' )
name_lgdat=$( echo $ctl_dir/lg.dat | sed 's/\//\\\//g' | sed 's/\_/\\\_/g' | sed 's/\./\\\./g' )
cp $ctl_dir/prepcodeml.ctl $codeml_dir
printf "\n" 
sed -i 's/ALN/'${name_aln}'/' $codeml_dir/*ctl
sed -i 's/TREE/'${name_tree}'/' $codeml_dir/*ctl
sed -i 's/lg\.dat/'${name_lgdat}'/' $codeml_dir/*ctl
