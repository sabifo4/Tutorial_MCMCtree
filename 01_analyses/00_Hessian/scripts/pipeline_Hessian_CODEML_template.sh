#!/bin/bash

#==========================================================================================#
# Contact Sandra Alvarez-Carretero for any doubts about this script: sandra.ac93@gmail.com #
#==========================================================================================#

# Set a `for` loop to run the tasks with as
# many alignments as defined by the user
for TASK_ID in `seq 1 NUM`
do

# ------------------------------------- #
# Creating file structure to run CODEML #
# ------------------------------------- # 

# 1. Find global dirs for paths
pipeline_dir=$( pwd )
main_dir=$( echo $pipeline_dir | sed 's/\/main\_protein\/..*/\/main\_protein\//' )
cd $main_dir/Hessian/$TASK_ID
home_dir=$( pwd ) 

# 2. Create specific log file
exec 3>&1> >(while read line; do echo "$line" >> $pipeline_dir/log.hessian.dir$TASK_ID.txt; done;) 2>&1
start=`date`
echo Job starts":" $start

# 3. Start analysis
echo The analyses will take place in directory $home_dir
printf "\n"
# Move to analysis dir
cd $home_dir
# Soft link the tmp* files here 
ln -s $home_dir/prepare_codeml/tmp0001.ctl .
ln -s $home_dir/prepare_codeml/tmp0001.trees .
ln -s $home_dir/prepare_codeml/tmp0001.txt .

# 4. Run CODEML
printf "\nRunning CODEML to calculate the Hessian and the gradient...\n"
codeml tmp0001.ctl

# 5. Close
printf "\n"
echo CODEML FINISHED"!"
printf "\n"
end=`date`
echo Job ends":" $end
done
