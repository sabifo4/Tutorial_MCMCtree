#!/bin/bash
#$ -cwd                    # Run the code from the current directory
#$ -V                      # Export environment to job
#$ -j y                    # Merge the standard output and standard error
#$ -l h_rt=240:00:00       # Limit each task to 10 days
#$ -l h_vmem=REQRAM        # Requested RAM
#$ -t 1-NUM                

#==========================================================================================#
# Contact Sandra Alvarez-Carretero for any doubts about this script: sandra.ac93@gmail.com #
#==========================================================================================#

# ------------------------------------- #
# Creating file structure to run BASEML #
# ------------------------------------- # 

# 1. Find global dirs for paths
pipeline_dir=$( pwd )
main_dir=$( echo $pipeline_dir | sed 's/WDNAME..*/WDNAME/' )
cd $main_dir/Hessian/$SGE_TASK_ID
home_dir=$( pwd ) 

# 3. Create specific log file
exec 3>&1> >(while read line; do echo "$line" >> $pipeline_dir/log.hessian.dir$SGE_TASK_ID.txt; done;) 2>&1
start=`date`
echo Job starts":" $start

# 4. Start analysis
echo The analyses will take place in directory $home_dir
printf "\n"
# Move to analysis dir
cd $home_dir
# Soft link the tmp* files here 
ln -s $home_dir/prepare_baseml/tmp0001.ctl .
ln -s $home_dir/prepare_baseml/tmp0001.trees .
ln -s $home_dir/prepare_baseml/tmp0001.txt .

# 5. Run BASEML
printf "\nRunning BASEML to calculate the Hessian and the gradient...\n"
$main_dir/BASEMLEXEC tmp0001.ctl

# 6. Close
printf "\n"
echo BASEML FINISHED"!"
printf "\n"
end=`date`
echo Job ends":" $end
