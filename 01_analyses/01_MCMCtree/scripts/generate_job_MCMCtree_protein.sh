#!/bin/bash

# Get args
dir=$1     # Alignment #1, #2, #3... The alignment being parsed at the moment
clock=$2   # `GBM` or `ILN`
ndat=$3    # 1, 2, 3... As many blocks as partitions in the alignment
pipeloc=$4 # Path to MCMCtree pipeline dir
runmcmc=$5 # Command to execute MCMCtree
nchains=$6 # Number of MCMCs to run

# Replace vars in template bash script for job array
cp pipeline_MCMCtree_protein_template.sh $pipeloc/$clock/pipeline_$clock".sh"
sed -i 's/DIR/'${dir}'/g' $pipeloc/$clock/pipeline_$clock".sh"
sed -i 's/CLK/'${clock}'/g' $pipeloc/$clock/pipeline_$clock".sh"
sed -i 's/NUMPARTS/'${ndat}'/' $pipeloc/$clock/pipeline_$clock".sh"
sed -i 's/CMDRUN/'${runmcmc}'/' $pipeloc/$clock/pipeline_$clock".sh"

if [[ $nchains -eq 1 ]]
then 
sed -i 's/for\ TASK\_ID\ in..*//' $pipeloc/$clock/pipeline_$clock".sh"
sed -i 's/^done//' $pipeloc/$clock/pipeline_$clock".sh"
sed -i 's/^do//' $pipeloc/$clock/pipeline_$clock".sh"
else 
sed -i 's/NUM/'${nchains}'/' $pipeloc/$clock/pipeline_$clock".sh"
fi