#!/bin/bash

# Get args
dir=$1     # Alignment #1, #2, #3... The alignment being parsed at the moment
pipeloc=$2 # Pipeline directory

# Replace vars in template bash script for job array
cp pipeline_Hessian_template.sh $pipeloc/pipeline_BASEML.sh
if [[ $dir -eq 1 ]]
then 
sed -i 's/for\ TASK\_ID\ in..*//' $pipeloc/pipeline_BASEML.sh
sed -i 's/^done//' $pipeloc/pipeline_BASEML.sh
sed -i 's/^do//' $pipeloc/pipeline_BASEML.sh
sed -i 's/\$TASK\_ID/1/' $pipeloc/pipeline_BASEML.sh
else 
sed -i 's/NUM/'${dir}'/' $pipeloc/pipeline_BASEML.sh
fi
