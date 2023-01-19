#!/bin/bash

# Set current working directory and move to
# directory with output files generated when
# sampling from the prior
curr_dir=$( pwd )
cd ../posterior 
home_dir=$( pwd )

# Get maximum number of chains ran
# and number of alignments
num_alns=$1
num_chains=$2
num_chains=$( echo $num_chains | sed 's/\,/\ /g' )

for dat in $num_alns
do
mkdir -p $dat/mcmc_files_GBM
for i in $num_chains
do

if [[ ! -f $dat/GBM/$i/mcmc.txt ]]
then 
printf "Sorry, no samples for run"$i" and data labelled as "$dat" under GBM ...\n"
printf "Data_"$dat"\tGBM\trun"$i"\n" >> "Not_collected_samples.tsv"
else 
printf "Parsing dat "$dat" for run"$i" under GBM... ... \n"
end=$( wc -l $dat/GBM/$i/mcmc.txt | sed 's/ ..*//' )
if [[ $i -eq 1 ]]
then
begin=1
else 
begin=2 
fi
sed -n ''${begin}','${end}'p' $dat/GBM/$i/mcmc.txt >> $dat/mcmc_files_GBM/mcmc.txt
#sed -n '1,'${end}'p' $dat/$i/GBM/mcmc.txt >> $dat/GBM/$i/mcmc_clean.txt 
fi

done 
done


# NOTE:
# After this script, copy a dummy aln, ctl file, and tree file 
# to generate the FigTree file using the option -1 !