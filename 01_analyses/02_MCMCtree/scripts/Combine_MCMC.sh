#!/bin/bash

curr_dir=$( pwd )
dat=$1
dirname=$2
seqchains=$3
clock=$4

mkdir -p $dirname
count=0
for i in $seqchains
do
	count=$(( count + 1 ))	
	if [[ ! -f $dat/$clock/$i/mcmc.txt ]]
	then 
		printf "Sorry, no samples for run"$i" under "$clock" ...\n"
		printf $clock"\trun"$i"\n" >> "Not_collected_samples.tsv"
	else 
		printf "Parsing dat for run"$i" under "$clock" ... ... \n"
		end=$( wc -l $dat/$clock/$i/mcmc.txt | sed 's/ ..*//' )
	if [[ $count -eq 1 ]]
		then
			begin=1
		else 
			begin=2 
		fi
		sed -n ''${begin}','${end}'p' $dat/$clock/$i/mcmc.txt >> $dirname/mcmc.txt
	fi
		
done  

# NOTE:
# After this script, you need to copy dummy aln, ctl file, and tree file 
# to generate the FigTree file using the option -1 !