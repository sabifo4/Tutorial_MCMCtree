#!/bin/bash

# =========================================================================== #
# This script goes through all the `*aln` files saved in                      #
# the `01_ordered_genes/$name_dir/` and carries out the following tasks:      #                               
#                                                                             #
#  * Prints a concatenated alignment and then 4 partitioned alignments in,    #
#    4 separate directories.                                                  #
#  * Then adds an extra blank line which will be used by the perl script      #
#    that is subsequently run to indicate the end of one gene alignment.      #
#                                                                             #
# You are expected to run this script from main dir for each data subset,     #
# e.g., `01_example2` as:                                                         #
#                                                                             #
# <path_to_script>/create_filestruct_partitions.sh                                  #
#                                                                             #
# =========================================================================== #
# Any questions/doubts/bugs, please send a message to:                        #
# Sandra Alvarez-Carretero, <sandra.ac93@gmail.com>                           #
# =========================================================================== #

# 1. Set counters (e.g., 13/4 = 3.25 ~ 3 genes per partition)
i=0
jinc=$2 # num_genes per partition
j=$2
count=0
count_per_parts=0
name_dir=$1 # tree

# 2. Get scripts dir and then move back to main dir
curr_dir=$3 # That should be the dir where you run this script from
scripts_dir=$( pwd ) # Where the alignments are
cd $curr_dir

# 3. Create `02_alignments` if not there,
# same with `out_logs`
if [ ! -d 02_alignments ] 
then 
mkdir 02_alignments 
fi 
if [ ! -d out_logs ] 
then 
mkdir out_logs
fi 

# 4. Put the corresponding number of fasta alignments in the
#    corresponding partition so later they can be input files 
#    for the `fasta-phylip-partitions` pipeline
for partition in `seq 1 4`
do 
	# Create folder for each partition 
	mkdir -p 02_alignments/{conc,part$partition}
		
	# If it is not the last partition
	if [ $partition != 4 ]
	then
		while [ $count_per_parts -lt $jinc ]
		do
			i=$(( i + 1 ))
			count=$(( count + 1 ))
			count_per_parts=$(( count_per_parts + 1 ))
			name_gene=$( echo `ls $name_dir/${i}/*fasta` | sed 's/..*\///' | sed 's/\.fasta//' )
			# Copy FASTA files to dirs
			cp $name_dir/${i}/*fasta 02_alignments/conc/gene$count"_conc_"$name_gene.fasta
			cp $name_dir/${i}/*fasta 02_alignments/part$partition/gene$count"_part"$partition"_"$name_gene.fasta
			echo Gene $count visited
			echo Gene $count visited >> out_logs/log_for_block_partition_$name_dir.txt
		done 
		
		# Update counters for genes
		printf "Total of genes for partition "$partition" in data subset: "$count_per_parts"\n"
		printf "Total of genes for partition "$partition" in data subset: "$count_per_parts"\n" >> out_logs/log_for_block_partition_$name_dir.txt
		# Count equals to i as there are no missing dirs !
		i=$( echo $count )
		j=$(( i + jinc ))
		echo Check that i is equal to $i and that count is equal to $count. Is j equal to $(( i + jinc )) and $(( count + jinc ))
		#j=$(( count + jinc ))
		count_per_parts=$( echo 0 )
	
	# Otherwise, just reach until the end of genes available 
	# in the directory
	else
		tot=$( ls $name_dir | wc -l )
		k=$(( tot ))
		# Set i as the last gene, check point in case something
		# went wrong before
		i=$( echo $count )
		while [ $i -lt $k ]
		#while [ $count -lt $k ]
		do
			i=$(( i + 1 ))
			count=$(( count + 1 ))
			name_gene=$( echo echo `ls $name_dir/${i}/*\.fasta` | sed 's/..*\///' | sed 's/\.fasta//' )
			# Copy FASTA files to dirs
			cp $name_dir/${i}/*fasta 02_alignments/conc/gene$count"_conc_"$name_gene.fasta
			cp $name_dir/${i}/*fasta 02_alignments/part$partition/gene$count"_part"$partition"_"$name_gene.fasta
			echo Gene $count visited
			echo Gene $count visited >> out_logs/log_for_block_partition_$name_dir.txt
			count_per_parts=$(( count_per_parts + 1 ))
		done 
		printf "Total of genes for partition "$partition": "$count_per_parts"\n"
		printf "Total of genes for partition "$partition": "$count_per_parts"\n" >> out_logs/log_for_block_partition_$name_dir.txt
	
	fi
	
	printf "Partition "${partition}" finishes !\n"

done 
