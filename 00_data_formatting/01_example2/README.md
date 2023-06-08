# Data filtering

For this example, we assume that the gene alignments and the corresponding tree files have been inferred following the best practices and making sure the best method has been chosen for such inference. Nevertheless, additional quality control checks can be carried out to further filter the given alignment and tree files. Therefore, before proceeding with timetree inference, we will make sure that:

1. We check for outliers in our dataset.
2. Our parsed files are converted to the correct format.

Now... Let's get started with our data filtering analyses!

## General overview of filtering steps

In this tutorial, we will use the following example data:

* [`Gene alignments`](00_raw_data/aln): inside this directory, you will find one file per gene alignment in FASTA format.
* [`Gene trees`](00_raw_data/trees): inside this directory, you will find one file per gene tree with branch lengths in Newick format.

Note that both files have the same name (except for the extension!) so it is easier to match the gene alignments with the corresponding gene trees. To filter our dataset, we will follow the next steps:

1. Run a relative branch length test to check whether there are any gene trees that should be discarded due to their length being larger than 60% of the total tree length.
2. Classify genes according to the corresponding evolutionary rate (i.e., rank from slow- to fast-evolving) and order them from slow- to fast-evolving.
3. Concatenate alignments from slow- to fast-evolving. Then, keep the unpartitioned datasets but additionally generate partitioned datasets (e.g., CPs) with four blocks.

Below, you can find one section for each step with the information needed to reproduce our results.

### Step 1: relative branch length test

According to [Springer & Gatesy, 2017](https://www.tandfonline.com/doi/full/10.1080/14772000.2017.1401016) and [dos Reis et al., 2012](https://pubmed.ncbi.nlm.nih.gov/22628470/), the relative branch length test can be useful to detect either misaligned or misidentified orthologs in gene alignments. We carried out this test following the same methodology described in [Álvarez-Carretero et al., 2021](https://t.co/W0RtKVHAjv).

I wrote the [`Filter_genes.R`](scripts/Filter_genes.R) to format, parse, and analyse this dataset. The first commands you will find have to do with data preparation and formatting. The commands used to run the relative branch length test can be found after the comment "Start search for outliers running the relative branch length test" (at the moment of writing, from lines 58-170).

Despite no genes failed the test, we found two potential outliers after our visual checks (i.e., look at the PDF files that have been output in the newly created `out_RData` directory inside [`00_raw_data`](00_raw_data)). Nevertheless, we will keep these genes for a bit longer until the next step.

## Step 2. Calculate the evolutionary rate

Now, we can continue with the next sections that are part of the script [`Filter_genes.R`](scripts/Filter_genes.R) so that we can estimate the mean evolutionary rate of the genes previously estimated to order the genes from slow- to fast-evolving. We will then generate a new directory to save the ordered genes according to this criterion.

In the previous step, we calculated the tree height of each gene tree (in substitutions/site per time unit; 100 Myr) and used the mean of the minimum and maximum ages used to calibrate the root of the tree as the mean root age (i.e., the time of divergence at the root of the phylogeny). To estimate the rate of each gene tree, we divided the gene tree height into the mean root age. You will find this calculation in lines 70 and 72 (at the time of writing) in the script [`Filter_genes.R`](scripts/Filter_genes.R).

When you run the commands from line 172 (at the time of writing), you will rank the genes from slow- to fast-evolving. After that, each gene alignment will be matched with the corresponding gene tree, and a copy of both files will be saved in a separate directory where genes will be ordered according to the estimated evolutionary rate. To ease subsequent analyses using job arrays in an HPC, each directory will be labelled from `1` to `n`, where `n` is the maximum number of gene alignments included for that specific dataset. The file structure looks as follows:

```txt
01_example2
  |- 00_raw_data
  |   |-aln
  |   |-out_logs
  |   |-out_RData
  |   |-trees
  |
  |- 01_ordered_genes
      |- [0-9]*
           |-[0-9].*.[fasta|tree]
```

You will find the `log_R_copy_ordered_genes_s2f.txt` log file saved in the newly created `out_logs` directory inside [`00_raw_data`](00_raw_data), which lists the gene trees from slow- to fast-evolving as well as the directory in which they have been now saved together with the corresponding gene alignment.

## Step 3. Prepare files for PAML

Now, we are ready to format the gene alignment files prior to run the [`fasta-phylip-partitions`](https://github.com/sabifo4/fasta-phylip-partitions) pipeline to generate the final alignments that we will use to run `MCMCtree`:

```sh
# Run from `01_example2`
# This will prepare the FASTA files for the alignments
chmod 775 ../../src/*sh
# This bash script takes three arguments: 
# $1: the name of the 
#     directory where the directories from `1` to `n` are
#     saved for each type of tree. In this case, we have 
#     the directory `01_ordered_genes`.
# $2: the number of genes to be used per partition
#     In this example, 13/4 = 3.25 ~ 3 genes per partition 
# $3: The ABSOLUTE path of the current working directory.
#
# NOTE: The directory you pass to otpion 1 can ONLY have 
# the directories from `1` to `n`, no other files!
curr_dir=$( pwd )
../../src/create_filestruct_partitions.sh 01_ordered_genes 3 $curr_dir
```

Once you have run the previous script, you will see the following file structure inside [`02_alignments`](02_alignments):

```text
02_alignments/
   |- conc/ 
   |   |- geneX_conc_<gene_name>.fasta
   |
   |- partY/  # A total of "Y" directories, in this case 4 as there are 4 partitions 
       |- geneX_partY_<gene_name>.fasta
```

In essence, the previous script has formatted the FASTA files and renamed them in such a way that the pipeline mentioned above will be able to concatenate all the genes within each subdirectory to generate the corresponding alignments:

* `conc`: all the genes here will be concatenated in a unique alignment ordered from slow- to fast-evolving.
* `partY`: for each partition, an alignment will be generated with the corresponding genes ordered from slow- to fast-evolving. In that way, there will be one alignment for all the genes under directory `part1`, another for those under directory `part2`, etc. Then, we will concatenate them in a unique alignment file, which will be the partitioned alignment file.

Now, we can use the [`fasta-phylip-partitions`](https://github.com/sabifo4/fasta-phylip-partitions) pipeline to generate the final "main" alignments with the following code snippet. To make things easier, we have included this pipeline in [the `src` directory in this repository too](../../src/fasta-phylip-partitions) without examples to save space:

```sh
# Run from `01_example2` and set curr dir
home_dir=$( pwd )
cd 02_alignments
curr_dir=$( pwd )
# Get the `species_names.txt` that I have already saved in `00_raw_data` for this example.
# This file lists the total number of species for which data have been collected.
# Some alignments can have all or some taxa or some of those, so we use this file
# to find out which taxa have missing sequences in the alignment.
# The `species_names.txt` file is very important for the pipeline to work, so 
# do NOT change this file name! We will also give permissions to the tool!
cp $home_dir/00_raw_data/species_names.txt $curr_dir/
chmod 775 ../../../src/fasta-phylip-partitions/src/Run_tasks.sh
chmod 775 ../../../src/fasta-phylip-partitions/src/Tools/*
# Loop over each directory to obtain individual concatenated
# alignments for each specific partition
for i in */    # `conc`, `part1`, `part2`, `part3`, `part4`
do 
name_dir=$( echo $i | sed 's/\///' )
printf "\n\n[[ PARSING GENES IN DIR "$name_dir" ]]\n\n"
# Now, run the pipeline. In essence, the first argument is the current 
# directory ("."), the second argument the tag ID for the job
# ("ex2_$name_dir"), and then a tag specifying if the alignment 
# needs to be partitioned into CPs, which we do want, and hence use 
# "partY".
# If you do not have this pipeline on your system, you have all the 
# informatoin you need and an extensive tutorial with more details about 
# how to use it [here](https://github.com/sabifo4/fasta-phylip-partitions/blob/main/README.md)
# NOTE: If you are running this code in a directory that you have synched to Dropbox or another 
# cloud storage and you have issues, just move the folder out of the synched directory and run the 
# command below again.
cd $curr_dir/$name_dir 
cp $curr_dir/species_names.txt .
../../../../src/fasta-phylip-partitions/src/Run_tasks.sh . ex2_$name_dir partY
done
# Move back to $curr_dir and remove copy of `species_names.txt` file
cd $curr_dir
rm species_names.txt
```

Once the code above has run to generate the alignments (both concatenated and partitioned alignments), we just need to concatenate the four partitions into a unique file. We will do this for the alignments without codon partitioning scheme applied and then with those alignments in which only the 1st and the 2nd codon positions have been kept:

```sh
# Go to `02_alignments` directory and run the 
# following code
home_dir=$( pwd )
cd $home_dir
# Go first with `conc` alignment
printf "\n\n[[ PARSING DIR \"conc\" ]]\n\n"
# Copy `conc` alignment in main dir
cp conc/phylip_format/02_concatenated_alignments/*aln $home_dir/ex2_conc.aln 
cp conc/phylip_format/02_concatenated_alignments/part12/*aln $home_dir/ex2_12CP_conc.aln 
# Now, generate alignment files with 4 blocks for 
# slow- to fast-evolving partitioning scheme
for i in `seq 1 4`
do
printf "\n\n[[ PARSING DIR \"part"$i"\" ]]\n\n"
cat part$i/phylip_format/02_concatenated_alignments/*aln >> $home_dir/ex2_4parts.aln
cp part$i/phylip_format/02_concatenated_alignments/*aln $home_dir/ex2_part$i.aln
cat part$i/phylip_format/02_concatenated_alignments/part12/*aln >> $home_dir/ex2_12CP_4parts.aln
cp part$i/phylip_format/02_concatenated_alignments/part12/*aln $home_dir/ex2_12CP_part$i.aln
done
cd $home_dir
```

Now, we will have the following file structure for each data subset directory:

```text
02_alignments/
  |- conc/
  |- part1/ 
  |- part2/ 
  |- part3/ 
  |- part4/ 
  |- ex2*aln

```

There are as many alignment files (i.e., files which name ends with `aln` in the file structure described above) as partitioning schemes have been devised:

* **Concatenation**: all genes concatenated in one block in the alignment file (i.e., `*conc.aln`).
* **Partitioning according to codon position**: we have an alignment file with one block with 1st+2nd CPs of all genes concatenated ( i.e., `*_12CP*.aln`). Another with one block with only 3rd CPs is available, but we are not using it for this analysis.
* **Partitioning according to evolutionary rate**: we have alignment files with four blocks in which genes are ordered from slow- to fast-evolving. These files have all the corresponding genes of the data subset concatenated and then partitioned following this scheme (i.e., `*_4parts*.aln`). In addition, we have also copied the individual partitions so the gradient and the Hessian can be calculated for each partition individually (i.e., `*_partX.aln`, where `X` is either 1, 2, 3, or 4, according to the partition number).

We can now count the missing taxa for each data subset with the script [count_missingdat.pl](../../src/count_missingdat.pl). We will run it using the next commands:

```sh
# Run from `01_example2`
home_dir=$( pwd )
cd 02_alignments
mkdir out_count_NA
for i in *part[0-9].aln *conc.aln
do
printf "Parsing alignment "$i"... ...\n"
printf "Parsing alignment "$i"... ...\n" > log_count_missdat_$i".txt"
perl ../../../src/count_missingdat.pl $i >> log_count_missdat_$i".txt"
mv log_count_missdat_$i".txt" out_count_NA
done
# Get a summary!
for i in *part[0-9].aln *conc.aln
do
name=$( echo $i | sed 's/\.aln//' )
printf "<< DIR "$name" >>\n" >> out_count_NA/$name"_countNA.tsv"
sed -n '1,2p' out_count_NA/$name"_avgmissdata.txt" >> out_count_NA/$name"_countNA.tsv"
sed -n '7,7p' out_count_NA/$name"_avgmissdata.txt" >> out_count_NA/$name"_countNA.tsv"
sed -n '9,9p' out_count_NA/$name"_avgmissdata.txt" >> out_count_NA/$name"_countNA.tsv"
printf "\n" >> out_count_NA/$name"_countNA.tsv"
done
cd $home_dir
```

Now that we have all the alignments generated for each partitioning scheme and for each data subset and we have counted the corresponding missing data, we can start our timetree inference!
