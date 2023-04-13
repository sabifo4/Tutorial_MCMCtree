# Timetree inference

>> **IMPORTANT NOTE**: The dummy alignment with protein data is yet to be included in this repository. The instructions below, however, can be followed to run analyses with protein data.

## 1. Setting the file structure to run `MCMCtree`

Given that we already generated the calibrated tree, we have the alignment file (or files, if you are running this tutorial with more than one gene), and we have just generated the `in.BV` file with the gradient and the Hessian... We have everything we need to run `MCMCtree`!

However, we first need to create the corresponding file structure, which will follow the same structure as the directories where the Hessian and the gradient were calculated (see [the previous tutorial to run `CODEML`](../01_Hessian/README_protein.md)). You will need to run the following code snippet:

```sh
# Run from `main_protein` dir on your local
# PC. Please change directories until
# you are there. Then run the following
# commands.
# Modify `num_aln` if you are running
# this tutorial with your own dataset.
# Note that, with large datasets, we
# strongly recommend you run at least
# 32 chains. This dataset is quite small
# so we will just run 5.
num_aln=1
num_chains=5
for i in `seq 1 $num_aln`
do
mkdir -p pipelines_MCMCtree/{GBM,ILN}
for j in `seq 1 $num_chains`
do
mkdir -p MCMCtree/$j/{GBM,ILN}/$i
done
done
```

>**IMPORTANT NOTE**: When sampling from the posterior, the likelihood is being calculated or approximated, depending on the `userdata` option you set in the control file to run `MCMCtree`. In other words, the larger the dataset, the more time it will take for `MCMCtree` to finish. For the dummy dataset of this tutorial, five chains are more than enough (we could have even had enough samples with 2 chains, which is the minimum amount of chains you must run to check for chain convergence). Nevertheless, sometimes you might find that five chains may not be enough once you run the MCMC diagnostics (we will see that later in the tutorial). If that's the case, you may want to increase the number of chains you are to run. However, that's a decision that you need to take based on the following:
>
> * **Types of computational resources you can use**: can you use as many nodes and partitions as you like in your HPCs or are you sharing and need to be cautious of not sending too many jobs so that other people can also run their analyses?
> * **Type of analysis you are running**: are you just exploring your data and want a rough result or have you already done that and want results to add in your manuscripts?
> * **Type of data**: which kind of dataset are you working with? Are there just 3 or 5 genes for 10 taxa or are you working with a huge phylogenomic alignment with 6K taxa and billions of base pairs?
> You might need to find a "sweet spot" for the amount of chains you are to run when sampling from the posterior in case you have computational limitations given your resources. You can always start by running fewer chains and, if you need more samples because most of your chains have not passed filters or seem to not have converged or the ESS is too low, you can always send more jobs to your HPCs to run more chains.

Now, we will use the same approach we used to generate the pipelines to run `CODEML`. In this case, however, the bash scripts have been modified accordingly. You will find them in the [`scripts` directory](scripts), from where you will copy them onto the `main_protein/scripts` directory that you will have already generated.

```sh
# Run from the `01_analyses/02_MCMCtree/scripts`
# dir on your local PC. Please change
# directories until you are there. Then run
# the following commands.
cp *protein*sh ../../../main_protein/scripts/
```

The `main_protein` directory  will now have these two extra directories with the corresponding subdirectories:

```text
main_protein
  |- MCMCtree
  |    |- [1-5]
  |         |- [GBM|ILN]
  |              |- [0-9] # As many as available datasets
  |                
  |- pipelines_MCMCtree
       |- [GBM|ILN]/
```

To prepare the script to run `MCMCtree`, you can run the bash scripts mentioned above using the commands below:

```sh
# Run from the `main_protein` dir on your local
# PC. Please change directories until
# you are there. Then run the following
# commands.
home_dir=$( pwd )
cd scripts
chmod 775 *sh
num_aln=1
for i in `seq 1 $num_aln`
do 
printf "Generating job array for dir "$i" and both clocks ... ...\n\n"
# The first argument corresponds to the number of 
# gene alignments, the second to the clock model, 
# the third to the number of partitions in the 
# alignment, the fourth to the directory where
# the pipeline # will be executed: `pipelines_MCMCtree`,
# the fifth the command to execute MCMCtree
# (i.e., `mcmctree` in this example), and the fifth
# to the number of chains used. Please modify the 
# command below if you are using your own dataset.
num_chains=5
./generate_job_MCMCtree_protein.sh $i GBM 1 $home_dir/pipelines_MCMCtree mcmctree $num_chains
./generate_job_MCMCtree_protein.sh $i ILN 1 $home_dir/pipelines_MCMCtree mcmctree $num_chains
done
```

Now, before running `MCMCtree` to sample from the posterior, we will run `MCMCtree` but sampling from the prior. In this way, we can verify whether the user-specified prior (the probability distributions that we want the software as specified by us in the calibrated nodes) and the effective prior (the probability distribution the software will use). Oftentimes, to deal with truncation issues, the effective prior might differ from the user-specified prior and not all possible ages framed by the distribution specified by the user will be included (see an extensive analysis about this effect in [dos Reis et al. 2015](https://pubmed.ncbi.nlm.nih.gov/26603774/)). To assess whether this occurs in our analysis, we will run first `MCMCtree` without the data (i.e., `MCMCtree` will be sampling from the prior distribution and will only use the information provided by the calibrations included in the user tree).

First, we will generate a directory where `MCMCtree` will run and where all the results when sampling from the prior will be saved:

```sh
# Run from `main_protein` dir on your local
# PC. Please change directories until
# you are there. Then run the following
# commands.
num_aln=1
num_chains=5
for i in `seq 1 $num_aln`
do
for j in `seq 1 $num_chains`
do
mkdir -p MCMCtree_prior/$j/CLK/$i
done
done
```

>**IMPORTANT NOTE**: When sampling from the prior, the likelihood is not being calculated or estimated. In other words, the most time-consuming part of the MCMC does not take place. In that way, you should be able to gather enough samples with fewer runs than those needed when sampling from the posterior. E.g., at least 2 runs and, perhaps, no more than 5 chains, if you keep using the settings we have specified in the template control files in this tutorial.

Then, we will copy the directory `pipelines_MCMCtree` and will generate a copy called `pipelines_MCMCtree_prior` :

```sh
# Run from `main_protein` dir on your local
# PC. Please change directories until
# you are there. Then run the following
# commands.
cp -R pipelines_MCMCtree pipelines_MCMCtree_prior 
```

We will modify the bash script that will be submitted as a job array so the `userdata` option in the control file is equal to `0`, which enables `MCMCtree` to sample from the prior instead of sampling from the posterior (i.e., the data alignment is ignored). We will also set the nucleotide model to `0` (i.e., JC69) given that no data are being used and hence this option is not enabled when sampling from the prior. We will also use the strict clock (i.e., `clock = 1`, the default option in the template control file) given that we are not using the data to estimate evolutionary rates. Last, we will change the path to where the results will be stored:

```sh
# Run from `pipelines_MCMCtree_prior` dir on your local
# PC. Please change directories until
# you are there. Then run the following
# commands.
# 
# Prepare directories to sample from the prior,
# only one needed as nucleotide subsitution models 
# are not used.
rm -r ILN 
mv GBM CLK
cd CLK
mv pipeline_GBM.sh pipeline_CLK.sh

# Modify bash script: options `usedata` and `model`
sed -i "s/..*usedata..*/sed \-i \'s\/usedata\.\.\*\/usedata \= 0\/\' \$home\_dir\/mcmctree\_\$dir\"_r\"\$TASK\_ID\"\.ctl\"\nsed \-i \'s\/model\.\.\*\/model \= 0\/\' \$home\_dir\/mcmctree\_\$dir\"_\"\$TASK\_ID\"\.ctl\"/" *sh
# Modify clock model 
sed -i "s/^fi/fi\nsed \-i \'s\/clock\.\.\*\/clock \= 1\/\' \$home\_dir\/mcmctree\_\$dir\"_r\"\$TASK\_ID\"\.ctl\"/" *sh

# Modify path to save results 
sed -i 's/MCMCtree/MCMCtree\_prior/' *sh
sed -i 's/\/GBM\//\/CLK\//' *sh
# Comment soft link
sed -i 's/ln \-s/\#ln \-s/' *sh
```

Now, you can check that the lines have been correctly modified:

```sh
# Run from `pipelines_MCMCtree_prior` dir on your local
# PC. Please change directories until
# you are there. Then run the following
# commands.
grep 'usedata' */*sh
grep 'model'  */*sh
grep 'MCMCtree_prior' */*sh
```

## 2. Analyses with `MCMCtree` when sampling from the prior

Now, we will be able to run `MCMCtree` first when sampling from the prior (i.e., no data used!) using the code snippet below:

```sh
# Run from `pipelines_MCMCtree_prior/CLK` dir on your local
# PC. Please change directories until
# you are there. Then run the following
# commands.
chmod 775 *sh
./pipeline_CLK.sh &
```

### 2.1. Setting the file structure

We will now create a `sum_analyses_prot` directory to analyse the `MCMCtree` output.
First, we need to come back to the `01_analyses/02_MCMCtree` directory and run the following code snippet:

```sh
# Run from `01_analyses/02_MCMCtree` dir on your local
# PC. Please change directories until
# you are there. Then run the following
# commands. If you are running this code with your
# own analyses, make sure that you have correctly
# defined `num_aln` and `num_chains` variables with
# the correct values!
# Note that we will generate some directories for
# when the analyses when sampling from the posterior
# are ready!
mkdir sum_analyses_prot
cd sum_analyses_prot
num_aln=1
num_chains=5
for i in `seq 1 $num_aln`
do
for j in `seq 1 $num_chains`
do
mkdir -p posterior/$i/{GBM,ILN}/$j
mkdir -p prior/$i/CLK/$j
done
done

# Now, you can start copying the
# data we only need!
# Run the commands below as they are written if you
# are still inside `01_analyses/02_MCMCtree/sum_analyses_prot`).
# If not, run them once you are in this directory.
for i in `seq 1 $num_aln`
do
for j in `seq 1 $num_chains` 
do
printf "\n[[ Copying run "$j" for dataset "$i" into the corresponding directories ]]\n\n"
cp ../../../main/MCMCtree_prior/$j/CLK/$i/mcmc.txt prior/$i/CLK/$j/
cp ../../../main/MCMCtree_prior/$j/CLK/$i/*ctl prior/$i/CLK/$j/
done
done
```

### 2.2. MCMC diagnostics

Now that we have the output files from the different MCMC runs in an organised file structure, we are ready to run MCMC diagnostics!

We are going to run the R script [`MCMC_diagnostics_prior_prot.R`](scripts/MCMC_diagnostics_prior_prot.R) and follow the detailed step-by-step instructions detailed in the script. In a nutshell, the protocol will be the following:

1. Load the `mcmc.txt` files generated after each run.
2. Generate a convergence plot with the unfiltered chains.
3. Find whether there are major differences between the time estimates sampled across the chains for the same node sin the 97.5% and the 2.5% quantiles. If so, flag and delete said chains.
4. If some chains have not passed the filters mentioned above, create an object with the chains that have passed the filter.
5. Generate a new convergence plot with those chains that passed filters.
6. Calculate the effective sample size for all model parameters and check whether chain convergence has been reached with the chains that have passed filters.

The MCMC diagnostics did not find any of the chains problematic (run script [`MCMC_diagnostics_prior_prot.R`](scripts/MCMC_diagnostics_prior_prot.R)). Therefore, we used an in-house bash script, [`Combine_MCMC.sh`](scripts/Combine_MCMC.sh), to concatenate all the `mcmc.txt` files for the 2 chains in a unique file.

```sh
# Run from `02_MCMCtree/scripts`
cp Combine_MCMC.sh ../sum_analyses_prot/prior
# One argument taken: number of chains
cd ../sum_analyses_prot/prior
## Variables needed
## arg1 --> name of directory where analyses have taken place (e.g., 1)
## arg2 --> mcmc_files_CLK, mcmc_files_GBM, mcmc_files_ILN, etc.
## arg3 --> `seq 1 36`, "1 2 5", etc. | depends on whether some chains were filtered out or not
## arg4 --> ILN, GBM, CLK
dataset=1
./Combine_MCMC.sh $dataset mcmc_files_CLK "`seq 1 5`" CLK
```

The script above will generate a directory called `mcmc_files_CLK` inside the `prior` directory, where the `mcmc.txt` with the concatenated samples will be saved. A template script to generate the `FigTree.tre` file with this `mcmc.txt` has been saved inside the [`dummy_ctl_files`](dummy_ctl_files) directory.  We now will create a dummy alignment with only 2 nucleotides to generate the `FigTree` files using the concatenated `mcmc.txt` files. In order to do that, we can run the [`Generate_dummy_aln.R`](scripts/Generate_dummy_aln.R). Once you run it, a new directory called `dummy_aln` will be created, which will contain the dummy alignment. We have also generated dummy control file with option `print = -1`, which will not run an MCMC but, instead, will use the input files (file with the dummy alignment, calibrated tree file, and concatenated `mcmc.txt` file) to generate a `FigTree.tre` file with the mean estimated divergence times and the corresponding mean CIs using all the samples collected during all the MCMCs.

```sh
# Run from `sum_analyses_prot/prior`
base_dir=$( pwd )
cd ../../dummy_ctl_files
ctl_dir=$( pwd )
cd ../../../00_data/01_inp_data/
tt_dir=$( pwd )
name_tt=`ls *_calib.tree`
cd $ctl_dir
cd ../dummy_aln
aln_dir=$( pwd )
name_aln=`ls *aln`
cd $base_dir
cd mcmc_files_CLK
printf "[[ Generating tree file for concatenated \"mcmc.txt\"  ... ... ]]\n"
cp $ctl_dir/*ctl .
name_mcmc=`ls *mcmc.txt`
sed_aln=$( echo $aln_dir"/"$name_aln | sed 's/\//\\\//g' | sed 's/_/\\_/g' |  sed 's/\./\\\./g' )
sed_tt=$( echo $tt_dir"/"$name_tt | sed 's/\//\\\//g' | sed 's/_/\\_/g' |  sed 's/\./\\\./g' )
sed -i 's/MCMC/'${name_mcmc}'/' *ctl
sed -i -e 's/ALN/'${sed_aln}'/' *ctl
sed -i 's/TREE/'${sed_tt}'/' *ctl
# Run now MCMCtree after having modified the global vars according to the path to these files
# Then, rename the output tree file so we can easily identify later which tree belongs to
# which dataset easier
mcmctree *ctl
printf "\n"
mv FigTree.tre FigTree_CLK.tree
cd $base_dir
```

The next step is to plot the user-specified prior VS the effective prior. We have written the R script [`Check_priors_effVSuser.R`](scripts/Check_priors_effVSuser.R) to generate these plots. Once this script has finished, you will see that a new directory `plots/effVSuser/exdata` will have been created. Inside this directory, you will find one directory for each individual dataset with indiviudal plots for each node. In addition, all these plots have been merged into a unique document as well (note: some plots may be too small to see for each node, hence why we have generated individual plots).

Now, once the MCMC diagnostics have finished, you can extract the final data that you can use to write a manuscript as it follows:

```sh
# Run from `02_MCMCtree`
mkdir sum_files_prior_prot
cp -R sum_analyses_prot/prior/mcmc_files_CLK/*tree sum_files_prior_prot/
cp -R sum_analyses_prot/prior/1/CLK/CLK/*all_mean*tsv sum_files_prior_prot/
cp -R plots/ESS_and_chains_convergence/*prior_prot*pdf sum_files_prior_prot/
cp -R plots/effVSuser sum_files_prior_prot/
```

## 3. Analyses with `MCMCtree` when sampling from the posterior

Now that we have verified that there are no issues between the user-specified prior and the effective prior, we can run `MCMCtree` when sampling from the posterior. We will do these analyses under the GBM and ILN relaxed-clock models using the code snippet below:

```sh
# The command below will start a process in which `MCMCtree`
# will run when sampling from the posterior under the GBM. 
# E.g., you will see a line in the format of `[1] 902` will
# be printed after the command `./pipeline_CLK.sh &`. The `&`  
# allows you to keep using the same terminal to run more commands.
# You can use the same terminal now to run the next commands 
# or open a new terminal to continue.

# Now, go to directory `main_protein/pipeline_MCMCtree/GBM` dir on your local
# PC and run the following command. Please change directories until
# you are there.
chmod 775 *sh
./pipeline_GBM.sh &

# Now, go to directory `pipeline_MCMCtree/ILN` dir on your local
# PC and run the following command. Please change directories until
# you are there.
chmod 775 *sh
./pipeline_ILN.sh &
```

### 3.1. Setting the file structure

We will now go to the previously created `sum_analyses_prot` directory to analyse the `MCMCtree` output.
First, we need to come back to the `01_analyses/02_MCMCtree` directory and run the following code snippet:

```sh
# It is assumed that you have already
# generated the `posterior` directory 
# inside the `02_MCMCtree/sum_analyses_prot`
# directory following the commands
# detailed in `section 2.1`.
# Run the commands below 
# inside `01_analyses/02_MCMCtree/sum_analyses_prot`.
# Modify the variables `num_aln` and `num_chains`
# according to your dataset if you are using a different
# one from the example dataset in this tutorial
num_aln=1
num_chains=5
for i in `seq 1 $num_aln`
do
for j in `seq 1 $num_chains` 
do
printf "\n[[ Copying run "$j" for dataset "$i" into the corresponding directories ]]\n\n"
cp ../../../main_protein/MCMCtree/$j/GBM/$i/mcmc.txt posterior/$i/GBM/$j
cp ../../../main_protein/MCMCtree/$j/GBM/$i/*ctl posterior/$i/GBM/$j
cp ../../../main_protein/MCMCtree/$j/ILN/$i/mcmc.txt posterior/$i/ILN/$j
cp ../../../main_protein/MCMCtree/$j/ILN/$i/*ctl posterior/$i/ILN/$j
done
done
```

### 3.2. MCMC diagnostics

Now that we have the output files from the different MCMC runs in an organised file structure, we are ready to check the chains for convergence!

We are going to run the R script [`MCMC_diagnostics_posterior_prot.R`](scripts/MCMC_diagnostics_posterior_prot.R) and follow the detailed step-by-step instructions detailed in the script, which are essentially the same ones used when analysing the chains when sampling from the prior. If all the checks above have been passed, then you are ready to concatenate the parameter values sampled across the chains that passed the filters.

```sh
# Run from `02_MCMCtree/scripts`
cp Combine_MCMC.sh ../sum_analyses_prot/posterior
# One argument taken: number of chains
cd ../sum_analyses_prot/posterior
## Variables needed
## arg1 --> name of directory where analyses have taken place (e.g., 1)
## arg2 --> mcmc_files_CLK, mcmc_files_GBM, mcmc_files_ILN, etc.
## arg3 --> `seq 1 36`, "1 2 5", etc. | depends on whether some chains were filtered out or not
## arg4 --> ILN, GBM, CLK
dataset=1
./Combine_MCMC.sh $dataset mcmc_files_GBM "`seq 1 5`" GBM
./Combine_MCMC.sh $dataset mcmc_files_ILN "`seq 1 5`" ILN
```

Once the scripts above have finished, a new directory called `mcmc_files_[GBM|ILN]` will be created inside `posterior/1/[GBM|ILN]`, respectively. To map the mean time estimates with the filtered chains, we need to copy a control file, the calibrated Newick tree, and the dummy alignment we previously generated in section 2 inside this directory:

```sh
# Run from `sum_analyses_prot/posterior` directory.
# Please change directories until
# you are there. Then run the following
# commands.
base_dir=$( pwd )
cd ../../dummy_ctl_files
ctl_dir=$( pwd )
cd ../../../00_data/01_inp_data/
tt_dir=$( pwd )
name_tt=`ls *_calib.tree`
cd $ctl_dir
cd ../dummy_aln
aln_dir=$( pwd )
name_aln=`ls *aln`
cd $base_dir
cd mcmc_files_GBM
printf "[[ Generating tree file for concatenated \"mcmc.txt\" for GBM ... ... ]]\n"
cp $ctl_dir/*ctl .
name_mcmc=`ls *mcmc.txt`
sed_aln=$( echo $aln_dir"/"$name_aln | sed 's/\//\\\//g' | sed 's/_/\\_/g' |  sed 's/\./\\\./g' )
sed_tt=$( echo $tt_dir"/"$name_tt | sed 's/\//\\\//g' | sed 's/_/\\_/g' |  sed 's/\./\\\./g' )
sed -i 's/MCMC/'${name_mcmc}'/' *ctl
sed -i -e 's/ALN/'${sed_aln}'/' *ctl
sed -i 's/TREE/'${sed_tt}'/' *ctl
# Run now MCMCtree after having modified the global vars according to the path to these files
# Then, rename the output tree file so we can easily identify later which tree belongs to
# which dataset easier
mcmctree *ctl
mv FigTree.tre FigTree_GBM.tree
printf "\n"
cd $base_dir/mcmc_files_ILN
printf "[[ Generating tree file for concatenated \"mcmc.txt\" in "$data"/"$i" for ILN ... ... ]]\n"
cp $ctl_dir/*ctl .
name_mcmc=`ls *mcmc.txt`
sed_aln=$( echo $aln_dir"/"$name_aln | sed 's/\//\\\//g' | sed 's/_/\\_/g' |  sed 's/\./\\\./g' )
sed_tt=$( echo $tt_dir"/"$name_tt | sed 's/\//\\\//g' | sed 's/_/\\_/g' |  sed 's/\./\\\./g' )
sed -i 's/MCMC/'${name_mcmc}'/' *ctl
sed -i -e 's/ALN/'${sed_aln}'/' *ctl
sed -i 's/TREE/'${sed_tt}'/' *ctl
# Run now MCMCtree after having modified the global vars according to the path to these files
# Then, rename the output tree file so we can easily identify later which tree belongs to
# which dataset easier
mcmctree *ctl
mv FigTree.tre FigTree_ILN.tree
printf "\n"
cd $base_dir
```

Now, once the MCMC diagnostics have finished, you can extract the final data that we used to write our manuscript as it follows:

```sh
# Run from `02_MCMCtree`
mkdir sum_files_post_prot
cp -R sum_analyses_prot/posterior/mcmc_files_*/*tree sum_files_post_prot/
cp -R sum_analyses_prot/posterior/1/*/*/*all_mean*tsv sum_files_post_prot/
cp -R plots/ESS_and_chains_convergence/*posterior_prot*pdf sum_files_post_prot/
```
