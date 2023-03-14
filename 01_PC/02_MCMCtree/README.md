# Timetree inference

## 1. Running `MCMCtree`

Given that we already generated the calibrated tree, we have the alignment file (or files, if you are running this tutorial with more than one gene), and we have just generated the `in.BV` file with the gradient and the Hessian... We have everything we need to run `MCMCtree`!

However, we first need to create the corresponding file structure, which will follow the same structure as the directories where the Hessian and the gradient were calculated (see [the previous tutorial to run `BASEML`](../01_Hessian/README.md)). You will need to run the following code snippet:

```sh
# Run from `main` dir on your local
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

Now, we will use the same approach we used to generate the pipelines to run `CODEML`. In this case, however, the bash scripts have been modified accordingly. You will find them in the [`scripts` directory](scripts), from where you will copy them onto the `main/scripts` directory that you will have alreay generated.

```sh
# Run from the `01_PC/02_MCMCtree/scripts`
# dir on your local PC. Please change
# directories until you are there. Then run
# the following commands.
cp *sh ../../../main/scripts/
```

The `main` directory  will now have these two extra directories with the corresponding subdirectories:

```text
main
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
# Run from the `main` dir on your local
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
./generate_job_MCMCtree.sh $i GBM 1 $home_dir/pipelines_MCMCtree mcmctree $num_chains
./generate_job_MCMCtree.sh $i ILN 1 $home_dir/pipelines_MCMCtree mcmctree $num_chains
done
```

Now, before running `MCMCtree` to sample from the posterior, we will run `MCMCtree` but sampling from the prior. In this way, we can verify whether the user-specified prior (the probability distributions that we want the software as specified by us in the calibrated nodes) and the effective prior (the probability distribution the software will use). Oftentimes, to deal with truncation issues, the effective prior might differ from the user-specified prior and not all possible ages framed by the distribution specified by the user will be included (see an extensive analysis about this effect in [dos Reis et al. 2015](https://pubmed.ncbi.nlm.nih.gov/26603774/)). To assess whether this occurs in our analysis, we will run first `MCMCtree` without the data (i.e., `MCMCtree` will be sampling from the prior distribution and will only use the information provided by the calibrations included in the user tree).

First, we will generate a directory where `MCMCtree` will run and where all the results when sampling from the prior will be saved:

```sh
# Run from `main` dir on your local
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

Then, we will copy the directory `pipelines_MCMCtree` and will generate a copy called `pipelines_MCMCtree_prior` :

```sh
# Run from `main` dir on your local
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

Now, we will be able to run the scripts! You may want to open two terminals in each directory (i.e., `MCMCtree_prior` and `MCMCtree`) and simultaneously execute the two bash scripts:

```sh
# Run from `pipelines_MCMCtree_prior/CLK` dir on your local
# PC. Please change directories until
# you are there. Then run the following
# commands.
chmod 775 *sh
./pipeline_CLK.sh

# While a terminal with the script above is
# running, please open another terminal and then 
# run from `pipeline_MCMCtree/GBM` dir on your local
# PC. Please change directories until
# you are there.
chmod 775 *sh
./pipeline_GBM.sh

# While a terminal with the script above is
# running, please open another terminal and then 
# run from `pipeline_MCMCtree/ILN` dir on your local
# PC. Please change directories until
# you are there.
chmod 775 *sh
./pipeline_ILN.sh
```

## 2. Formatting the file structure to analyse the results

We will create a results directory to analyse the `MCMCtree` output.
Now, we will leave the `main` directory and come back to the `01_PC/02_MCMCtree` directory.

```sh
# Run from `01_PC/02_MCMCtree` dir on your local
# PC. Please change directories until
# you are there. Then run the following
# commands.
num_aln=1
for i in `seq 1 $num_aln`
do
for j in `seq 1 32`
do
mkdir -p posterior/$i/{GBM,ILN}/$j
mkdir -p prior/$i/CLK/$j
done
done

# Now, you can start copying the
# data we only need!
# It is assumed that you have already
# generated the `posterior` and
# `prior` directories following the 
# first set of commands above.
# Change `num_aln` if you are running
# this tutorial with a different dataset
# so it matches your requirements and
# do not change directories (i.e., you
# are still inside `01_PC/02_MCMCtree`).
num_aln=1
num_chains=5
for i in `seq 1 $num_aln`
do
for j in `seq 1 $num_chains` 
do
printf "\n[[ Moving run "$j" for dataset "$i" into the corresponding directories ]]\n\n"
cp ../../main/MCMCtree_prior/$j/CLK/$i/mcmc.txt prior/$i/CLK/$j/
cp ../../main/MCMCtree_prior/$j/CLK/$i/*ctl prior/$i/CLK/$j/
cp ../../main/MCMCtree/$j/GBM/$i/mcmc.txt posterior/$i/GBM/$j
cp ../../main/MCMCtree/$j/GBM/$i/*ctl posterior/$i/GBM/$j
cp ../../main/MCMCtree/$j/ILN/$i/mcmc.txt posterior/$i/ILN/$j
cp ../../main/MCMCtree/$j/ILN/$i/*ctl posterior/$i/ILN/$j
done
done
```

## 3. Generating summary output files

Now that we have the output files from the different MCMC runs in an organised file structure, we are ready to check the chains for convergence!

We are going to run the R script [`Calculate_ESS_posterior.R`](scripts/Calculate_ESS_posterior.R) and follow the detailed step-by-step instructions detailed in the script. In a nutshell, the protocol will be the following:

1. Load the `mcmc.txt` files generated after each run.
2. Generate a convergence plot with the unfiltered chains.
3. Find whether there are major differences between the time estimates sampled across the chains for the same node sin the 97.5% and the 2.5% quantiles. If so, flag and delete said chains.
4. Create an object with the chains that have passed the filter above.
5. Generate a new convergence plot with only the filtered chains.
6. Calculate the effective sample size for all model parameters and check whether chain convergence has been reached.

If all the checks above have been passed, then you are ready to concatenate the parameter values sampled across the chains that passed the filters. Once concatenated, the final tree topology with the mean time estimates can be generated. In order to do this, we will use the `Combine_MCMC_prior.sh` , `Combine_MCMC_GBM_posterior.sh`, and  `Combine_MCMC_ILN_posterior.sh` bash scripts. You can run them using the following commands:

```sh
# Run from `01_PC/02_MCMCtree/scripts` dir on your local
# PC. Please change directories until
# you are there. Then run the following
# commands.
# The scripts have two arguments: 
# [ IF NO CHAINS WERE DELETED ]
# Number of alignments and the number of chains
# run. 
# In this example, we have only one dataset 
# labelled `1` and all chains are to be kept
# when sampling from and the posterior, so 
# the first arg is `1` and the second `5`. 
#
# [ IF SOME CHAINS WERE FILTERED OUT ]
# Note that, if we had few chains, we should
# write the number and separate them with a
# comma (e.g., "2,5,7,8", without spaces and 
# without quotation marks). We should then 
# run the scripts that end with `_filt`, 
# e.g., `Combine_MCMC_prior_filt.sh`.
chmod 775 *sh
./Combine_MCMC_prior.sh 1 5
./Combine_MCMC_posterior_GBM.sh 1 5
./Combine_MCMC_posterior_ILN.sh 1 5
```

Once the scripts above have finished, a new directory called `mcmc_files_[CLK|GBM|ILN]` will be created inside `prior/1`, `posterior/1/[GBM|ILN]`, respectively. To map the mean time estimates with the filtered chains, we need to copy a control file, the calibrated newick tree, and a dummy alignment inside this directory. You can also copy the alignment you used but, if you are to save space, you can just generate a dummy alignment that contains an alignment with only two nucleotides for all taxa and a PHYLIP header that indicates so (i.e., second element a `2`). If you want to create the dummy alignment, please run the [`generate_dummy_aln.R`](scripts/Generate_dummy_aln.R), which will create a new directory called `dummy_aln` and, inside, the dummy alignment file called `dummy_aln.aln`.

Now, we will copy the control file and modify it so it can generate the corresponding final timetrees for the `prior` and `posterior` analyses:

```sh
# Run from `01_PC/02_MCMCtree` dir on your local
# PC. Please change directories until
# you are there. Then run the following
# commands.
home_dir=$( pwd )
num_aln=1
for i in `seq 1 $num_aln`
do
cp ../01_Hessian/control_files/prepbaseml.ctl prior/$i/mcmc_files_CLK/mcmctree_tree.ctl
cp ../01_Hessian/control_files/prepbaseml.ctl posterior/$i/mcmc_files_GBM/mcmctree_tree.ctl
cp ../01_Hessian/control_files/prepbaseml.ctl posterior/$i/mcmc_files_ILN/mcmctree_tree.ctl
name_tree=$( echo ls ../../00_data/01_inp_data/*_calib.tree | sed 's/..*\///' | sed 's/\//\\\//g' | sed 's/\_/\\\_/g' | sed 's/\./\\\./g' )
cd ../../00_data/01_inp_data
path_to_tree=$( pwd | sed 's/\//\\\//g' | sed 's/\_/\\\_/g' | sed 's/\./\\\./g' )
cd $home_dir/dummy_aln
name_aln=$( echo dummy\_aln\.aln )
path_to_aln=$( pwd | sed 's/\//\\\//g' | sed 's/\_/\\\_/g' | sed 's/\./\\\./g' )
cd $home_dir
sed -i 's/ALN/'${path_to_aln}'\/'${name_aln}'/' p*/$i/mcmc_files_*/*ctl
sed -i 's/TREE/'${path_to_tree}'\/'${name_tree}'/' p*/$i/mcmc_files_*/*ctl
sed -i 's/print\ \=\ 1/print\ \=\ \-1/' p*/$i/mcmc_files_*/*ctl
sed -i 's/usedata..*/usedata\ \=\ 0/' p*/$i/mcmc_files_*/*ctl
done
```

> **NOTE**: If the commands below do not work,  PAML v4.9h has this option implemented. If you are using the latest release of PAML and the commands below do not work, please raise an issue on [the PAML GitHub repository](https://github.com/abacus-gene/paml) regarding option `print = -1` for `MCMCtree`.

Now, you are ready to run the control file inside `prior/1/mcmc_files_CLK`, `posterior/1/mcmc_files_GBM`, and `posterior/1/mcmc_files_ILN` to generate the final timetree! You should just navigate to these directories and execute `MCMCtree` using the command that you use on your PC. E.g.:

```sh
# Run from `01_PC/02_MCMCtree/prior/1/`
# and then from the corresponding
# directories when `MCMCtree` was
# run when sampling from the posterior.
# If your executable file is called
# `mcmctree`, run the following
# command. Otherwise, modify the 
# command to run `MCMCtree` on your PC
# accordingly
cd prior/1/mcmc_files_CLK
mcmctree *ctl
cd ../../../posterior/1/mcmc_files_GBM/
mcmctree *ctl
cd ../mcmc_files_GBM/
mcmctree *ctl
```

And that should be it: you now have the final timetree under each relaxed-clock model!

## EXERCISE (hard)

How would you now plot the effective prior versus the user-specified prior to make sure that there is no conflict between the calibrations specified by the user and the ones used by `MCMCtree`?
