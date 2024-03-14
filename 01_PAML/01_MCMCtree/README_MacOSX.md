# Bayesian inference of species divergences

## 1. Setting the file structure to run `MCMCtree`

Given that we already generated the calibrated tree, the partitioned alignment file, and we have just generated the `in.BV` file... We have everything we need to run `MCMCtree`!

First, we will create the file structure required for the timetree inference analyses using the following code snippet:

```sh
# Run from `example_dating` dir.
# Please change directories until
# you are there. Then, run the following
# commands.
num_aln=1 
num_chains=6 # num chains we will run
for j in `seq 1 $num_chains`
do
for i in `seq 1 $num_aln`
do
mkdir -p MCMCtree/$j/{GBM,ILN}/$i
done
done
# pipelines dir
for i in `seq 1 $num_aln`
do
mkdir -p pipelines_MCMCtree/$i/{GBM,ILN}
done
```

The `example_dating` directory will now have these two extra directories with the corresponding subdirectories:

```text
example_dating
  |- MCMCtree
  |    |- [1-6] #6 chains
  |         |- [GBM|ILN]
  |              |- [1] # 1 dataset, 1 dir
  |                
  |- pipelines_MCMCtree
       |- [1] # One for each dataset
            |- [GBM|ILN]/
```

>**IMPORTANT NOTE**: When sampling from the posterior, the likelihood is being calculated or approximated, depending on the `userdata` option you set in the control file to run `MCMCtree`. The larger the dataset, the more time it will take for `MCMCtree` to finish regardless of the option you use, although the approximation will be faster than the exact likelihood calculation.

Now, we will copy our in-house bash scripts to run `MCMCtree` and corresponding template files with which we will parse our data (available under [`scripts`](scripts)) to the `example_dating` directory:

```sh
# Run from `01_PAML/01_MCMCtree/scripts`
cp *MCMCtree*sh ../../../example_dating/scripts
```

You can now go back to your `example_dating` directory and run the `generate_job_MCMCtree.sh` (if cluster) or `generate_job_MCMCtree_PC_MacOSX.sh` (if local PC) script, one of our in-house bash scripts mentioned above, using the commands below:

```sh
# Run from the `example_dating` dir.
# Please change directories until
# you are there. Then run the following
# commands.
home_dir=$( pwd )
cd scripts
chmod 775 *sh
num_aln=1
num_chains=6
# Arg 1: Number of gene alignments being analysed.
# Arg 2: Clock model (e.g., "GBM", "ILN", or "CLK).
# Arg 3: Number of partitions in the alignment file.
# Arg 4: Path to the pipeline directory.
# Arg 5: Command to execute MCMCtree (e.g., `mcmctree`, `mcmctree_2000`, etc.)
# Arg 6: Number of chains run.
# Arg 7: Name of working directory (e.g., `example_dating`)
# Arg 8: Boolean, enable duplication option? 0: no, 1: yes
# Arg 9: Boolean, PAML exported to the path? `Y` or `N`.
#        If `N`, the executable file will be required to be in the home dirctory,
#        i.e., directory which name you type as `Arg3`.
#
for i in `seq 1 $num_aln`
do
./generate_job_MCMCtree_PC_MacOSX.sh $i GBM 1 $home_dir/pipelines_MCMCtree mcmctree $num_chains example_dating 0 Y
./generate_job_MCMCtree_PC_MacOSX.sh $i ILN 1 $home_dir/pipelines_MCMCtree mcmctree $num_chains example_dating 0 Y
done
# Check paths were correctly incorporated!
cd ../pipelines_MCMCtree
grep '^dir=\$' */*/*sh
grep 'MCMCtree' */*/*sh
grep 'alignments' */*/*sh
grep 'Hessian' */*/*sh
grep 'ctl_dir' */*/*sh
```

Now, before running `MCMCtree` to sample from the posterior, we will run `MCMCtree` but sampling from the prior. In this way, we can verify whether the user-specified priors constraining some node ages in our fixed tree topology (i.e., probability distributions that we want the dating program to use) and the effective prior (i.e., probability distributions the dating program will use after building the joint prior distribution, which considers not only the user-specified priors but also additional priors such as the birth-death process for uncalibrated nodes, the rate prior, etc.). Oftentimes, truncation issues may arise when the effective priors are in disagreement with the user-specified priors (see an extensive study about this effect in [dos Reis et al. 2015](https://pubmed.ncbi.nlm.nih.gov/26603774/)). To assess the effect of truncation prior to analysing our dataset, we will run `MCMCtree` without the data (i.e., `MCMCtree` will be sampling from the prior distribution, and thus sequence data will not be used).

First, we will generate a directory where `MCMCtree` will run when sampling from the prior:

```sh
# Run from `example_dating` dir.
# Please change directories until
# you are there. Then run the following
# commands.
num_aln=1
num_chains=6
for j in `seq 1 $num_chains`
do
for i in `seq 1 $num_aln`
do
mkdir -p MCMCtree_prior/$j/CLK/$i
done
done
```

>**IMPORTANT NOTE**: When sampling from the prior, the likelihood is not being calculated or estimated. In other words, the most time-consuming part of the MCMC does not take place. In that way, you should be able to gather enough samples with fewer runs than those needed when sampling from the posterior.

Then, we will copy the directory `pipelines_MCMCtree` and will generate a copy called `pipelines_MCMCtree_prior`:

```sh
# Run from `example_dating` dir.
# Please change directories until
# you are there. Then run the following
# commands.
cp -R pipelines_MCMCtree pipelines_MCMCtree_prior
cd pipelines_MCMCtree_prior
```

We will modify the bash script that will be submitted as a job array so that the `userdata` option in the control file is equal to `0`, which enables `MCMCtree` to sample from the prior instead of sampling from the posterior (i.e., the alignment file is ignored). In that way, the rest of the setting concerning the evolutionary model will not be enabled by `MCMCtree` as the sequence data are not being used. Last, we will change the path to where the results will be stored:

```sh
# Run from `pipelines_MCMCtree_prior`.
# Please change directories until
# you are there. Then run the following
# commands.
# 
# Prepare directories to sample from the prior,
# only one needed as nucleotide subsitution models 
# are not used.
#
num_aln=1
home_dir=$( pwd )
for i in `seq 1 $num_aln`
do
cd $home_dir/$i
rm -r ILN 
mv GBM CLK
cd CLK
mv pipeline_GBM.sh pipeline_CLK.sh
done
cd $home_dir

# Modify bash script: options `usedata` and `model`
sed -i '' "s/..*usedata..*/sed \-i \'\' \'s\/usedata\.\.\*\/usedata \= 0\/\' \$home\_dir\/mcmctree\_\$dir\"_r\"\$SGE\_TASK\_ID\"\.ctl\"\nsed \-i \'\' \'s\/model\.\.\*\/model \= 4\/\' \$home\_dir\/mcmctree\_\$dir\"_r\"\$SGE\_TASK\_ID\"\.ctl\"/" */*/*sh
# Modify clock model 
sed -i '' "s/^fi/fi\nsed \-i \'\' \'s\/clock\.\.\*\/clock \= 1\/\' \$home\_dir\/mcmctree\_\$dir\"_r\"\$SGE\_TASK\_ID\"\.ctl\"/" */*/*sh

# Modify path to save results 
sed -i '' 's/\/GBM\//\/CLK\//' */*/*sh
sed -i '' 's/MCMCtree/MCMCtree\_prior/' */*/*sh 
# Comment soft link
sed -i '' 's/ln \-s/\#ln \-s/' */*/*sh
# Not relevant in this example but, if you had
# more chains in the posterior than in the prior,
# please modify the corresponding lines in your script!
```

Now, you can check that the lines have been correctly modified:

```sh
# Run from `pipelines_MCMCtree_prior`.
# Please change directories until
# you are there. Then run the following
# commands.
grep '^dir=' */*/*sh
grep 'usedata' */*/*sh
grep 'model'  */*/*sh
grep 'MCMCtree' */*/*sh
```

## 2. Analyses with `MCMCtree` when sampling from the prior

### Running `MCMCtree` (prior)

Now, we will be able to run `MCMCtree` first when sampling from the prior (i.e., no data used!) using the code snippet below:

```sh
# Run from `pipelines_MCMCtree_prior/1/CLK`.
# Please change directories until
# you are there. Then run the following
# commands.
chmod 775 *sh
./pipeline_CLK.sh &
```

### Setting the file structure to analyse `MCMCtree` output - prior

We will now create a `sum_analyses` directory to analyse the `MCMCtree` output. First, we need to create a file structure compatible with the tools you will subsequently use to summarise the data:

```sh
# Run everything from `example_dating` dir
num_chains=6
num_datasets=1
mkdir -p tmp_to_transfer/00_prior
cd tmp_to_transfer
for j in `seq 1 $num_datasets`
do
for i in `seq 1 $num_chains`
do
mkdir -p 00_prior/CLK/$j/$i/
# Now, copy the files that are required for sum stats.
# We have run 6 chains for analyses sampling from the prior
printf "\n[[ Copying run "$i" for analyses sampling from the prior -- dataset "$j" ]]\n\n"
cp ../MCMCtree_prior/$i/CLK/$j/mcmc.txt 00_prior/CLK/$j/$i
cp ../MCMCtree_prior/$i/CLK/$j/*ctl 00_prior/CLK/$j/$i
cp ../MCMCtree_prior/$i/CLK/$j/SeedUsed 00_prior/CLK/$j/$i
done
grep 'Species tree for FigTree' -A1 ../MCMCtree_prior/$j/CLK/1/out.txt | sed -n '2,2p' > 00_prior/node_tree_$j.tree
done
```

Now, you can transfer these temporary directories to your `01_PAML/01_MCMCtree` directory on your local PC, e.g., using `rsync`:

```sh
# Run from `01_PAML/01_MCMCtree` dir. 
# Please change directories until
# you are there. Then run the following
# commands. If you are running this code with your
# own analyses, make sure that you have correctly
# defined `num_aln` and `num_chains` variables with
# the correct values!
# Note that we will generate some directories for
# when the analyses when sampling from the posterior
# are ready!
cp -R ../../../example_dating/tmp_to_transfer/00_prior .
```

### MCMC diagnostics - prior

Now that we have the output files from the different MCMC runs in an organised file structure, we are ready to run MCMC diagnostics!

First, however, we need to generate a file with calibration information that is compatible with the subsequent scripts. For that purpose, we can use our in-house R script [`Merge_node_labels.R](scripts/Merge_node_labels.R), which will generate one calibration file for each dataset analysed in case there are differences in the tree topology/ies and the node labels which age is being constrained.

Now, we can run the R script [`MCMC_diagnostics_prior.R`](scripts/MCMC_diagnostics_prior.R) and follow the detailed step-by-step instructions detailed in the script. In a nutshell, the protocol will be the following:

1. Load the `mcmc.txt` files generated after each run.
2. Generate a convergence plot with the unfiltered chains.
3. Find whether there are major differences between the time estimates sampled across the chains for the same node sin the 97.5% and the 2.5% quantiles. If so, flag and delete said chains.
4. If some chains have not passed the filters mentioned above, create an object with the chains that have passed the filters.
5. Generate a new convergence plot with those chains that passed the filters.
6. Calculate Rhat, tail-ESS, and bulk-ESS to check whether chain convergence has been reached with the chains that have passed filters.

The MCMC diagnostics did not find any of the chains problematic after running [our in-house R script `MCMC_diagnostics_prior.R`](scripts/MCMC_diagnostics_prior.R). Therefore, we used [our in-house bash script `Combine_MCMC_MacOSX.sh`](scripts/Combine_MCMC_MacOSX.sh) to concatenate the `mcmc.txt` files for the 6 chains that passed the filters in a unique file.

```sh
# Run from `01_MCMCtree/scripts`
chmod 775 *sh
cp Combine_MCMC_MacOSX.sh ../sum_analyses/00_prior
# One argument taken: number of chains
cd ../sum_analyses/00_prior
## Variables needed
## arg1 --> name of directory where analyses have taken place (e.g., CLK, GBM, ILN)
## arg2 --> output dir: mcmc_files_CLK, mcmc_files_GBM, mcmc_files_ILN, etc.
## arg3 --> "`seq 1 36`", "1 2 5", etc. | depends on whether some chains were filtered out or not
## arg4 --> clock model used: ILN, GBM, CLK
## arg5 --> number of samples specified to collect in the control file to run `MCMCtree`
dataset=$( echo CLK )
num_dat=1
name_dat=( 'example' )
count=-1 #  counter will start at 0 in first iteration!
for i in `seq 1 $num_dat`
do
count=$(( count + 1 ))
./Combine_MCMC_MacOSX.sh $dataset/$i mcmc_files_${name_dat[count]}_CLK "`seq 1 6`" CLK 20000
done
```

The script above will generate three directories (one for each dataset) inside the `00_prior` directory, where the `mcmc.txt` with the concatenated samples will be saved. In addition, a directory called `mcmcf4traces` will also be generated so that formatted MCMC files compatible with programs such as `Tracer` can be used to check for chain convergence. A template script to generate the `FigTree.tre` file with this `mcmc.txt` has been saved inside the [`dummy_ctl_files`](dummy_ctl_files) directory.

We will now create a dummy alignment with only 2 AAs to generate the `FigTree` files using the concatenated `mcmc.txt` files. In order to do that, we can run the [`Generate_dummy_aln.R`](scripts/Generate_dummy_aln.R). Once you run it, a new directory called `dummy_aln` will be created, which will contain the dummy alignment.

We have also generated dummy control file with option `print = -1`, which will not run an MCMC but, instead, will use the input files (file with the dummy alignment, calibrated tree file, and concatenated `mcmc.txt` file) to generate a `FigTree.tre` file with the mean estimated divergence times and the corresponding mean CIs using all the samples collected during all the MCMCs.

```sh
# Run from `sum_analyses/00_prior`
name_dat=( 'example' )
num_dat=1
count=-1
for i in `seq 1 $num_dat`
do
count=$(( count + 1 ))
printf "\n[[ Analysing dataset "${name_dat[count]}" ]]\n"
base_dir=$( pwd )
cd ../../dummy_ctl_files
ctl_dir=$( pwd )
cd ../../../00_data_formatting/01_inp_data/
tt_dir=$( pwd )
name_tt=`ls ${name_dat[count]}"_calib_MCMCtree.tree"`
cd $ctl_dir
cd ../dummy_aln
aln_dir=$( pwd )
name_aln=`ls *aln`
cd $base_dir
cd mcmc_files_${name_dat[count]}"_CLK"
printf "[[ Generating tree file for concatenated \"mcmc.txt\"  ... ... ]]\n"
cp $ctl_dir/*ctl .
name_mcmc=`ls *mcmc.txt`
sed_aln=$( echo $aln_dir"/"$name_aln | sed 's/\//\\\//g' | sed 's/_/\\_/g' |  sed 's/\./\\\./g' )
sed_tt=$( echo $tt_dir"/"$name_tt | sed 's/\//\\\//g' | sed 's/_/\\_/g' |  sed 's/\./\\\./g' )
sed -i '' 's/MCMC/'${name_mcmc}'/' *ctl
sed -e -i '' 's/ALN/'${sed_aln}'/' *ctl
sed -i '' 's/TREE/'${sed_tt}'/' *ctl
# Run now MCMCtree after having modified the global vars according to the path to these files
# Then, rename the output tree file so we can easily identify later which tree belongs to
# which dataset easier
mcmctree *ctl
mv FigTree.tre FigTree_${name_dat[count]}"_CLK.tree"
printf "\n"
cd $base_dir
done
```

The next step is to plot the calibration density (commonly referred to as "user-specified prior") VS the marginal density (also known as "effective prior"). We used our in-house R script [`Check_priors_calVSmarg.R`](scripts/Check_priors_calVSmarg.R) to generate these plots. If you are to run this script with other datasets, however, make sure that your "hard bounds" are not `0.000` in the `Calibnodes_*csv` files and, instead, they are `1e-300` (i.e., while 1e-300 is rounded to `0.000` in the `MCMCtre` output, which can be used to generate the csv files aforementioned, we need `1e-300` to plot distributions in R). To make sure this was not affecting our csv files, we ran the following code snippet:

```sh
# Run from `01_MCMCtree/calib_files`
sed -i 's/0\.000/1e\-300/g' *csv
```

Once this script has finished, you will see that a new directory `plots/calVSmarg` will have been created. Inside this directory, you will find one directory for each individual dataset with individual plots for each node. In addition, all these plots have been merged into a unique document as well (note: some plots may be too small to see for each node, hence why we have generated individual plots).

Now, once the MCMC diagnostics have finished, you can extract the final data that you can use to write a manuscript as it follows:

```sh
# Run from `01_MCMCtree`
mkdir sum_files_prior
cp -R sum_analyses/00_prior/mcmc_files*CLK/*CLK*tree sum_files_prior/
cp -R sum_analyses/00_prior/CLK/*/*/*all_mean*tsv sum_files_prior/
cp -R plots/ESS_and_chains_convergence/*prior*pdf sum_files_prior/
cp -R plots/calVSmarg sum_files_prior/
```

## 3. Analyses with `MCMCtree` when sampling from the posterior

### Running `MCMCtree` (posterior)

Now that we have verified that there are no issues between the user-specified prior and the effective prior, we can run `MCMCtree` when sampling from the posterior. We will do these analyses under the GBM and ILN relaxed-clock models using the code snippet below:

```sh
# Now, go to directory `example_dating/pipelines_MCMCtree/1/GBM` dir 
# and run the following command. Please change directories until
# you are there.
chmod 775 *sh
./pipeline_GBM.sh &
cd ../ILN
chmod 775 *sh
./pipeline_ILN.sh &
```

### Setting the file structure to analyse `MCMCtree` output - posterior

We will now create a directory inside the `sum_analyses` directory to analyse the `MCMCtree` output as we did with the results obtained when sampling from the prior:

```sh
# Run from `example_dating` directory
# We have run 6 chains for analyses sampling from the posterior.
# Therefore, `i` will go form 1 to 6
cd tmp_to_transfer
num_chains=6
num_datasets=1
# The `01_posterior` should already exist from the previous
# analyses. If not, it will be created during the `for` loop
for j in `seq 1 $num_datasets`
do
for i in `seq 1 $num_chains` 
do
mkdir -p 01_posterior/{GBM,ILN}/$j/$i
printf "\n[[ Copying run "$i" for analyses sampling from the posterior -- dataset "$j" ]]\n\n"
cp ../MCMCtree/$i/GBM/$j/mcmc.txt 01_posterior/GBM/$j/$i
cp ../MCMCtree/$i/GBM/$j/*ctl* 01_posterior/GBM/$j/$i
cp ../MCMCtree/$i/GBM/$j/SeedUsed 01_posterior/GBM/$j/$i
cp ../MCMCtree/$i/ILN/$j/mcmc.txt 01_posterior/ILN/$j/$i
cp ../MCMCtree/$i/ILN/$j/*ctl 01_posterior/ILN/$j/$i
cp ../MCMCtree/$i/ILN/$j/SeedUsed 01_posterior/ILN/$j/$i
done
done
```

Now, we will copy this file structure to directory `01_PAML/01_MCMCtre`:

```sh
# Run from `01_PAML/01_MCMCtree` dir on your local
# PC. Please change directories until
# you are there. Then run the following
# commands. If you are running this code with your
# own analyses, make sure that you have correctly
# defined `num_aln` and `num_chains` variables with
# the correct values!
# Note that we will generate some directories for
# when the analyses when sampling from the posterior
# are ready!
cd sum_analyses
cp -R ../../../example_dating/tmp_to_transfer/01_posterior/ .
```

### MCMC diagnostics - posterior

Now that we have the output files from the different MCMC runs in an organised file structure, we are ready to check for chain convergence!

We are going to run the R script [`MCMC_diagnostics_posterior.R`](scripts/MCMC_diagnostics_posterior.R) and follow the detailed step-by-step instructions detailed in the script, which are essentially the same ones used when analysing the chains when sampling from the prior. Given that no problems have been found with any of the chains we ran, we are ready to concatenate the parameter values sampled across the 6 independent chains we ran:

```sh
# Run from `01_MCMCtree/scripts`
cp Combine_MCMC_MacOSX.sh ../sum_analyses/01_posterior
# One argument taken: number of chains
cd ../sum_analyses/01_posterior
## Variables needed
## arg1 --> name of directory where analyses have taken place (e.g., CLK, GBM, ILN)
## arg2 --> output dir: mcmc_files_CLK, mcmc_files_GBM, mcmc_files_ILN, etc.
## arg3 --> "`seq 1 36`", "1 2 5", etc. | depends on whether some chains were filtered out or not
## arg4 --> clock model used: ILN, GBM, CLK
## arg5 --> number of samples specified to collect in the control file to run `MCMCtree`
dirname_1=GBM
dirname_2=ILN
num_dat=1
name_dat=( 'example' )
count=-1 #  counter will start at 0 in first iteration!
for i in `seq 1 $num_dat`
do
count=$(( count + 1 ))
./Combine_MCMC_MacOSX.sh $dirname_1/$i mcmc_files_${name_dat[count]}_GBM "`seq 1 6`" GBM 20000
./Combine_MCMC_MacOSX.sh $dirname_2/$i mcmc_files_${name_dat[count]}_ILN "`seq 1 6`" ILN 20000
done
```

Once the scripts above have finished, directories called `mcmc_files_part_[GBM|ILN]` and `mcmcf4traces` will be created inside `01_posterior/`. To map the mean time estimates with the filtered chains, we need to copy a control file, the calibrated Newick tree, and the dummy alignment we previously generated when analysing the results when sampling from the prior:

```sh
# Run from `sum_analyses_prot/01_posterior` directory.
# Please change directories until
# you are there. Then run the following
# commands.
name_dat=( 'example' )
num_dat=1
count=-1
for i in `seq 1 $num_dat`
do
printf "\n[[ ANALYSING DATASET "${name_dat[count]}" ]]\n\n"
count=$(( count + 1 ))
base_dir=$( pwd )
cd ../../dummy_ctl_files
ctl_dir=$( pwd )
cd ../../../00_data_formatting/01_inp_data/
tt_dir=$( pwd )
name_tt=`ls ${name_dat[count]}"_calib_MCMCtree.tree"`
cd $ctl_dir
cd ../dummy_aln
aln_dir=$( pwd )
name_aln=`ls *aln`
cd $base_dir
cd mcmc_files_${name_dat[count]}"_GBM"
printf "[[ Generating tree file for concatenated \"mcmc.txt\" for GBM ... ... ]]\n"
cp $ctl_dir/*ctl .
name_mcmc=`ls *mcmc.txt`
sed_aln=$( echo $aln_dir"/"$name_aln | sed 's/\//\\\//g' | sed 's/_/\\_/g' |  sed 's/\./\\\./g' )
sed_tt=$( echo $tt_dir"/"$name_tt | sed 's/\//\\\//g' | sed 's/_/\\_/g' |  sed 's/\./\\\./g' )
sed -i '' 's/MCMC/'${name_mcmc}'/' *ctl
sed -e -i '' 's/ALN/'${sed_aln}'/' *ctl
sed -i '' 's/TREE/'${sed_tt}'/' *ctl
# Run now MCMCtree after having modified the global vars according to the path to these files
# Then, rename the output tree file so we can easily identify later which tree belongs to
# which dataset easier
mcmctree *ctl
mv FigTree.tre FigTree_${name_dat[count]}"_GBM.tree"
printf "\n"
cd $base_dir/mcmc_files_${name_dat[count]}"_ILN"
printf "[[ Generating tree file for concatenated \"mcmc.txt\" in "$data"/"$i" for ILN ... ... ]]\n"
cp $ctl_dir/*ctl .
name_mcmc=`ls *mcmc.txt`
sed_aln=$( echo $aln_dir"/"$name_aln | sed 's/\//\\\//g' | sed 's/_/\\_/g' |  sed 's/\./\\\./g' )
sed_tt=$( echo $tt_dir"/"$name_tt | sed 's/\//\\\//g' | sed 's/_/\\_/g' |  sed 's/\./\\\./g' )
sed -i '' 's/MCMC/'${name_mcmc}'/' *ctl
sed -e -i '' 's/ALN/'${sed_aln}'/' *ctl
sed -i '' 's/TREE/'${sed_tt}'/' *ctl
# Run now MCMCtree after having modified the global vars according to the path to these files
# Then, rename the output tree file so we can easily identify later which tree belongs to
# which dataset easier
mcmctree *ctl
mv FigTree.tre FigTree_${name_dat[count]}"_ILN.tree"
printf "\n"
cd $base_dir
done
```

Now, once the MCMC diagnostics have finished, we can run our [in-house R script](scripts/Check_priors_VS_posteriors.R) to plot the posterior distributions against the prior distributions, which can help to better assess how informative the data are and whether there are any serious contradictions between the prior and the posterior distributions.

Last, you can extract the final data that we used to write our manuscript as it follows:

```sh
# Run from `01_MCMCtree`
mkdir sum_files_post
cp -R sum_analyses/01_posterior/mcmc_files_*/FigTree*tree sum_files_post/
cp -R sum_analyses/01_posterior/*/*/*/*all_mean*tsv sum_files_post/
cp -R plots/priorVSpost*pdf sum_files_post/
cp -R plots/ESS_and_chains_convergence/*post*pdf sum_files_post/
```
