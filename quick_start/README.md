# Timetree inference with `PAML` programs

## Input data

We will be using a modified version of the primate mitochondrial dataset available from the [`examples` directory in the `PAML` GitHub repository](https://github.com/abacus-gene/paml/blob/master/examples) (i.e., [original input sequence file](https://github.com/abacus-gene/paml/blob/master/examples/DatingSoftBound/mtCDNApri123.txt) and [original input tree file](https://github.com/abacus-gene/paml/blob/master/examples/DatingSoftBound/mtCDNApri.trees)). Please note that I made some changes to these files:

* You will see that I changed the calibration notation in [our tree file](00_inp_data/mtCDNApri.trees) so that it is clearer which bounds are being used to constrain node ages (i.e., 2 soft-bound calibrations and 1 upper-bound calibration). Note that I incorporated the root age constraint (i.e., upper-bound calibration) in the Newick tree as recommended when using `MCMCtree` (i.e., using option `RootAge` in the control file is somewhat discouraged).
* I have generated one file for each of the 3 alignment blocks available in the [original sequence file](https://github.com/abacus-gene/paml/blob/master/examples/DatingSoftBound/mtCDNApri123.txt). We will be estimating the branch lengths, the gradient, and the Hessian with `BASEML` for each of the alignment blocks, which is what `MCMCtree` will use to approximate the likelihood calculation ([dos Reis and Yang, 2011](https://academic.oup.com/mbe/article/28/7/2161/1051613)). These files are named `mtCDNApri_part[1-3].txt` and saved in directory [00_inp_data](00_inp_data).
* Lastly, I used the [control file in the `PAML` GitHub repository](https://github.com/abacus-gene/paml/blob/master/examples/DatingSoftBound/mcmctree.ctl) as a template to create our [own template control file for this tutorial](01_PAML/mcmctree_orig.ctl).

You will not need to further process the aforementioned input files in directory [`00_inp_data`](00_inp_data/). Nevertheless, you can later try the [tutorial on reproducible timetree inference](../README.md#tutorial-on-reproducibility), which starts [from data formatting](../00_data_formatting), if you want to see common filtering steps users tend to follow before running `PAML` programs.

Below, you will find the steps you need to follow to run `PAML` programs `BASEML` and `MCMCtree` for timetree inference using the approximate likelihood calculation.

## `PAML` programs

### Running `BASEML`

We will first prepare the input files required for `BASEML`, which you can do by running `MCMCtree` with settings `usedata = 3`. Note that, once the input files are ready, we will run `BASEML` three times: one for each alignment block.

We will use the code below to generate three copies of the [`mcmctree_orig.ctl` file](01_PAML/mcmctree_orig.ctl) (i.e., one for each alignment block), which we will modify with command `sed` to tailor the settings for the analyses we want to run:

```sh
## Run from directory `quick_start`
# Please use your terminal to navigate to 
# this directory, then copy the commands below
# on your terminal
cd 01_PAML
for i in `seq 1 3`
do
mkdir -p 00_BASEML/$i
cp mcmctree_orig.ctl 00_BASEML/$i/mcmctree_baseml_part$i.ctl
# Modify the path to the input sequence file,
# the alignment blocks
sed -i 's/seqfile..*/seqfile = \.\.\/\.\.\/.\.\/00\_inp\_data\/mtCDNApri_part'${i}'\.txt/' 00_BASEML/$i/mcmctree_baseml_part$i.ctl
# Modify the path to the input tree file
sed -i 's/treefile\ \=\ /treefile\ \=\ \.\.\/\.\.\/.\.\/00\_inp\_data\//' 00_BASEML/$i/mcmctree_baseml_part$i.ctl
# Specify `usedata = 3`, which will enable the settings
# for generating the `in.BV` file for approximating
# likelihood calculation
sed -i 's/usedata..*/usedata\ \=\ 3/' 00_BASEML/$i/mcmctree_baseml_part$i.ctl
# Modify option `ndata` to make sure you will have
# only one alignment block (i.e., remember that you are
# creating one control file per alignment block within 
# this `for` loop!)
sed -i 's/ndata..*/ndata\ \=\ 1/' 00_BASEML/$i/mcmctree_baseml_part$i.ctl
done
```

Now, we have everything we need to run `BASEML` in directories `00_BASEML/1`, `00_BASEML/2`, and `00_BASEML/3`:

```sh
## Run from `01_PAML`
# You should still be in this directory 
# but, if not, please change directories until
# you are there. Then, run the following
# commands.
cd 00_BASEML
home_dir=$( pwd )
for i in `seq 1 3`
do
cd $home_dir/$i
printf "[[ Analysing alignment block "$i" ]]\n"
# Run `MCMCtree` to prepare input files
# for `BASEML` (i.e., you will not run
# timetree inference now, `MCMCtree` is
# just in charge of creating the input
# files for `BASEML`!)
mcmctree *ctl
done
cd $home_dir
```

Note that, when we ran the commands above, we were not interested in running `BASEML` or `MCMCtree` until "the end". We just wanted to execute `MCMCtree` with option `usedata = 3` so that this program generates the `tmp000*` files that `BASEML` will then need as input files to estimate the branch lengths, the gradient, and the Hessian. We carry out this analysis in two steps so that we can replace option `method = 0` with `method = 1` in the `tmp0001.ctl` file that will be output with the commands above. As explained in the [`PAML` documentation](https://github.com/abacus-gene/paml/blob/master/doc/pamlDOC.pdf) (at the time of writing, page 56), the iteration algorithm enabled when setting `method = 1` is much more efficient with large datasets than the algorithm enabled when setting `method = 0`. While this is not a very large dataset, it is good practice to divide these steps into two so that you can always check the `BASEML` settings in the automatically created `tmp0001.ctl` files!

Once we see the `tmp0001.ctl` file generated, we can therefore kill the job by pressing `ctrl+C`. In this case, the dataset is small, and so you will not even have time to press `ctrl+C`! When analysing large datasets, you may be stuck for a while after you see the following lines, which is when you can kill the job as the `tmp*` files will appear:

```text
*** Locus 1 ***
running baseml tmp0001.ctl
BASEML in paml version 4.10.7, June 2023
ns = 7          ls = 205
Reading sequences, sequential format..
Reading seq # 7: gibbon
Sequences read..

205 site patterns read, 3331 sites
Counting frequencies..

      252 bytes for distance
    39360 bytes for conP
     8200 bytes for fhK
  8000000 bytes for space
```

After you have killed the job, you can run the following commands to make sure that the correct settings to run `BASEML` are enabled:

```sh
## Run from the `00_BASEML`.
# You should still be in this directory 
# but, if not, please change directories until
# you are there. Then, run the following
# commands.
sed -i 's/method\ \=\ 0/method\ \=\ 1/' */tmp0001.ctl
grep 'method = 1' */tmp0001.ctl | wc -l # You should get as many as datasets you have
grep 'alpha' */tmp0001.ctl   # You should see `fix_alpha = 0` and `alpha = 0.5`
grep 'ncatG' */tmp0001.ctl   # You should see `ncatG = 5`
grep 'model' */tmp0001.ctl   # You should see `model = 4` (i.e., HKY model)
```

Now, we can run `BASEML`!

```sh
## Run from `00_BASEML`
# You should still be in this directory 
# but, if not, please change directories until
# you are there. Then, run the following
# commands.
home_dir=$( pwd )
for i in `seq 1 3`
do
cd $home_dir/$i
mkdir $home_dir/$i/baseml
cp $home_dir/$i/tmp0001* $home_dir/$i/baseml
cd  $home_dir/$i/baseml
rm  tmp0001.out
baseml *ctl > $home_dir/$i/baseml/log.txt
done
cd $home_dir
```

<details>
<summary><b>TIP FOR ANALYSES WITH LARGE DATASETS</b></summary>
<br>

<i>If you have a large dataset, you may want to use a job array to run the commands above: one task for each alignment block, which can take from hours to days depending on the site patterns of your alignment! You can [check the pipelines that we already provide as part of the tutorial about reproducible timetree inference](../01_PAML/00_BASEML/scripts/) if you need some examples!</i>

</details>
<br>

The branch lengths, the gradient, and the Hessian will be saved in output file `rst2`. What we now need to do is concatenate the content of the three `rst2` files we have just obtained (one for each alignment block) in a unique file: **the so-called `in.BV` file**. It is of utmost importance that we append them in the same order as the alignment blocks that were used to generate them appear in the partitioned input sequence file (i.e., the content of the `rst2` file output when reading the first alignment block will be the first block in the `in.BV`, the content of the `rst2` output file generated when reading the second alignment block will be second, etc.). We can do generate our `in.BV` file by using the following commands:

```sh
## Run from `00_BASEML`
# You should still be in this directory 
# but, if not, please change directories until
# you are there. Then, run the following
# commands.
home_dir=$( pwd )
for i in `seq 1 3`
do
printf "[[ Adding rst2 file for alignment block "$i "]]\n"
if [[ $i -eq 1 ]]
then
# If it is the first alignment block, create
# the `in.BV` file, do not append
cat $i/baseml/rst2 > $home_dir/in.BV
# After the first block is added, the rest
# of the content can be appended
printf "\n" >> $home_dir/in.BV
else
# For the rest of the alignment blocks, 
# just append
cat $i/baseml/rst2 >> $home_dir/in.BV
printf "\n" >> $home_dir/in.BV
fi
done
cd $home_dir
```

Now, we have everything we need to run `MCMCtree`!

### Running `MCMCtree`

We are going to run `MCMCtree` when sampling from the prior (i.e., no data are used, useful to check whether there are problems between the calibration densities specified by the user and the corresponding marginal densities inferred by the program) and the posterior (i.e., when our data are used!).

We will run 6 chains when sampling from the prior and 6 chains when sampling from the posterior:

```sh
## Run from `00_BASEML`.
# You should still be in this directory 
# but, if not, please change directories until
# you are there. Then, run the following
# commands.
cd ../
for i in `seq 1 6`
do
# Create file structure for the analyses we will carry out
mkdir -p 01_MCMCtree/{00_prior/CLK,01_posterior/{ILN,GBM}}/$i
printf "[[ Copy files for chain "$i" ]]\n"
# Copy files for prior
cp mcmctree_orig.ctl 01_MCMCtree/00_prior/CLK/$i/mcmctree_$i.ctl
# Modify the path to the input sequence file:
# the name of the partitioned alignment was there,
# but the path has changed. We are using a relative
# path!
sed -i 's/seqfile\ \=\ /seqfile = \.\.\/\.\.\/\.\.\/\.\.\/\.\.\/00\_inp\_data\//' 01_MCMCtree/00_prior/CLK/$i/mcmctree_$i.ctl
# Same as above, but now modify the path to the
# input tree file
sed -i 's/treefile\ \=\ /treefile\ \=\ \.\.\/\.\.\/\.\.\/\.\.\/\.\.\/00\_inp\_data\//' 01_MCMCtree/00_prior/CLK/$i/mcmctree_$i.ctl
# Modify option `usedata` so that `MCMCtree` samples
# from the prior (i.e., no data)
sed -i 's/usedata..*/usedata\ \=\ 0/' 01_MCMCtree/00_prior/CLK/$i/mcmctree_$i.ctl
# The analysis will always run faster when no `sigma_2` values are
# sampled, so we specify `clock = 1`. If you specified a
# relaxed-clock model, you would get estimates for the rate
# variation in the `mcmc.txt` file apart from `mu` values
# (one for each alignment block if there is a partitioned
# alignment)
sed -i 's/clock..*/clock\ \=\ 1/' 01_MCMCtree/00_prior/CLK/$i/mcmctree_$i.ctl
# Copy files for posterior (ILN)
cp mcmctree_orig.ctl 01_MCMCtree/01_posterior/ILN/$i/mcmctree_$i.ctl
# Same as above: modify the path to where the input sequence
# file is
sed -i 's/seqfile\ \=\ /seqfile = \.\.\/\.\.\/\.\.\/\.\.\/\.\.\/00\_inp\_data\//' 01_MCMCtree/01_posterior/ILN/$i/mcmctree_$i.ctl
# Same as above: modify the path to where the input tree
# file is
sed -i 's/treefile\ \=\ /treefile\ \=\ \.\.\/\.\.\/\.\.\/\.\.\/\.\.\/00\_inp\_data\//' 01_MCMCtree/01_posterior/ILN/$i/mcmctree_$i.ctl
# Specify `usedata = 2` and the path to the `in.BV` that
# we have just generated, which has the branch lengths,
# gradient, and Hessian estimated by `BASEML` that 
# `MCMCtree` will use to approximate the likelihood
# calculation during the MCMC
sed -i 's/usedata..*/usedata\ \=\ 2 \.\.\/\.\.\/\.\.\/\.\.\/00\_BASEML\/in\.BV/' 01_MCMCtree/01_posterior/ILN/$i/mcmctree_$i.ctl
# Set `clock = 2` so that the independent-rates
# log-normal model is enabled
sed -i 's/clock..*/clock\ \=\ 2/' 01_MCMCtree/01_posterior/ILN/$i/mcmctree_$i.ctl
# Copy files for posterior (GBM)
cp mcmctree_orig.ctl 01_MCMCtree/01_posterior/GBM/$i/mcmctree_$i.ctl
# Same as above: modify the path to where the input sequence
# file is
sed -i 's/seqfile\ \=\ /seqfile = \.\.\/\.\.\/\.\.\/\.\.\/\.\.\/00\_inp\_data\//' 01_MCMCtree/01_posterior/GBM/$i/mcmctree_$i.ctl
# Same as above: modify the path to where the input tree
# file is
sed -i 's/treefile\ \=\ /treefile\ \=\ \.\.\/\.\.\/\.\.\/\.\.\/\.\.\/00\_inp\_data\//' 01_MCMCtree/01_posterior/GBM/$i/mcmctree_$i.ctl
# Same as above, set `usedata = 2`
sed -i 's/usedata..*/usedata\ \=\ 2 \.\.\/\.\.\/\.\.\/\.\.\/00\_BASEML\/in\.BV/' 01_MCMCtree/01_posterior/GBM/$i/mcmctree_$i.ctl
# Set `clock = 3` so that the geometric Brownian 
# motion model (autocorrelated rates) is enabled
sed -i 's/clock..*/clock\ \=\ 3/' 01_MCMCtree/01_posterior/GBM/$i/mcmctree_$i.ctl
done
```

Now, we can run `MCMCtree`! You can either run it on your PC or decide whether you want to prepare a job array to run these analyses on a cluster (we will not discuss the latter in this quick start tutorial, but feel free to look at the job arrays and other pipelines we use throughout the [tutorial on reproducible timetree inference](../01_PAML/01_MCMCtree/scripts/)):

```sh
## Run from `01_PAML`
# You should still be in this directory 
# but, if not, please change directories until
# you are there. Then, run the following
# commands.
cd 01_MCMCtree
home_dir=$( pwd )
for i in `seq 1 6`
do
cd $home_dir/00_prior/CLK/$i
printf "[[ Running MCMCtree for chain "$i" | prior ]]\n"
# Run `MCMCtree` while you see the screen output but
# you also save it in a file called `log_mcmc$i"_prior.txt"`
# You will run this analysis when sampling from the prior,
# the posterior under GBM, and the posterior under ILN
mcmctree *ctl 2>&1 | tee log_mcmc$i"_prior.txt"
printf "[[ Running MCMCtree for chain "$i" | posterior ILN ]]\n"
cd $home_dir/01_posterior/ILN/$i
mcmctree *ctl 2>&1 | tee log_mcmc$i"_postILN.txt"
printf "[[ Running MCMCtree for chain "$i" | posterior GBM ]]\n"
cd $home_dir/01_posterior/GBM/$i
mcmctree *ctl 2>&1 | tee log_mcmc$i"_postGBM.txt"
done
cd $home_dir
# We can now create a file that will be handy for subsequent
# steps: node numbers as given by `MCMCtree` 
# Only one tree is enough as we have the same tree topology!
grep 'Species tree for FigTree' -A1 00_prior/CLK/1/out.txt | sed -n '2,2p' > 00_prior/node_tree.tree
```

> **IMPORTANT**: when analysing your own datasets, you should always first run `MCMCtree` when sampling from the prior, proceed to [carry out MCMC diagnostics as explained below](./README.md#prior), and make sure that there are no problems between the calibration densities you specified and the marginal densities `MCMCtree` inferred. If you observed serious discrepancies, you would need to go back to your calibrations and check whether you made a mistake or you need to adjust them until the marginal densities really represent your belief about the fossil record. Then, once everything looks alright, you can run `MCMCtree` when sampling from the posterior, and then [run the corresponding MCMC diagnostics as explained below](./README.md#posterior). Nevertheless, we are running `MCMCtree` both when sampling from the prior and the posterior so that we can have the output files ready for both MCMC diagnostics while (hopefully!) completing this tutorial on time for this session :) You will see that this workflow (i.e., `prior --> checks --> prior again if checks failed --> check again --> posterior if everything is fine`) is the one you shall follow when running [the tutorial on reproducible timetree inference](../README.md#tutorial-on-reproducibility).

### MCMC diagnostics

#### Prior

Now that we have the output files from the different MCMC runs in an organised file structure, we are ready to run MCMC diagnostics!

<details>
<summary><b>How did we generate the input calibration files?</b></summary>
<br>

<i>Unfortunately, we will not have time to go through these scripts but, when you finish the tutorial, feel free to check how to run [the in-house R script](01_PAML/scripts/Include_calibrations_MCMCtree.R), which uses file [`calibrations.txt`](00_inp_data/calibs/raw_calibs/calibrations.txt) to generate the rest of the intermediate output files that you shall see inside [`raw_calibs`](00_inp_data/calibs/raw_calibs/) and [`tree_display`](00_inp_data/tree_display/).

If you have run `MCMCtree` and generated output file [`node_tree.tree`](01_PAML/01_MCMCtree/00_prior/node_tree.tree), you can then run the second [in-house R script `Merge_node_labels.R`](01_PAML/scripts/Merge_node_labels.R), which will output the [`Calibnodes_mtcdnapri.csv` file](00_inp_data/calibs/inp_calibs/Calibnodes_mtcdnapri.csv) that you will need for the subsequent steps of the MCMC diagnostics!</i>

</details>
<br>

We are going to run the R script [`MCMC_diagnostics_prior.R`](01_PAML/scripts/MCMC_diagnostics_prior.R) and follow the detailed step-by-step instructions detailed in the script. In a nutshell, the protocol will be the following:

1. Load the `mcmc.txt` files generated after each run.
2. Generate a convergence plot with the unfiltered chains.
3. Find whether there are major differences between the time estimates sampled across the chains for the same node sin the 97.5% and the 2.5% quantiles. If so, flag and delete said chains.
4. If some chains have not passed the filters mentioned above, create an object with the chains that have passed the filters.
5. Generate a new convergence plot with those chains that passed the filters.
6. Calculate Rhat, tail-ESS, and bulk-ESS to check whether chain convergence has been reached with the chains that have passed filters.

The MCMC diagnostics did not find any of the chains problematic after running [our in-house R script `MCMC_diagnostics_prior.R`](01_PAML/scripts/MCMC_diagnostics_prior.R). Therefore, we used [our in-house bash script `Combine_MCMC.sh`](01_PAML/scripts/Combine_MCMC.sh) to concatenate all the `mcmc.txt` files for the 6 chains in a unique file.

```sh
## Run from `01_MCMCtree`
# You should still be in this directory 
# but, if not, please change directories until
# you are there. Then, run the following
# commands.
cd ../scripts
cp Combine_MCMC.sh ../01_MCMCtree/00_prior
# One argument taken: number of chains
cd ../01_MCMCtree/00_prior
## Variables needed
## arg1 --> name of directory where analyses have taken place (e.g., CLK, GBM, ILN)
## arg2 --> output dir: mcmc_files_CLK, mcmc_files_GBM, mcmc_files_ILN, etc.
## arg3 --> "`seq 1 36`", "1 2 5", etc. | depends on whether some chains were filtered out or not
## arg4 --> clock model used: ILN, GBM, CLK
## arg5 --> number of samples specified to collect in the control file to run `MCMCtree`
## arg6 --> 'Y' to generate a directory with files compatible with programs such as `Tracer` to visually
##          inspect traceplots, convergence plots, etc. 'N' otherwise
## arg7 --> if arg6 is 'Y', arg7 needs to have a name for the `mcmcf4traces_<name>` that will be
##          generated. If `arg6` is equal to 'N', then please write `N` too.
path_to_data=$( echo CLK )
num_dat=1
name_dat=( 'mtcdnapri' ) # if you had more dataset, you would add them here!
count=-1 #  counter will start at 0 in first iteration!
for i in `seq 1 $num_dat`
do
count=$(( count + 1 ))
# Run the in-house script to generate concatenated
# `mcmc.txt` file and the individual `mcmc.txt` files
# ready to visually inspect in e.g., Tracer
./Combine_MCMC.sh $path_to_data mcmc_files_${name_dat[count]}_CLK "`seq 1 6`" GBM 20000 Y ${name_dat[count]}_CLK
done
```

The script above will generate a directory called `mcmc_files_<name_dataset>_CLK` inside the `00_prior` directory, where the `mcmc.txt` with the concatenated samples will be saved. In addition, a directory called `mcmcf4traces_<namedataset>_CLK` will also be generated so that formatted MCMC files compatible with programs such as `Tracer` can be used to check for chain convergence. A template script to generate the `FigTree.tre` file with this `mcmc.txt` has been saved inside the [`dummy_ctl_files`](01_PAML/dummy_ctl_files) directory.

We will now create a dummy alignment with only 2 nucleotides to quickly run `MCMCtree` with option `print = -1`. This setting will basically (i) ignore all the settings regarding the evolutionary model and the MCMC, (ii) read the `mcmc.txt` file which path is set in option `mcmcfile`, (iii) and summarise all the samples in such file to generate a timetree. To create the dummy alignment, we will run the [in-house R script `Generate_dummy_aln.R`](01_PAML/scripts/Generate_dummy_aln.R). Once you run it, a new directory called `dummy_aln` will be created, which will contain the dummy alignment. Then, we are ready to run `MCMCtree` with option `print = -1`! Note that the `mcmc.txt` file will have all the samples collected by the chains that passed previous filters during MCMC diagnostics.

```sh
## Run from `00_prior`
# You should still be in this directory 
# but, if not, please change directories until
# you are there. Then, run the following
# commands.
name_dat=( 'mtcdnapri' ) # you would list more datasets here
num_dat=1
count=-1
for i in `seq 1 $num_dat`
do
printf "\n[[ ANALYSING DATASET "${name_dat[count]}" ]]\n\n"
count=$(( count + 1 ))
base_dir=$( pwd )
cd ../../dummy_ctl_files
# Get path to directory where the ctl file is
ctl_dir=$( pwd )
cd ../../00_inp_data/tree_display/
# Get path to directory where the tree file is
# and file name
tt_dir=$( pwd )
name_tt=`ls ${name_dat[count]}"_calib_MCMCtree.tree"`
# Go to directory with dummy sequence file,
# get directory name and file name
cd $ctl_dir
cd ../dummy_aln
aln_dir=$( pwd )
name_aln=`ls *aln`
cd $base_dir
# Go to directory where the concatenated
# `mcmc.txt` file is and start preparing
# the directory to run `MCMCtree` with
# option `print = -1`
cd mcmc_files_${name_dat[count]}_CLK
printf "[[ Generating tree file for concatenated \"mcmc.txt\"  ... ... ]]\n"
cp $ctl_dir/*ctl .
name_mcmc=`ls *mcmc.txt`
sed_aln=$( echo $aln_dir"/"$name_aln | sed 's/\//\\\//g' | sed 's/_/\\_/g' |  sed 's/\./\\\./g' )
sed_tt=$( echo $tt_dir"/"$name_tt | sed 's/\//\\\//g' | sed 's/_/\\_/g' |  sed 's/\./\\\./g' )
sed -i 's/MCMC/'${name_mcmc}'/' *ctl
sed -i -e 's/ALN/'${sed_aln}'/' *ctl
sed -i 's/TREE/'${sed_tt}'/' *ctl
# Run now `MCMCtree` after having modified
# the global vars according to the path to
# these files. Then, rename the output tree
# file so we can easily identify later which
# tree belongs to which dataset easier
mcmctree *ctl
mv FigTree.tre FigTree_${name_dat[count]}"_CLK_95HPD.tree"
cd $base_dir
done
```

We now have our timetree inferred with all the samples collected by all the chains that passed the filters during MCMC diagnostics (when sampling from the prior)! The next step is to plot the calibration densities VS the marginal densities to verify whether there are any serious clashes that may arise because of truncation or problems with the fossil calibrations used. We will use the [in-house R script `Check_priors_margVScalib.R`](01_PAML/scripts/Check_priors_margVScalib.R) to generate these plots.

Once this script has finished, you will see that a new directory `plots/margVScalib/<name_dataset>` will have been created. Inside this directory, you will find one directory for each individual dataset with individual plots for each node. In addition, all these plots have been merged into a unique document as well (note: some plots may be too small to see for each node, hence why we have generated individual plots).

Now, once the MCMC diagnostics have finished, you can extract the final data that you can use to write a manuscript as it follows (easy way of accessing important files without having to navigate the file structure!):

```sh
## Run from `00_prior`
# You should still be in this directory 
# but, if not, please change directories until
# you are there. Then, run the following
# commands.
cd ../../
mkdir sum_files_prior
cp -R 01_MCMCtree/00_prior/mcmc_files*CLK/*CLK*tree sum_files_prior/
cp -R 01_MCMCtree/00_prior/CLK/*/*all_mean*tsv sum_files_prior/
cp -R plots/ESS_and_chains_convergence/*prior*pdf sum_files_prior/
cp -R plots/margVScalib sum_files_prior/
```

#### Posterior

Now it is time to analyse the samples we collected when running `MCMCtree` with our data!

We will run the R script [`MCMC_diagnostics_posterior.R`](01_PAML/scripts/MCMC_diagnostics_posterior.R) and follow the detailed step-by-step instructions detailed in the script, which are essentially the same ones you used when analysing the samples collected when sampling from the prior. Given that no problems have been found with any of the chains we ran, we are ready to concatenate the parameter values sampled across the 6 independent chains we ran:

```sh
## Run from `01_PAML`
# You should still be in this directory 
# but, if not, please change directories until
# you are there. Then, run the following
# commands.
cd scripts
cp Combine_MCMC.sh ../01_MCMCtree/01_posterior
# One argument taken: number of chains
cd ../01_MCMCtree/01_posterior
## Variables needed
## arg1 --> name of directory where analyses have taken place (e.g., CLK, GBM, ILN)
## arg2 --> output dir: mcmc_files_CLK, mcmc_files_GBM, mcmc_files_ILN, etc.
## arg3 --> "`seq 1 36`", "1 2 5", etc. | depends on whether some chains were filtered out or not
## arg4 --> clock model used: ILN, GBM, CLK
## arg5 --> number of samples specified to collect in the control file to run `MCMCtree`
## arg6 --> 'Y' to generate a directory with files compatible with programs such as `Tracer` to visually
##          inspect traceplots, convergence plots, etc. 'N' otherwise
## arg7 --> if arg6 is 'Y', arg7 needs to have a name for the `mcmcf4traces_<name>` that will be
##          generated. If `arg6` is equal to 'N', then please write `N` too.
path_to_data_GBM=GBM
path_to_data_ILN=ILN
num_dat=1
name_dat=( 'mtcdnapri' ) # if you had more dataset, you would add them here!
count=-1 #  counter will start at 0 in first iteration!
for i in `seq 1 $num_dat`
do
count=$(( count + 1 ))
./Combine_MCMC.sh $path_to_data_GBM mcmc_files_${name_dat[count]}_GBM "`seq 1 6`" GBM 20000 Y ${name_dat[count]}_GBM
./Combine_MCMC.sh $path_to_data_ILN mcmc_files_${name_dat[count]}_ILN "`seq 1 6`" ILN 20000 Y ${name_dat[count]}_ILN
done
```

Once the scripts above have finished, directories called `mcmc_files_part_[GBM|ILN]` and `mcmcf4traces` will be created inside `01_posterior/`. To map the mean time estimates with the filtered chains, we need to copy a control file, the calibrated Newick tree, and the dummy alignment we previously generated when analysing the results when sampling from the prior:

```sh
## Run from `01_posterior`
# You should still be in this directory 
# but, if not, please change directories until
# you are there. Then, run the following
# commands.
name_dat=( 'mtcdnapri' )
num_dat=1
count=-1
for i in `seq 1 $num_dat`
do
printf "\n[[ ANALYSING DATASET "${name_dat[count]}" ]]\n\n"
count=$(( count + 1 ))
base_dir=$( pwd )
cd ../../dummy_ctl_files
# Get path to directory where the ctl file is
ctl_dir=$( pwd )
cd ../../00_inp_data/tree_display/
# Get path to directory where the tree file is
# and file name
tt_dir=$( pwd )
name_tt=`ls ${name_dat[count]}"_calib_MCMCtree.tree"`
# Go to directory with dummy sequence file,
# get directory name and file name
cd $ctl_dir
cd ../dummy_aln
aln_dir=$( pwd )
name_aln=`ls *aln`
cd $base_dir
# Go to directory where the concatenated
# `mcmc.txt` file is and start preparing
# the directory to run `MCMCtree` with
# option `print = -1`
cd mcmc_files_${name_dat[count]}"_GBM"
printf "[[ Generating tree file for concatenated \"mcmc.txt\" for GBM ... ... ]]\n"
cp $ctl_dir/*ctl .
name_mcmc=`ls *mcmc.txt`
sed_aln=$( echo $aln_dir"/"$name_aln | sed 's/\//\\\//g' | sed 's/_/\\_/g' |  sed 's/\./\\\./g' )
sed_tt=$( echo $tt_dir"/"$name_tt | sed 's/\//\\\//g' | sed 's/_/\\_/g' |  sed 's/\./\\\./g' )
sed -i 's/MCMC/'${name_mcmc}'/' *ctl
sed -i -e 's/ALN/'${sed_aln}'/' *ctl
sed -i 's/TREE/'${sed_tt}'/' *ctl
# Run now `MCMCtree` after having modified
# the global vars according to the path to
# these files. Then, rename the output tree
# file so we can easily identify later which
# tree belongs to which dataset easier
mcmctree *ctl
mv FigTree.tre FigTree_${name_dat[count]}"_GBM.tree"
printf "\n"
# Now, do the same but for the analyses under ILN
cd $base_dir/mcmc_files_${name_dat[count]}"_ILN"
printf "[[ Generating tree file for concatenated \"mcmc.txt\" in "$data"/"$i" for ILN ... ... ]]\n"
cp $ctl_dir/*ctl .
name_mcmc=`ls *mcmc.txt`
sed_aln=$( echo $aln_dir"/"$name_aln | sed 's/\//\\\//g' | sed 's/_/\\_/g' |  sed 's/\./\\\./g' )
sed_tt=$( echo $tt_dir"/"$name_tt | sed 's/\//\\\//g' | sed 's/_/\\_/g' |  sed 's/\./\\\./g' )
sed -i 's/MCMC/'${name_mcmc}'/' *ctl
sed -i -e 's/ALN/'${sed_aln}'/' *ctl
sed -i 's/TREE/'${sed_tt}'/' *ctl
mcmctree *ctl
mv FigTree.tre FigTree_${name_dat[count]}"_ILN.tree"
printf "\n"
cd $base_dir
done
```

Once you have your final timetrees estimated under the two different relaxed-clock models (yay!), we can run our [in-house R script](01_PAML/scripts/Check_priors_VS_posteriors.R) to compare various distributions: marginal densities, calibration densities, and posterior time densities! These plots can help assess how informative the data are and check whether there are serious differences between the marginal densities and the posterior time densities.

Lastly, you can extract the final data that we used to write our manuscript as it follows:

```sh
## Run from `01_posterior`
# You should still be in this directory 
# but, if not, please change directories until
# you are there. Then, run the following
# commands.
cd ../../
mkdir sum_files_post
cp -R 01_MCMCtree/01_posterior/mcmc_files_*/FigTree*tree sum_files_post/
cp -R 01_MCMCtree/01_posterior/*/*/*all_mean*tsv sum_files_post/
cp -R plots/priorVSpost*pdf sum_files_post/
cp -R plots/ESS_and_chains_convergence/*post*pdf sum_files_post/
```

And... This is the end of the tutorial! Time to discuss the results and ask questions :)

----
