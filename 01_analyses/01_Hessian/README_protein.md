# Calculate the gradient and the Hessian with `CODEML`

>> **IMPORTANT NOTE**: The dummy alignment with protein data is yet to be included in this repository. The instructions below, however, can be followed to run analyses with protein data.

Before running `MCMCtree`, we need to calculate the gradient and the Hessian so we can use the approximate likelihood during timetree inference to save computational time ([dos Reis and Yang, 2011](https://academic.oup.com/mbe/article/28/7/2161/1051613)).

## 1. Pick rate prior

We will use a vague gamma distribution for the dataset considering the tree height (molecular distance in substitutions per site) and the divergence time at the root of the phylogeny (in time unit). As the [tree file in NEXUS](../../00_data/00_raw_data/tree_ML.nexus) format has information about the branch lengths, we can load this file in `R` to estimate the tree height. We also have defined the root calibration in the calibrated tree file that you just generated, which is a rough idea of the age of the root of the phylogeny based on the fossil record. This calibration suggests that the mean root age is 583 Myr (i.e., the first parameter of the ST distribution used to calibrate the root is $\alpha=5.83$, and so we could use either 583 Myr [time unit = 1 Myr] or 5.83 in 100 Myr time unit).

By setting a vague shape ($\alpha=2$) for the gamma distribution, we can account for the uncertainty on the mean rate. If we had more knowledge on the mean rate, however, we should use a narrower prior with a larger $\alpha$ that better represents our prior information.

Now, we have all the information we need to calculate the $\beta$ parameter for the Gamma distribution that will be the prior on the estimated mean evolutionary rate. We have written the [R script `calculate_rateprior.R`](scripts/calculate_rateprior.R) to carry out all the tasks mentioned above. You can open this file in RStudio to find out the value of $\beta$ and plot the final prior on the rates. A summary of what you will find in the script is described below:

```text
First, we know that the molecular distance (tree height, distance from the root to present time) is equal to the mean evolutionary rate (in substitutions per site per year) times the age of the divergence time at the root (in time unit, which we can define later). If we have estimated our phylogeny, and therefore have estimated the branch lengths, we will be able to estimate the tree height. The units of the tree height will be the following:

tree_height = rate * root_age --> units_tree_height = subst/site/y * y = subst/site

One way of estimating the tree height is by using the R function `phytools::nodeHeights`. The maximum height calculated by this function corresponds to the length from the root to the heighest tip. 

After estimating the tree height of our phylogeny (in subst/site) and considering the age of the root based on fossils (time unit = 1 Myr), we can get a rough estimate of the mean rate depending. We will calculate the mean rate using two different time units:

(a) Time unit = 1 Myr (mean root age in Myr)    --> mean_rate_1Myr = tree_height / root_age = (subst/site) / (Myr) = subst/site per time unit (time unit = 1 Myr = 10^6 years) --> (subst/site)/10^6 years.
(b) Time unit = 100 Myr (mean root age in 100 Myr) --> mean_rate_100Myr = tree_height / root_age = (subst/site) / (Myr) =  subst/site per time unit (time unit = 100 Myr = 10^8 years) --> (subst/site)/10^8 years

We also know that the mean of the gamma distribution is our parameter of interest: the mean evolutionary rate. Therefore:

mean_G = mean_rate = alpha / beta 

According to the two cases stated above:

(a) Time unit = 1 Myr: mean_rate_1Myr = alpha / beta --> beta = alpha / mean_rate = 2/mean_rate_1Myr  
(b) Time unit = 100 Myr: mean_rate_100Myr = alpha / beta --> beta = alpha / mean_rate = 2/mean_rate_100Myr

The beta parameter when time unit = 1 Myr is normally too large, so we shall focus on the mean rate calculated when time unit = 100 Myr. In that way, the calibrated tree should be in this time unit too (i.e., do not forget to scale the calibrations accordingly!). 
```

If you run the [R script `calculate_rateprior.R`](scripts/calculate_rateprior.R), you will see how all the steps described above take place and a new PDF file with the prior distribution to be used will be generated in a new directory called `out_RData`.

As part of this tutorial, we have included a [template control file](control_files/prepcodeml.ctl) with the $\alpha$ and $\beta$ parameters (as defined using the R script above) for the gamma distribution as a prior on the rates. Note that several options will be subsequently modified to fit the analysis with this dataset (i.e., you will see some options that have flags in capital letters, which will be replaced with the correct value for said option).

Please note that, if you are adapting this tutorial to analyse a different dataset, you should change the options of the control file to fit your dataset. E.g.:

* Adjust `ndata` to the number of alignment blocks (i.e., amount of data partitions) in your alignment file.
* Adjust the prior on the rate, `rgene_gamma`, according to the value of $\alpha$ and $\beta$ that you have calculated following the steps above.
* Adjust the prior on $\sigma^2$ according to your dataset. Large $\sigma^2$ (e.g., 0.2 or 0.5) means serios violations of the clock and small values of $\sigma^2$ (e.g., 0.01) mean that the clock almost holds. For shallow trees, the clock may not be very wrong and a prior with mean 0.05 can be specified, such as `sigma2_gamma 2 40​`. For deep phylogenies, the clock may be expected to be seriously violated and you can specify a larger mean as a prior such as 0.2 using something like `sigma2_gamma 2 10​`.

## 2. Set up the file structure

Before running `MCMCtree` using the approximate likelihood calculation to speed up timetree inference, we first need to calculate the gradient and the Hessian. We will use `CODEML` for this purpose!

The file structure we will use is the following:

```text
main_protein/
  |
  |- alignments/
  |    |- X/ # Directory for alignment X -- have as many directories as alignments
  |       
  |- control_files/ # Pre-defined control file with flags to be later replaced with specific settings
  |
  |- Hessian/
  |    |- X # Directory for alignment X -- have as many directories as alignments
  |          
  |- pipelines_Hessian # Directory where the pipeline to run `CODEML` will be executed
  |
  |- scripts # Scripts used to prepare control files to run `CODEML
  |
  |- trees
      |- calibrated   # Directory with the calibrated tree for `MCMCtree`
      |- uncalibrated # Directory with the uncalibrated tree for `CODEML`
```

> **IMPORTANT NOTE:** The example dataset used in this tutorial is quite small and you should be able to run the analyses described in this tutorial quite quick on your own PC. Nevertheless, if you are to use large phylogenomic datasets, you should use an HPC. We are currently finalising an adaptation of this tutorial to be used in the HPC, but will be available soon. For the aim of this workshop, however, you can just run all the analyses on your PC, so you can continue reading below!

To create the `main_protein` file structure, please run the following commands:

```sh
# Run the following commands from the
# directory `01_analyses`
cd ../
mkdir main_protein 
cd main_protein
# Important to set `num_aln` in case 
# you decide to use this tutorial later to 
# analse more than one alignment
num_aln=1
for i in `seq 1 $num_aln`
do
mkdir -p alignments/$i
mkdir -p Hessian/$i/prepare_codeml
mkdir -p pipelines_Hessian
mkdir control_files
mkdir -p trees/{uncalibrated,calibrated}
mkdir scripts
done
```

Once the file structure is created, we can now populate it with the input files we have generated some minutes ago: alignment file, tree files, and control file. E.g.:

```sh
# Run from `00_data/01_inp_data`
# Please change directories until 
# you are there. Then, run the 
# following commands.
#
# Copy input alignment and tree files 
cp *phy ../../main_protein/alignments/1
cp *_uncalib.tree ../../main_protein/trees/uncalibrated/
cp *_calib.tree ../../main_protein/trees/calibrated/
# Change to where the control file is and transfer the 
# template control file
cd ../../01_analyses/01_Hessian/control_files
cp prepcodeml.ctl ../../../main_protein/control_files
cp lg.dat ../../../main_protein/control_files
```

While the alignment file, the tree files, and the control file (with the already correct prior on the rates, MCMC settings, and flags to later replace with the correct values) have already been generated and only need to be copied onto their corresponding directories as shown above, we need to generate other input files to estimate the Hessian and the gradient: the input control files for `CODEML`.

To do this in a reproducible manner, you can use the [script `generate_prepcodeml.sh`](scripts/generate_prepcodeml.sh), which you can find in the [`01_analyses/01_Hessian/scripts` directory](scripts). You should copy this bash script in the `main_protein/scripts` directory previously created:

```sh
# Run from the `01_analyses/01_Hessian/scripts`
# dir on your local PC. Please change
# directories until you are there. Then,
# run the following commands.
cp generate_prepcodeml.sh ../../../main_protein/scripts
```

The [`generate_prepcodeml.sh` script](scripts/generate_prepcodeml.sh) needs one argument: the amount of alignment files used. As we are just using one alignment file, we will use `1` as the argument:

```sh
# Run from `main_protein/scripts` on your local
# PC. Please change directories until
# you are there. Then, run the following
# commands.
chmod 775 *sh
# In this case, there is only one 
# alignment and hence the `for` loop
# does not seem adequate as it would go
# from `1` to `1` (i.e., no need). 
# Nevertheless, we have incorporated this
# code snippet in case you are to 
# reuse this tutorial with more alignments
# and `num_aln` increases
num_aln=1
for i in `seq 1 $num_aln`
do
./generate_prepcodeml.sh $i
done
```

To make sure that all the paths have been properly extracted, you can run the following code snippet:

```sh
# Run from `main_protein/Hessian` dir on your local
# PC. Please change directories until
# you are there. Then, run the following
# commands.
grep 'seqfile' */prepare_codeml/*ctl
grep 'treefile' */prepare_codeml/*ctl
grep 'aaRatefile' */prepare_codeml/*ctl
```

## 2. `CODEML`

### Preparing input files

Now that we have the input files (alignment and tree files) and the instructions to run `CODEML` (control file), we will be manually running `MCMCtree` inside each `prepare_codeml` directory (see file structure above) in a special mode that launches `CODEML` for the purpose want: calculating the gradient and the Hessian.

```sh
# Run `MCMCtree` from
# `main_protein/Hessian/1/prepare_codeml`
# dir on your local PC. 
#
# We execute `MCMCtree` by running
# the command `mcmctree` as this is
# the name our compiled version has.
# Nevertheless, please change this 
# command to the one you are using in
# case you have given another name 
# to the executable file that runs 
# `MCMCtree`.
#
# Please change directories until
# you are in `main_protein/Hessian/1/prepare_codeml`.
# If you have more than one alignment, 
# please run these two commands inside
# each `main_protein/Hessian/<num_aln>/prepare_codeml`.
# Do not run this in a `for` loop -- there are
# no issues with small datasets, but larger
# datasets need more attention and running a 
# `for` loop might make it impossible
cd 1/prepare_codeml
mcmctree prepcodeml.ctl
```

First, you will see that `MCMCtree` starts parsing the first locus. Then, you will see something like the following printed on your screen:

```text
*** Locus 1 ***
running codeml tmp0001.ctl

AAML in paml version 4.10.6, November 2022
ns = 4          ls = 48
Reading sequences, sequential format..
Reading seq # 4: sp4
Sequences read..

48 site patterns read, 609 sites
Counting frequencies..
```

As soon as you see the last line, you will see that the `tmp000X*` files will have been created, and hence you can stop this run by typing `ctrl+C`. Given that the example dataset is small and will finish quick, you will see that you may not even have time to stop the run as it will end before. Nevertheless, you need to do this with larger datasets. Therefore, please do not stop the run until the `tmp000X*` files are created in case you are using large datasets.

**NOTE**: Remember that we do not want to run `MCMCtree` now, we just want to obtain the `tmp000X*` files to run `CODEML` to estimate the gradient and the Hessian. You can check that all the needed files have been created by running the following command:

```sh
# Run from the `main_protein/Hessian` dir on your local
# PC. Please change directories until
# you are there. Then, run the following
# commands.
grep 'seqfile' */*/tmp0001.ctl | wc -l # You should get as many datasets as you have, in this case 1
```

In addition, you need to make sure that option `method = 1` is enabled, which will speed up the computation of the Hessian and the gradient:

```sh
# Run from the `main_protein/Hessian` dir on your local
# PC. Please change directories until
# you are there. Then, run the following
# commands.
sed -i 's/method\ \=\ 0/method\ \=\ 1/' */*/tmp0001.ctl
grep 'method = 1' */*/tmp0001.ctl | wc -l # You should get as many as datasets you have
```

### Running `CODEML`

Now that we have the control file ready to run `CODEML` as well as the required input files, we can run `CODEML`!

We have created a template bash script with flags
(i.e., `pipeline_Hessian_CODEML_template.sh` in the [`scripts` directory](scripts)), which will be replaced with the appropriate values by another bash script (`generate_job_CODEML.sh`, also saved in the [`scripts` directory](scripts)). Please note that the second bash script will edit the template bash script according to the data alignment/s that will be analysed. Now, we just need to copy them onto the `main_protein/scripts` directory:

```sh
# Run from `main_protein` dir on your local
# PC. Please change directories until
# you are there. Then, run the following
# commands.
home_dir=$( pwd )
cp ../01_analyses/01_Hessian/scripts/generate_job_CODEML.sh $home_dir/scripts
cp ../01_analyses/01_Hessian/scripts/pipeline_Hessian_CODEML_template.sh $home_dir/scripts
cd $home_dir/scripts
chmod 775 *sh
num_aln=1
for i in `seq 1 $num_aln`
do
# The first argument corresponds to the number of 
# gene alignments and the second to the path to 
# where the pipeline will be executed: `pipelines_Hessian`
./generate_job_CODEML.sh $num_aln $home_dir/pipelines_Hessian
done
```

Now, we just need to go to the `pipelines_Hessian` directory and run the script that will have been generated using the commands above:

```sh
# Run from `main_protein/pipelines_Hessian` dir on your local
# PC. Please change directories until
# you are there. Then, run the following
# commands.
#
# If you list the content of this directory,
# you will see the pipeline you will need 
# to execute in a bash script called
# `pipelines_CODEML.sh`
ls *
# Now, execute this bash script
chmod 775 *sh
./pipeline_CODEML.sh
```

Once `CODEML` finishes, we are ready to generate the `in.BV` files that we will later use when running `MCMCtree` using the approximate likelihood calculation:

```sh
# Run from dir `main_protein/Hessian` dir on your local
# PC. Please change directories until
# you are there. Then, run the following
# commands.
num_aln=1
for i in `seq 1 $num_aln`
do 
printf "\nGenerating in.BV files for dir "$i"  ... ...\n\n"
mv $i/rst2 $i/in.BV
done
```

You can now move on and [start the next tutorial to run `MCMCtree`](../02_MCMCtree/README_protein.md)!

----

Once the job finishes, if you want to count the constant sites in your alignment and the length of the alignment/s, you can run the following commands:

```sh
# Run from `main_protein/Hessian` dir on your local
# PC. Please change directories until
# you are there. Then run the following
# commands.
#
# NOTE: If you are running this tutorial with
# your own data, please remember to change
# `num_aln` and `num_taxa` accordingly.
num_aln=1
num_taxa=4
home_dir=$( pwd )
for i in `seq 1 $num_aln`
do
cd $i
printf "\nPrinting site patterns for dataset "$i":\n"
grep 'constant' *out
# Move to the alignments directory now
cd ../../alignments/$i/
printf "\nPrinting the length of dataset "$i":\n"
grep '^'$num_taxa'' *phy
cd $home_dir
done
```
