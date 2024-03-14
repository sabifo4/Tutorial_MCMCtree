# `BASEML` analysis

## 1. Pick rate prior

We will use a vague gamma distribution centered on a mean evolutionary rate estimated by considering the tree height (molecular distance in substitutions per site) and the mean age for the root of our phylogeny (time unit = 100 Mya). As [our rooted phylogeny](../../00_data_formatting/00_raw_data/trees/) has information about the branch lengths, we can use [our R in-house script](scripts/calculate_rateprior.R) to calculate the tree heights. We also have a calibration to constrain the root age, which we will use as an approximate age for the root of our phylogeny to estimate the mean evolutionary rate.

By setting a vague shape ($\alpha=2$) for the gamma distribution that we will use as a rate prior, we can account for the uncertainty surrounding our mean rate estimate. If we had more knowledge on the mean rate, however, we should use a narrower prior with a larger $\alpha$ that better represents our prior information.

Now, we have all the information we need to calculate the $\beta$ parameters for the gamma distributions. We have written the [R script `calculate_rateprior.R`](scripts/calculate_rateprior.R) to carry out all the tasks mentioned above. You can open this file in RStudio to find out the value of $\beta$ and plot the final prior on the rates. A summary of what you will find in the script is described below:

```text
First, we know that the molecular distance (tree height, distance from the root to present time) is equal to the mean evolutionary rate (in substitutions per site per time unit) times the age of the divergence time at the root (in time unit, which we can define later). If we have estimated our phylogeny, and therefore have estimated the value of the branch lengths, we will be able to calculate the tree height. The units of the tree height will be the following:

tree_height = rate * root_age --> units_tree_height = subst/site/time * time = subst/site

There are various functions we can use to calculate the tree heigt. We have chosen the R function `phytools::nodeHeights`. The maximum height calculated by this function corresponds to the length from the root to the heighest tip.

After calculating the tree height of our phylogeny (in subst/site) and considering the age of the root based on the fossil record or geological events (time unit = 1 Ga = 100 Mya = 1e9 years), we can get a rough estimate for the mean rate:

Time unit = 100 Mya --> mean_rate = tree_height / root_age = (subst/site) / (Mya) = subst/site/Mya (time unit = 100 Mya) 

We also know that the mean of the gamma distribution that we want to use as rate prior is our parameter of interest: the mean evolutionary rate. Therefore:

mean_Gamma = mean_rate = alpha / beta 
Time unit = 100 Mya: mean_rate = alpha / beta --> beta = alpha / mean_rate = 2 / mean_rate

The calibrated tree needs to incorporate the age constraints in the same time unit used to infer the mean evolutionary rate and establish the rate prior (i.e., do not forget to scale the calibrations accordingly if needed!). 
```

If you run the [R script `calculate_rateprior.R`](scripts/calculate_rateprior.R), you will see how all the steps described above take place and a new PDF file with the prior distribution to be used will be generated in a new directory called `out_RData`.

We have then updated our [template control files](control_files) with the $\alpha$ and $\beta$ parameters (as defined using the R script above) for the gamma distributions as rate priors. Note that several options in this control file will be subsequently modified to fit the analysis with this dataset (i.e., you will see some options that have flags in capital letters, which will be replaced with the correct value for said option). Given how shallow this tree is, the clock may be expected to be seriously violated, and thus we have fixed a mean for the `sigma2` parameter (i.e., variation in the clock) as 0.1 using a gamma prior with $\alpha=1$ and $\beta=10$: `sigma2_gamma 1 10​`.

## 2. Set up the file structure

Before running `MCMCtree` using the approximate likelihood calculation to speed up timetree inference, we first need to calculate the vector of branch lengths, the gradient (vector), and the Hessian (matrix). We will use `BASEML` for this purpose as our dataset is an amino acid alignment.

The file structure we will use is the following:

```text
example_dating/
  |
  |- alignments/
  |    |- X/ # Directory for alignment X -- have as many directories as alignments
  |       
  |- control_files/ # Pre-defined control file with flags to be later replaced with specific settings
  |
  |- Hessian/
  |    |- X # Directory for alignment X -- have as many directories as alignments
  |          
  |- pipelines_Hessian # Directory where the pipeline to run `BASEML` will be executed
  |
  |- scripts # Scripts used to prepare control files to run `BASEML
  |
  |- trees
      |- calibrated   # Directory with the calibrated tree for `MCMCtree`
      |- uncalibrated # Directory with the uncalibrated tree for `BASEML`
```

To create the `example_dating` file structure, we run the following commands:

```sh
# Run the following commands from 
# directory `00_BASEML`
mkdir -p tmp/example_dating
cd tmp/example_dating 
num_aln=1
for i in `seq 1 $num_aln`
do
mkdir -p alignments/$i
mkdir -p Hessian/$i/prepare_baseml
mkdir -p trees/{uncalibrated/$i,calibrated/$i}
mkdir -p control_files/$i
done
mkdir -p pipelines_Hessian
mkdir scripts
```

Once the file structure is created, we can now populate it with the input files we have generated some minutes ago: alignment files, tree files, and control files. We will now relocate the location of this `example_dating` directory so that our pipelines work:

```sh
# Run from `tmp/example_dating`
# Copy alignment
cp ../../../../00_data_formatting/01_inp_data/example_aln.phy alignments/1/
# Now, transfer calibrated tree
cp ../../../../00_data_formatting/01_inp_data/example_calib_MCMCtree.tree trees/calibrated/1
# Transfer uncalibrated tree
cp ../../../../00_data_formatting/01_inp_data/tree_example_uncalib.tree trees/uncalibrated/1
# Next, copy control file
cp ../../control_files/prepbaseml.ctl control_files/1
# Last, copy the in-house bash scripts with our pipeline
cp ../../scripts/*sh scripts/
cd ../
mv example_dating ../../../
cd ../
rm -r tmp
```

Now, we need to generate other input files to estimate the Hessian and the gradient: the input control files for `BASEML`. To do this in a reproducible manner, you can use the [script `generate_prepbaseml.sh`](scripts/generate_prepbaseml.sh), which you can find in the [`01_PAML/00_BASEML/scripts`](01_PAML/00_BASEML/scripts) and which you should have in your `example_dating` directory. Now,go to the `example_dating` directory and run the next code snippet, where you will execute this script. Specifically, the [`generate_prepbaseml.sh` script](scripts/generate_prepbaseml.sh) needs one argument: the directory name where the alignments are saved: `1`, `2`, and `3` in our case!

```sh
# Run from `example_dating/scripts`
# Please change directories until
# you are there. Then, run the following
# commands.
chmod 775 *sh
# In this case, there are three alignments, so
# we can execute our script within a loop
num_aln=1
for i in `seq 1 $num_aln`
do
./generate_prepbaseml.sh $i
done
```

To make sure that all the paths have been properly extracted, you can run the following code snippet:

```sh
# Run from `example_dating/Hessian`
# Please change directories until
# you are there. Then, run the following
# commands.
grep 'seqfile' */prepare_baseml/*ctl
grep 'treefile' */prepare_baseml/*ctl
```

## 3. Run `BASEML`

### Preparing input files

Now that we have the input files (alignment and tree files) and the instructions to run `BASEML` (control file) in `example_dating`, we will be manually running `MCMCtree` inside each `prepare_baseml` directory (see file structure above) in a special mode that launches `BASEML` for the sole purpose we want: to infer the vectors and matrix required to approximate the likelihood calculation.

```sh
# Run `MCMCtree` from
# `example_dating/Hessian/1/prepare_baseml`.
# Please change directories until
# you are in there.
# The first command to change directories 
# will work if you are still in 
# `main/Hessian`, otherwise ignore and 
# move to such directory with the command
# that best suits your current directory.
# If you had more than one alignment, you
# could write a `for` loop or access each
# dir individually.
dir=1
cd $dir/prepare_baseml
mcmctree prepbaseml*ctl # You may have other aliases to run `MCMCtree`, 
                        # so run this command accordingly!
```

First, you will see that `MCMCtree` starts parsing the first locus. Then, you will see something like the following printed on your screen (some sections may change depending on the PAML version you have installed on your cluster!):

```text
*** Locus 1 ***
running baseml tmp0001.ctl
BASEML in paml version 4.10.7, June 2023
ns = 4          ls = 48
Reading sequences, sequential format..
Reading seq # 4: sp4
Sequences read..

48 site patterns read, 609 sites
Counting frequencies..

       72 bytes for distance
     4608 bytes for conP
     1536 bytes for fhK
  8000000 bytes for space
```

As soon as you see the last line, you will see that various `tmp000X*` files will have been created, and hence you can stop this run by typing `ctrl+C` on the terminal that you have used to run such command -- the dummy alignment is so short that you will not have time to even do this, but this trick is useful when having large phylogenomic datasets!

Once you have done this, you can check that the control file you will later need has been created:

```sh
# Run from the `example_dating/Hessian`
# Please change directories until
# you are there. Then, run the following
# commands.
grep 'seqfile' */*/tmp0001.ctl | wc -l # You should get as many datasets as you have, in this case 1
```

Note that, when we ran the commands above, we were not interested in running `BASEML` or `MCMCtree`. We just wanted to execute `MCMCtree` with option `usedata = 3` so that it generates the `tmp000*` files that `BASEML` will later need to estimate the branch lengths, the gradient, and the Hessian. We do this analysis in two steps given that there may be restrictions in the HPC you will be using in the future (or are using at the moment) that do not allow to run `BASEML` + `MCMCtree` in one unique job within a reasonable amount of time (i.e., this applies to genomic datasets in case you wanted to scale your analyses with larger datasets, you could see how quickly this runs with a dummy dataset!). In addition, we want to modify some settings in the control file that is automatically generated when enabling `usedata = 3` so that they match what we want to do for our inference. In a nutshell, this is what you will be doing:

1. Run `MCMCtree` to generate the `tmp000*` files.
2. Modify the `tmp0001.ctl` file according to the settings we want to enable to analyse our dataset with `BASEML`.
3. Run `BASEML` using the `tmp000*` files so that it estimates the branch lengths, the gradient, and the Hessian and saves them in a file called `rst2`.
4. Generate the final `in.BV` file for our dataset, which will be later used by `MCMCtree`.

Once all `tmp000*` files are generated for all alignments, we need to make sure that the correct evolutionary model has been enabled (i.e., `model = 4`, `ncatG=4` for HKY85+G4) and that option `method = 1` is enabled, which will speed up the computation of the Hessian and the gradient. We can run the next code snippet to very that the four requirements aforementioned are met:

```sh
# Run from the `example_dating/Hessian`
# Please change directories until
# you are there. Then, run the following
# commands.
sed -i 's/method\ \=\ 0/method\ \=\ 1/' */*/tmp0001.ctl
grep 'method = 1' */*/tmp0001.ctl | wc -l # You should get as many as datasets you have
grep 'alpha' */*/tmp0001.ctl   # You should see `fix_alpha = 0` and `alpha = 0.5`
grep 'ncatG' */*/tmp0001.ctl   # You should see `ncatG = 4`
grep 'model' */*/tmp0001.ctl   # You should see `model = 3` (i.e., empirical+F model)
```

### Executing `BASEML`

We can now run `BASEML` given that we have the control file ready as well as all the required input files!

We have created a template bash script with flags (i.e., see script  `pipeline_Hessian_BASEML_template_PC.sh` in the [`scripts` directory](01_PAML/00_BASEML/scripts)), which will be replaced with the appropriate values by another bash script (i.e.,`generate_job_BASEML_PC.sh`, also saved in the [`scripts` directory](01_PAML/00_BASEML/scripts)). Please note that the second bash script will edit the template bash script according to the data alignment/s that will be analysed. We had already copied these scripts to the `example_dating` directory when setting our file structure. Therefore, we just need to execute the following code snippet there:

```sh
# Run from `example_dating` dir.
# Please change directories until
# you are there. Then, run the following
# commands.
home_dir=$( pwd )
cd scripts
chmod 775 *sh
num_aln=1
# Arg1: Number of alignments
# Arg2: Path to the pipeline directory
# Arg3: Name of the working directory (i.e., `example_dating` in this analysis)
# Arg4: Name of the executable file for BASEML. E.g., `baseml4.10.7`, `baseml`, etc.
# Arg5: Boolean, PAML exported to the path? `Y` or `N`.
#       If `N`, the executable file will be required to be in the home dirctory,
#       i.e., directory which name you type as `Arg3`.
./generate_job_BASEML_PC.sh $num_aln $home_dir/pipelines_Hessian example_dating baseml Y
```

Next, we will go to the `pipelines_Hessian` directory and run the script that will have been generated using the commands above:

```sh
# Run from `example_dating/pipelines_Hessian`.
# Please change directories until
# you are there. Then, run the following
# commands.
#
# If you list the content of this directory,
# you will see the pipeline you will need 
# to execute in a bash script called
# `pipeline_Hessian.sh`
ll *
# Now, execute this bash script
chmod 775 *sh
./pipeline_Hessian.sh & # Include the `&` to run this job in the background!
```

Once `BASEML` finishes, we are ready to generate the `in.BV` file that we will later use when running `MCMCtree` to approximate the likelihood calculation:

```sh
# Run from dir `example_dating/Hessian/`
# Please change directories until
# you are there. Then, run the following
# commands.
num_aln=1
for i in `seq 1 $num_aln`
do
printf "\nGenerating in.BV files for dir "$i" ... ...\n\n"
cp $i/rst2 $i/in.BV
done
```

We can now proceed to timetree inference with `MCMCtree`! [You can click this link to move to the next `README.md` file](../01_MCMCtree/README.md)!
