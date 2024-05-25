#-------------------#
# CLEAN ENVIRONMENT #
#-------------------#
rm( list = ls( ) )

#-----------------------------------------------#
# LOAD PACKAGES, FUNCTIONS, AND SET ENVIRONMENT #
#-----------------------------------------------#
# This package lets you find automatically the path to a specific location
# in your file structure
# If you have not installed this package, you will need to install it. 
# You can uncomment the following line to do this:
#install.packages( "rstudioapi" )
library( rstudioapi ) 
# Get the path to current open R script and find main dir
path_to_file <- getActiveDocumentContext()$path
wd <- paste( dirname( path_to_file ), "/", sep = "" )
setwd( wd )
# This script and the `Functions.R` script are saved under
# a directory called `scripts`. We will remove this part from 
# the path to find our home directory
home_dir <- gsub( pattern = "scripts/", replacement = "", x = wd )
# Load main script with all functions required below
source( file = "../../../src/Functions.R" )

#-------------------------------------------------------------#
# DEFINE GLOBAL VARIABLES -- modify according to your dataset #
#-------------------------------------------------------------#
# First, we will define global variables that we will keep using throughout this
# script.

# 1. Label the file with calibrations. If you have tested different calibrations
# and have more than one file with the corresponding calibrations, give as 
# many labels as files you have.
dat <- "mtcdnapri"

# 2. Number of divergence times that have been estimated. One trick to find
# this out quickly is to subtract 1 to the number of species. In this case,
# there are 7 taxa (7), so the number of internal nodes
# is `n_taxa-=7-1=6`.
# Another way to verify this is by opening the `mcmc.txt` file and check the
# header. The first element after `Gen` will have the format of `t_nX`, where
# X will be an integer (i.e., 8). Subtract two to this number 
# (i.e., 8-2=6) and this will be your number of divergence times that are 
# parameters of the MCMC. Please modify the number below so it fits to the 
# dataset you are using. 
num_divt <- 6

# 3. Total number of samples that you collected after generating the
# final `mcmc.txt` files with those from the chains that passed the filters. 
# You can check these numbers in scripts `MCMC_diagnostics_posterior.R` and
# `MCMC_diagnostics_prior.R`. E.g., `sum_post_QC$<name_dataset>$total_samples`
# or `sum_prior_QC$<name_dataset>$total_samples`
# CLK: The number of samples is 120006; you need to specify one less
# GBM: The number of lines is 120006; so you need to specify one less
# ILN: The number of lines is 120006; so you need to specify one less
#
# NOTE: If you had more than one dataset, you would add another vector of three
# values with the samples for CLK, GBM, and ILN to create `def_samples`
# E.g. two datasts: c( c( 120005, 120005, 120005), c( 120005, 120005, 120005) )
def_samples <- c( 120005, 120005, 120005 )

# 4. Quantile percentage that you want to set By default, the variable below is 
# set to 0.975 so the 97.5% and 2.5% quantiles (i.e., 95%CI). If you want to
# change this, however, just modify the value.
perc <- 0.975

# 5. Number of columns in the `mcmc.txt` that are to be deleted as they do not 
# correspond to sample values for divergence times (i.e., the entries are not 
# names following the format `t_nX`). To figure out this number quickly, you 
# can open the `mcmc.txt` file, read the header, and count the number of `mu*`
# and `sigma2*` elements. Do not count the `lnL` value when looking at 
# `mcmc.txt` files generated when sampling from the posterior -- this is 
# automatically accounted for in the in-house R functions that you will 
# subsequently use. E.g., assuming an MCMC ran under a relaxed-clock model with  
# no partitions, we would see `mu` and `sigma2` columns. Therefore, the variable  
# would be set to `delcol_post <- 2`. Please modify the value/s below 
# (depending on having one or more datasets) according to the `mcmc.txt` file
# generated when sampling from the posterior (`delcol_obj`). When running
# from the prior and `clock = 1`, you will only see `mu*` columns but, if you
# ran it with options `clock = 2` or `clock = 3`, you shall also see `sigma2*`
# columns.
##> NOTE: If you ran `MCMCtree` with `clock = 2` or `clock = 3` when
##> sampling from the prior, you will also need to count the `sigma2*`
##> columns! We ran `clock = 1` so that the analyses ran quicker, and thus
##> we only have `mu*` columns.
delcol_obj <- rep( c( 1, 6, 6 ), length(dat) )
                   # prior (1 col), postGBM (6 cols), postILN (6 cols)

# 6. Path to the directory where the concatenated `mcmc.txt` file has been 
# generated. Note that, if you have run more than one chain in `MCMCtree` for
# each hypothesis tested, you are expected to have generated a concatenated 
# `mcmc.txt` file with the bash script `Combine_MCMC_prior.sh` or any similar 
# approaches.
num_dirs  <- length( dat )*3 # Keep "3" if three models tested: CLK, GBM, ILN
paths_dat <- vector( mode = "character", length = num_dirs )
start <- stop <- 0
for( i in 1:length( dat ) ){
  start <- stop + 1
  stop  <- stop + 3
  #cat( "start: ", start, "\nstop:", stop, "\n" )
  paths_dat[start:stop] <- c( paste( home_dir,
                                     "01_MCMCtree/00_prior/mcmc_files_",
                                     dat[i], "_CLK", sep = "" ),
                              paste( home_dir,
                                     "01_MCMCtree/01_posterior/mcmc_files_",
                                     dat[i], "_GBM", sep = "" ),
                              paste( home_dir,
                                     "01_MCMCtree/01_posterior/mcmc_files_",
                                     dat[i], "_ILN", sep = "" ) )
}

# 7. Load a semicolon-separated file with info about calibrated nodes. Note that
# this file is output by script `Merge_node_labels.R`. A summary of its content
# in case you are to generate your own input files:
#
# Each column needs to be separated with semicolons and an extra blank line
# after the last row with calibration information needs to be added. If the
# extra blank is not added, R will complain and will not load the file!
# If you add a header, please make sure you name the column elements as 
# `Calib;node;Prior`. If not, the R function below will deal with the header,
# but make sure you set `head_avail = FALSE` when running `read_calib_f` 
# function below. An example of the content of this file is given below:
#
# ```
# Calib;node;Prior
# ex_n5;5;ST(5.8300,0.0590,0.1120,109.1240)
# ex_n7;7;B(4.1200,4.5200,0.0250,0.0250)
#
# ```
#
# The first column will have the name of the calibration/s that can help you
# identify which node belongs to which calibration. The second column is the
# number given to this node by`MCMCtree` (this information is automatically
# found when you run the script `Merge_node_labels.R`, otherwise you will need
# to keep checking the output file `node_trees.tree` to figure out which node
# is which). The third column is the calibration used for that node in
# `MCMCtree` format.
# 
# [[ NOTES ABOUT ALLOWED CALIBRATION FORMATS]]
#
# Soft-bound calibrations: 
#  E.g.1: A calibration with a minimum of 0.6 and a maximum of 0.8 would with  
#         the default tail probabilities would have the following equivalent 
#         formats:
#         >> B(0.6,0.8) | B(0.6,0.8,0.025,0.025)
#  E.g.2: A calibration with a minimum of 0.6 and a maximum of 0.8 would with  
#         the pL=0.001 and pU=0.025 would have the following format. Note that, 
#         whenever you want to modify either pL or pU, you need to write down 
#         the four  parameters in the format of "B(min,max,pL,pU)":
#         >> B(0.6,0.8,0.001,0.025)
#
# Lower-bound calibrations: 
#  E.g.1: A calibration with a minimum of 0.6 and the default parameters for
#         p = 0.1, c = 1, pL = 0.025:
#         >> L(0.6) | L(0.6,0.1,1,0.025)
#  E.g.2: A calibration with a hard minimum at 0.6, and so pL = 1e-300. 
#         Note that, whenever you want to modify either pL or pU, you need to  
#         write down the four parameters in the format of "L(min,p,c,pL)":
#         >> L(0.6,0.1,1,1e-300)
#
# Upper-bound calibrations: 
#  E.g.1: A calibration with a maximum of 0.8 and the default parameters for
#         pU = 0.025:
#         >> U(0.8) | U(0.8,0.025)
#  E.g.2: A calibration with a hard maximum at 0.8, and so pU = 1e-300. 
#         Note that, if you want to modify pU, you need to write down the two
#         parameters in the format of "U(max,pU)":
#         >> U(0.8,1e-300)
#
# ST distributions: 
#  The format accepted has four parameters: xi (location, mean root age), 
#  omega (scale), alpha (shape), nu (df). Accepted format: 
#  >> ST(5.8300,0.0590,0.1120,109.1240)
#
# SN distributions: 
#  The format accepted has three parameters: xi (location, mean root age), 
#  omega (scale), alpha (shape). Accepted format: 
#  >> SN(5.8300,0.0590,0.1120)  
#
#
# The next command executes the `read_calib_f` in-house function, which reads
# your input files (semicolon-separated files). The path to this directory is 
# what the argument `main_dir` needs. The argument `f_names` requires the name 
# of the file/s that you have used. Argument `dat` requires the same global 
# object that you have created at the beginning of the script.
dat_ff <- list.files( path = "../../00_inp_data/calibs/inp_calibs/",
                      pattern = "*csv",
                      full.names = FALSE )
##> CHECK
# If there was a `margVScalib` generated, do not use it! If not, this means
# that your `csv` file has already the required format to proceed
is_margcsv <- grep( pattern = "margVScalib", x = dat_ff )
if( length( is_margcsv ) > 0 ){
  dat_ff <- dat_ff[-is_margcsv]
}
##> END CHECK
calib_nodes <- vector( mode = "list", length = length(dat) )
names( calib_nodes ) <- dat
for( i in 1:length(dat_ff) ){
  calib_nodes[[ i ]] <- read_calib_f( main_dir = paste( "../../00_inp_data/calibs/inp_calibs/",
                                                        sep = "" ),
                                      f_names = dat_ff[i],
                                      dat = dat[i], head_avail = TRUE )
}

#-----------#
# LOAD DATA #
#-----------#
# Load mcmc files for all datasets
mcmc_obj <- vector( "list", num_dirs )
count <- 0
for( i in dat ){
  count <- count + 1
  names( mcmc_obj )[count] <- paste( "CLK_", i, sep = "" )
  count <- count + 1
  names( mcmc_obj )[count] <- paste( "GBM_", i, sep = "" )
  count <- count + 1
  names( mcmc_obj )[count] <- paste( "ILN_", i, sep = "" )
}
prior <- rep( c( TRUE, FALSE, FALSE ), length( dat ) )
count <- 0
for( i in c( paths_dat ) ){
  count <- count + 1
  cat( "[[ Parsing file for dataset", names( mcmc_obj )[count], " ]]\n" )
  mcmc_obj[[count]] <- load_dat( mcmc = paste( i, "/mcmc.txt", sep = "" ),
                                 delcol = delcol_obj[count], perc = perc,
                                 def_samples = def_samples[count],
                                 prior = prior[count] )
}

#---------------------------#
# PLOTS: prior VS posterior #
#---------------------------#
# Plot calibrated nodes
plot_priorVSpost( dat = dat, calib_nodes = calib_nodes,
                  mcmc_obj = mcmc_obj, home_dir = home_dir )
  



