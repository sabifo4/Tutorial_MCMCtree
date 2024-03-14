#-------------------#
# CLEAN ENVIRONMENT #
#-------------------#
rm( list = ls( ) )

#---------------------------------------------------#
# SET WORKING DIRECTORY AND LOAD IN-HOUSE FUNCTIONS #
#---------------------------------------------------#
# Load colour-friendly package
# If you have not installed this package, you will need to install it. 
# You can uncomment the following line to do this:
#install.packages( "ColorBlindness" )
library( colorBlindness )
# This package lets you find automatically the path
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

#--------------------------------#
# DEFINE USER'S GLOBAL VARIABLES #
#--------------------------------#
# First, we will define global variables that we will keep using throughout this
# script.

# 1. Label the file with calibrations. If you have tested different calibrations
# and have more than one file with the corresponding calibrations, give as 
# many labels as files you have.
dat <- c( "example" )

# 2. Number of divergence times that have been estimated. One trick to find
# this out quickly is to subtract 1 to the number of species. In this case,
# there are 4 taxa, so the number of internal nodes
# is `n_taxa-=4-1=3`.
# Another way to verify this is by opening the `mcmc.txt` file and check the
# header. The first element after `Gen` will have the format of `t_nX`, where
# X will be an integer (i.e., 5). Subtract two to this number 
# (i.e., 5-2=3) and this will be your number of divergence times that are 
# parameters of the MCMC. Please modify the number below so it fits to the 
# dataset you are using. 
num_divt <- 3

# 3. Number of samples that you specified in the `MCMCtree` control file to 
# collect. NOTE: you may have not collect them all, but do not worry!
# The number of lines for CLK, GBM, and ILN was 120006, so you need to 
# specify one less
def_samples <- c( 120005, 120005, 120005 ) 

# 4. Quantile percentage that you want to set By default, the variable below is 
# set to 0.975 so the 97.5% and 2.5% quantiles. If you want to change this,
# however, just modify the value.
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
# would be set to `delcol = 2`. Please modify the values below according to your  
# the `mcmc.txt` file generated when sampling from the prior (`delcol_prior`) 
# dataset for and when sampling from the posterior (`delcol_posterior`).
delcol_obj <- rep( c( 1, 2, 2 ), length( dat ) ) # prior (1), posterior (2), posterior (2)

# Path to the directory where the concatenated `mcmc.txt` file has been 
# generated. Note that, if you have run more than one chain in `MCMCtree` for
# each hypothesis tested, you are expected to have generated a concatenated 
# `mcmc.txt` file with the bash script `Combine_MCMC_prior.sh` or any similar 
# approaches.
num_dirs  <- length( dat )*3
paths_dat <- vector( mode = "character", length = num_dirs )
start <- stop <- 0
for( i in 1:length( dat ) ){
  start <- stop + 1
  stop  <- stop + 3
  #cat( "start: ", start, "\nstop:", stop, "\n" )
  paths_dat[start:stop] <- c( paste( home_dir,
                                     "sum_analyses/00_prior/mcmc_files_",
                                     dat[i], "_CLK", sep = "" ),
                              paste( home_dir,
                                     "sum_analyses/01_posterior/mcmc_files_",
                                     dat[i], "_GBM", sep = "" ),
                              paste( home_dir,
                                     "sum_analyses/01_posterior/mcmc_files_",
                                     dat[i], "_ILN", sep = "" ) )
}

# Load semicolon-separated file/s with info about calibrated nodes. Note that
# each column needs to be separated with semicolons and an extra blank line
# after the last row with calibration information needs to be added (i.e., files
# need to have an extra blank line so R does not complain when reading them). 
# If you add a header, please make sure you name the column elements as 
# `Calib;node;Prior`. If not, the R function below will deal with the header. 
# An example of the format you need to follow to summarise the calibration info
# for each node is the following:
#
# ```
# Calib;node;Prior
# ex_n5;5;ST(5.8300,0.0590,0.1120,109.1240)
# ex_n7;7;B(4.1200,4.5200,0.0250,0.0250)
#
# ```
#
# The first column should have the name of the calibration (e.g., Afrotheria, 
# Laurasiatheria, etc.) as it will help you identify which plot belongs to which
# calibration. The second column is the node used in MCMCtree. The third column
# is the calibration used for that node in MCMCtree format.
# 
# [[ NOTES ABOUT ALLOWED CALIBRATION FORMATS]]
#
# Soft-bound calibrations: 
#  E.g.1: A calibration with a minimum of 0.6 and a maximum of 0.8 would with  
#         the default tail probabilities would have the following equivalent 
#         formats:
##        >> B(0.6,0.8) | B(0.6,0.8,0.025,0.025)
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
##        pU = 0.025:
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
# your input files (semicolon-separated files). Please save all your calibration 
# files (if more than one) in the same directory. The path to this directory is 
# what the argument `main_dir` needs. The argument `f_names` requires the name 
# of the file/s that you have used. Argument `dat` requires the same global 
# object that you have created at the beginning of the script. If your input  
# files have a header, please keep `head_avail = TRUE`. Otherwise, change this  
# to FALSE.
dat_ff <- list.files( path = "../calib_files", pattern = "csv",
                      full.names = FALSE )
calib_nodes <- read_calib_f( main_dir = paste( home_dir, "calib_files/",
                                               sep = "" ),
                             f_names = dat_ff,
                             dat = dat, head_avail = TRUE )

#-----------#
# LOAD DATA #
#-----------#
# Load mcmc files for all datasets
mcmc_obj <- vector( "list", num_dirs )
names( mcmc_obj ) <- rep( c( "CLK", "GBM", "ILN" ), length(dat) )
prior <- rep( c( TRUE, FALSE, FALSE ), length(dat) )
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
# Plot calibrated nodes without cross-bracing
for( k in 1:length(calib_nodes) ){
  name_labs   <- calib_nodes[[ k ]][,2]
  labs_in_csv <- paste( "t_n", calib_nodes[[ k ]][,2], sep = "" )
  # Start plotting!
  for( i in 1:length(labs_in_csv) ){
    pdf( file = paste( home_dir, "plots/priorVSpost_", dat[k], "_", 
                       calib_nodes[[ k ]][i,1], "_", name_labs[i], "_plot.pdf",
                       sep = "" ), 
         paper = "a4r", width = 0, height = 0 )
    # Get calibration density in correct formatting (4 nums)
    is_B    <- grep( x = calib_nodes[[ k ]][i,3], pattern = "B\\(" )
    is_U    <- grep( x = calib_nodes[[ k ]][i,3], pattern = "U\\(" )
    is_L    <- grep( x = calib_nodes[[ k ]][i,3], pattern = "L\\(" )
    is_ST   <- grep( x = calib_nodes[[ k ]][i,3], pattern = "ST\\(" )
    is_flag <- grep( x = calib_nodes[[ k ]][i,3], pattern = "flag" )
    if( length( is_B ) == 1 && length( is_flag ) == 0 ){
      # 4 numbers
      tmp_allcal <- as.numeric( stringr::str_split( string = gsub( pattern = "B\\(|\\)", 
                                                                   replacement = "",
                                                                   x = calib_nodes[[ k ]][i,3] ),
                                                    pattern = "," )[[1]] )
    }
    if( length( is_U ) == 1 && length( is_flag ) == 0 ){
      # 2 numbers
      tmp_allcal <- as.numeric( stringr::str_split( string = gsub( pattern = "U\\(|\\)", 
                                                                   replacement = "",
                                                                   x = calib_nodes[[ k ]][i,3] ),
                                                    pattern = "," )[[1]] )
    }
    if( length( is_L ) == 1 && length( is_flag ) == 0 ){
      # 2 numbers
      tmp_allcal <- as.numeric( stringr::str_split( string = gsub( pattern = "L\\(|\\)", 
                                                                   replacement = "",
                                                                   x = calib_nodes[[ k ]][i,3] ),
                                                    pattern = "," )[[1]] )
    }
    if( length( is_ST ) == 1 && length( is_flag ) == 0 ){
      # 4 numbers
      tmp_allcal <- as.numeric( stringr::str_split( string = gsub( pattern = "ST\\(|\\)", 
                                                                   replacement = "",
                                                                   x = calib_nodes[[ k ]][i,3] ),
                                                    pattern = "," )[[1]] )
    }
    if( length( is_flag ) == 1 ){
      tmp_allcal <- ""
    }
    # Get max/min x and y
    ##> NOTE: Only the first number is used, because the other number/s in 
    ##> `dup_dat[[i]` have the same fossil info!
    tn_matrix <- which( colnames(mcmc_obj$GBM$divt) %in% labs_in_csv[i] )
    max_y     <- max( c( max( density( mcmc_obj$GBM$divt[[ tn_matrix[1] ]], adj = 1 )$y ),
                         max( density( mcmc_obj$ILN$divt[[ tn_matrix[1] ]], adj = 1 )$y ),
                         max( density( mcmc_obj$CLK$divt[[ tn_matrix[1] ]], adj = 1 )$y ) ) )
    if( length( is_B ) == 1 && length( is_flag ) == 0 ){
      max_x <- max( c( max( density( mcmc_obj$GBM$divt[[ tn_matrix[1] ]], adj = 1 )$x ),
                       max( density( mcmc_obj$ILN$divt[[ tn_matrix[1] ]], adj = 1 )$x ),
                       max( density( mcmc_obj$CLK$divt[[ tn_matrix[1] ]], adj = 1 )$x ),
                       tmp_allcal[2]) )
      min_x <- min( c( min( density( mcmc_obj$GBM$divt[[ tn_matrix[1] ]], adj = 1 )$x ),
                       min( density( mcmc_obj$ILN$divt[[ tn_matrix[1] ]], adj = 1 )$x ),
                       min( density( mcmc_obj$CLK$divt[[ tn_matrix[1] ]], adj = 1 )$x ),
                       tmp_allcal[1] ) )
    }
    if( length( is_U ) == 1 && length( is_flag ) == 0 ){
      max_x <- max( c( max( density( mcmc_obj$GBM$divt[[ tn_matrix[1] ]], adj = 1 )$x ),
                       max( density( mcmc_obj$ILN$divt[[ tn_matrix[1] ]], adj = 1 )$x ),
                       max( density( mcmc_obj$CLK$divt[[ tn_matrix[1] ]], adj = 1 )$x ),
                       c( tmp_allcal[1] + 1 ) ) )
      min_x <- min( c( min( density( mcmc_obj$GBM$divt[[ tn_matrix[1] ]], adj = 1 )$x ),
                       min( density( mcmc_obj$ILN$divt[[ tn_matrix[1] ]], adj = 1 )$x ),
                       min( density( mcmc_obj$CLK$divt[[ tn_matrix[1] ]], adj = 1 )$x ) ) )
    }
    if( length( is_L ) == 1 && length( is_flag ) == 0 ){
      max_x <- max( c( max( density( mcmc_obj$GBM$divt[[ tn_matrix[1] ]], adj = 1 )$x ),
                       max( density( mcmc_obj$ILN$divt[[ tn_matrix[1] ]], adj = 1 )$x ),
                       max( density( mcmc_obj$CLK$divt[[ tn_matrix[1] ]], adj = 1 )$x ),
                       c( tmp_allcal[1]) + 2 ) )
      min_x <- min( c( min( density( mcmc_obj$GBM$divt[[ tn_matrix[1] ]], adj = 1 )$x ),
                       min( density( mcmc_obj$ILN$divt[[ tn_matrix[1] ]], adj = 1 )$x ),
                       min( density( mcmc_obj$CLK$divt[[ tn_matrix[1] ]], adj = 1 )$x ),
                       tmp_allcal[1] ) )
    }
    if( length( is_ST ) == 1 ){
      max_x <- max( c( max( density( mcmc_obj$GBM$divt[[ tn_matrix[1] ]], adj = 1 )$x ),
                       max( density( mcmc_obj$ILN$divt[[ tn_matrix[1] ]], adj = 1 )$x ),
                       max( density( mcmc_obj$CLK$divt[[ tn_matrix[1] ]], adj = 1 )$x ) ) )
      min_x <- min( c( min( density( mcmc_obj$GBM$divt[[ tn_matrix[1] ]], adj = 1 )$x ),
                       min( density( mcmc_obj$ILN$divt[[ tn_matrix[1] ]], adj = 1 )$x ),
                       min( density( mcmc_obj$CLK$divt[[ tn_matrix[1] ]], adj = 1 )$x ) ) )
    }
    if( length( is_flag ) == 1 ){
      max_x <- max( c( max( density( mcmc_obj$GBM$divt[[ tn_matrix[1] ]], adj = 1 )$x ),
                       max( density( mcmc_obj$ILN$divt[[ tn_matrix[1] ]], adj = 1 )$x ),
                       max( density( mcmc_obj$CLK$divt[[ tn_matrix[1] ]], adj = 1 )$x ) ) )
      min_x <- min( c( min( density( mcmc_obj$GBM$divt[[ tn_matrix[1] ]], adj = 1 )$x ),
                       min( density( mcmc_obj$ILN$divt[[ tn_matrix[1] ]], adj = 1 )$x ),
                       min( density( mcmc_obj$CLK$divt[[ tn_matrix[1] ]], adj = 1 )$x ) ) )
    }
    # Get name of calibration, which corresponds to the first calibration
    # in the list
    tmp_name_calib <- calib_nodes[[ k ]][i,1]
    # Start plotting calibration density, marginal density, and posterior densities for
    # each cross-braced node -- easy for comparisons!
    # Get temporary matching column in `mcmc_obj` to `j`
    plot( density( unlist( mcmc_obj$GBM$divt[ tn_matrix ] ), adj = 1 ),
          main = paste( "Post VS Marginal density - ",
                        labs_in_csv[i], " | ", tmp_name_calib, sep = "" ),
          col = "blue", ylim = c(0,max_y), xlim = c(min_x, max_x) )
    lines( density( unlist( mcmc_obj$CLK$divt[ tn_matrix ] ), adj = 1 ), col = "brown" )
    lines( density( unlist( mcmc_obj$ILN$divt[ tn_matrix ] ), adj = 1 ),
           col = "darkolivegreen3" )
    if( length( is_B ) == 1 && length( is_flag ) == 0 ){
      curve( mcmc3r::dB( x, tL = tmp_allcal[1], tU = tmp_allcal[2],
                         pL = tmp_allcal[3], pU = tmp_allcal[4] ),
             from = tmp_allcal[1]-0.5, to = tmp_allcal[2]+0.5,
             n = 1e5, add = TRUE, col = "black" )
    }
    if( length( is_U ) == 1 && length( is_flag ) == 0 ){
      curve( mcmc3r::dU( x, tU = tmp_allcal[1], pU = tmp_allcal[2] ),
             from = min_x, to = tmp_allcal[1]+1,
             n = 1e5, add = TRUE, col = "black" )
    }
    if( length( is_L ) == 1 && length( is_flag ) == 0 ){
      curve( mcmc3r::dL( x, tL = tmp_allcal[1], pL = tmp_allcal[2] ),
             from = tmp_allcal[1]-1, to = max_x,
             n = 1e5, add = TRUE, col = "black" )
    }
    if( length( is_ST ) == 1 && length( is_flag ) == 0 ){
      curve( sn::dst( x, xi = tmp_allcal[1], omega = tmp_allcal[2],
                      alpha = tmp_allcal[3], nu = tmp_allcal[4]),
             from = tmp_allcal[1]-0.5, to = tmp_allcal[1]+0.5,
             n = 1e5, add = TRUE, col = "black" )
    }
    legend( "topright", legend = c( paste( "Calibration density\n",
                                          calib_nodes[[ k ]][i,3], "\n",
                                          sep = "" ),
                                   "Marginal density",
                                   "Post-GBM", "Post-ILN"),
            lwd = 1, bty = "n", cex = 1,
            col = c( "black", "brown", "blue", "darkolivegreen3" ) )
    cols_vec <- c( "blue", "red", "purple", "darkgreen", "brown" )
    # Close pdf
    dev.off()
  }
}


