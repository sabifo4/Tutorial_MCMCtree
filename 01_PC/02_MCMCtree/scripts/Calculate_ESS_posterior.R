#-------------------#
# CLEAN ENVIRONMENT #
#-------------------#
rm( list = ls( ) )

#-----------------------------------------------#
# LOAD PACKAGES, FUNCTIONS, AND SET ENVIRONMENT #
#-----------------------------------------------#
# Load package needed to automatically find the path to the "scripts"
# directory
library( rstudioapi )
scripts_dir   <- gsub( pattern = "scripts..*", replacement = "scripts/",
                       x = getActiveDocumentContext()$path )
# Load the file with all the functions used throughout this script
source( file = paste( scripts_dir, "Functions.R", sep = "" ) )
# Set home directory and output directory for ESS and convergence tests
home_dir      <- set_homedir()$home
outchecks_dir <- set_homedir()$ESS
# By now, set the working dirctory to `home_dir`
setwd( home_dir )

#--------------------------------#
# DEFINE USER'S GLOBAL VARIABLES #
#--------------------------------#
# First, we will define global variables that we will keep using throughout this
# tutorial. The values used here correspond to the example files for the tutorial
# so, if you are running this script with your own data, please modify them
# accordingly.

# Number of chains
num_chains <- 5 

# Number of divergence times. One trick to find this out quickly is to open the
# `mcmc.txt` file and check the header. The first element after `Gen` will have
# the format of `t_nX`, where X will be an integer (i.e., 5). Subtract two to 
# this number and this will be your number of divergence times that are parameters
# of the MCMC. Please modify the number below so it fits to the dataset you are using.
num_divt <- 3

# If you have a file called "mcmc_clean" file (i.e., the last line has all the
# entries because the `MCMCtree` job finished), then set the variable to TRUE.
# Otherwise, please type `FALSE`. If the latter, if there are any incomplete 
# lines, they will be ignored.
# Please modify the number below so it fits to the dataset you are using.
clean <- FALSE

# Quantile percentage that you want to set By default, the variable below is 
# set to 0.975 so the 97.5% and 2.5% quantiles. If you want to change this,
# however, just modify the value.
perc <- 0.975

# Number of columns in the `mcmc.txt` that are to be deleted as they do not 
# correspond to sample values for divergence times (i.e., the entries are not 
# names following the format `t_nX`). To figure out this number quickly, you 
# can open the `mcmc.txt` file, read the header, and count the number of `mu*`
# and `sigma2*` elements. Do not count the `lnL` value when looking at `mcmc.txt`
# files generated when sampling from the posterior -- this is automatically
# accounted for in the in-house R functions that you will subsequently use. 
# E.g., sssuming an MCMC ran under a relaxed-clock model with no partitions, 
# we would see `mu` and `sigma2` columns. Therefore, the variable would be set 
# to `delcol = 2`. Please modify the values below according to your dataset for 
# the `mcmc.txt` file generated when sampling from the prior (`delcol_prior`) 
# and when sampling from the posterior (`delcol_posterior`).
delcol_prior <- 1 # There is one column: `mu`
delcol_post  <- 2 # There are two column: `mu`, `sigma2`

# Path to the directory of each data alignment. In this case, we will have the 
# paths to one dataset analysed under the prior and the posterior. If that 
# changed, you will need to change the path accordingly depending on your dataset
# You may have as many paths as datasets you have analysed.
path_prior    <- paste( home_dir, "prior/1/CLK/", sep = "" )
path_post_GBM <- paste( home_dir, "posterior/1/GBM/", sep = "" )
path_post_ILN <- paste( home_dir, "posterior/1/ILN/", sep = "" )

#-----------#
# LOAD DATA #
#-----------#

#\\\\\\\#
# PRIOR #
#-------#
# 1. Obtain object with all chains and generate first convergence plot
prior_CLK_sum  <- find_prob_MCMC( num_dirs = num_chains, delcol = delcol_prior, 
                                  data_dir = path_prior,
                                  num_divt = num_divt, dataset = "prior_out",
                                  node_calib = paste( outchecks_dir,
                                                      "Only_calibnodes.csv",
                                                      sep = "" ),
                                  perc = perc, clean = clean,
                                  out_dat = path_prior, prior = TRUE )

half_chains <- num_chains/2
prior_CLK_sum_half1 <- apply( X = prior_CLK_sum$mean[1:half_chains,], MARGIN = 2,
                              FUN = mean )
prior_CLK_sum_half2 <- apply( X = prior_CLK_sum$mean[(half_chains+1):num_chains,],
                              MARGIN = 2, FUN = mean )
pdf( paste( outchecks_dir, "plots/Convergence_plot_priorCLK.pdf", sep = "" ),
     paper = "a4" )
plot_convergence( name_dir_subt = "Prior - CLK", mean_divt1 = prior_CLK_sum_half1,
                  mean_divt2 = prior_CLK_sum_half2, num_runs = num_chains )
dev.off()

#> CHECK: Looking for q97.5% and q2.5% issues across the chains. Basically, this 
#> function compares the divergence times estimated for each node across all 
#> chains. The difference between time estimates is then computed. Differences
#> larger than the threshold set by the user will be flagged, and hence the 
#> corresponding chain will be flagged as problematic.
#> 
#> If you are working with deep phylogenies, you may want to be more generous 
#> with the threshold (e.g., >0.4). Otherwise, a threshold lower than 0.4 can 
#> be too stringent and many chains will be flagged as "problematic". Bear in mind
#> that should not bee too vague (e.g., you would accept a different of +-threshold
#> in time estimates at the 97.5% and 2.5% quantile) as the larger the threshold 
#> the larger the differences between the estimates collected across chains for 
#> the same nodes.
#> For shallow datasets, you may try a threshold `th = 0.2` so you can better
#> refine which chains you keep. Once you run the ESS checks, you can make sure
#> that Rhat is not larger than 1.05 too.
#> If your analyses has returned flagged chains, please check the output file
#> `check_chains.txt` that will be generated before deciding whether a chain is
#> to be kept or deleted from your analaysis.
#> 
#> 1. Set threshold and run preliminary comparison anlaysis across MCMC runs
th <- 0.1
sum_quantiles <- check_quantiles( dat = prior_CLK_sum, threshold = th,
                                  num_chains = num_chains, compare_all = TRUE )
#> 2. Find out which is the one that should be used as "main chain" against which
#> the rest should be compared (i.e., the one with less differences in q97.5% and q2.5%)
chain_ind <- which( sum_quantiles$sum_chains[,1] %in% min( sum_quantiles$sum_chains[,1] ) )
main_ch   <- sum_quantiles$sum_chains[chain_ind,2]
if( length( main_ch ) > 1 ){
  # Use by default one of the chains that ran first
  # i.e., starting from chain labelled "1"
  main_ch <- min( main_ch )
}
#> 3. Find out which chains should be removed by using this main chain
#>    In this case, the main chain is chain-32
check_quantiles( dat = prior_CLK_sum, threshold = th, main_chain = main_ch,
                 num_chains = num_chains, compare_all = FALSE, outdir = path_prior )
#>
#> Chains to be deleted: none.

# 2. Get filtered summary object and generate convergence plots with filtered dataset
##> COMMENTED AS NO CHAINS HAD TO BE DELETED.
# filt_chains <- c( )
# prior_CLK_FILT_sum  <- find_prob_MCMC( num_dirs = filt_chains, delcol = delcol_prior, 
#                                        data_dir = path_prior,
#                                        num_divt = num_divt, dataset = "prior_FILT_out",
#                                        node_calib = paste( outchecks_dir,
#                                                            "Only_calibnodes.csv",
#                                                            sep = "" ),
#                                        perc = perc, clean = clean,
#                                        out_dat = path_prior, prior = TRUE )
# num_filt_chains  <- length( filt_chains )
# half_filt_chains <- as.integer( num_filt_chains/2 )
# if( half_filt_chains == 1 ){
#   prior_CLK_FILT_sum_half1 <- prior_CLK_FILT_sum$mean[1:half_filt_chains,]
# }else{
#   prior_CLK_FILT_sum_half1 <- apply( X = prior_CLK_FILT_sum$mean[1:half_filt_chains,],
#                                      MARGIN = 2, FUN = mean )
# }
# if( c(half_filt_chains+1) == 1 ){
#   prior_CLK_FILT_sum_half2 <- prior_CLK_FILT_sum$mean[c(half_filt_chains+1):num_filt_chains,]
# }else{
#   prior_CLK_FILT_sum_half2 <- apply( X = prior_CLK_FILT_sum$mean[c(half_filt_chains+1):num_filt_chains,],
#                                      MARGIN = 2, FUN = mean )
# }
# 
# pdf( paste( outchecks_dir, "plots/Convergence_plot_priorCLK_filt.pdf", sep = "" ), paper = "a4" )
# plot_convergence( name_dir_subt = "Prior - CLK", mean_divt1 = prior_CLK_FILT_sum_half1,
#                   mean_divt2 = prior_CLK_FILT_sum_half2, num_runs = num_filt_chains )
# dev.off()


#\\\\\\\\\\\\\\\\\#
# POSTERIOR - GBM #
#-----------------#
# 1. Obtain object with all chains and generate first convergence plot
posterior_GBM_sum  <- find_prob_MCMC( num_dirs = num_chains, delcol = delcol_post, 
                                      data_dir = path_post_GBM,
                                      num_divt = num_divt, dataset = "post_GBM_out",
                                      node_calib = paste( outchecks_dir,
                                                          "Only_calibnodes.csv",
                                                          sep = "" ),
                                      perc = perc, clean = clean,
                                      out_dat = path_post_GBM, prior = FALSE )

half_chains <- num_chains/2
posterior_GBM_sum_half1 <- apply( X = posterior_GBM_sum$mean[1:half_chains,], MARGIN = 2,
                              FUN = mean )
posterior_GBM_sum_half2 <- apply( X = posterior_GBM_sum$mean[(half_chains+1):num_chains,],
                              MARGIN = 2, FUN = mean )
pdf( paste( outchecks_dir, "plots/Convergence_plot_postGBM.pdf", sep = "" ),
     paper = "a4" )
plot_convergence( name_dir_subt = "Posterior - GBM", mean_divt1 = posterior_GBM_sum_half1,
                  mean_divt2 = posterior_GBM_sum_half2, num_runs = num_chains )
dev.off()

#> CHECK: Looking for q97.5% and q2.5% issues across the chains. Basically, this 
#> function compares the divergence times estimated for each node across all 
#> chains. The difference between time estimates is then computed. Differences
#> larger than the threshold set by the user will be flagged, and hence the 
#> corresponding chain will be flagged as problematic.
#> 
#> If you are working with deep phylogenies, you may want to be more generous 
#> with the threshold (e.g., >0.4). Otherwise, a threshold lower than 0.4 can 
#> be too stringent and many chains will be flagged as "problematic". Bear in mind
#> that should not bee too vague (e.g., you would accept a different of +-threshold
#> in time estimates at the 97.5% and 2.5% quantile) as the larger the threshold 
#> the larger the differences between the estimates collected across chains for 
#> the same nodes.
#> For shallow datasets, you may try a threshold `th = 0.2` so you can better
#> refine which chains you keep. Once you run the ESS checks, you can make sure
#> that Rhat is not larger than 1.05 too.
#> If your analyses has returned flagged chains, please check the output file
#> `check_chains.txt` that will be generated before deciding whether a chain is
#> to be kept or deleted from your analaysis.
#> 
#> 1. Set threshold and run preliminary comparison anlaysis across MCMC runs
th <- 0.1
sum_quantiles <- check_quantiles( dat = posterior_GBM_sum, threshold = th,
                                  num_chains = num_chains, compare_all = TRUE )
#> 2. Find out which is the one that should be used as "main chain" against which
#> the rest should be compared (i.e., the one with less differences in q97.5% and q2.5%)
chain_ind <- which( sum_quantiles$sum_chains[,1] %in% min( sum_quantiles$sum_chains[,1] ) )
main_ch   <- sum_quantiles$sum_chains[chain_ind,2]
if( length( main_ch ) > 1 ){
  # Use by default one of the chains that ran first
  # i.e., starting from chain labelled "1"
  main_ch <- min( main_ch )
}
#> 3. Find out which chains should be removed by using this main chain
#>    In this case, the main chain is chain-32
check_quantiles( dat = posterior_GBM_sum, threshold = th, main_chain = main_ch,
                 num_chains = num_chains, compare_all = FALSE, outdir = path_post_GBM )
#>
#> Chains to be deleted: NONE
#> END CHECK

# 2. Get filtered summary object and generate convergence plots with filtered dataset
# ##> COMMENTED AS NO CHAINS HAD TO BE DELETED.
# filt_chains <- c( )
# post_GBM_FILT_sum  <- find_prob_MCMC( num_dirs = filt_chains, delcol = delcol_post, 
#                                        data_dir = path_post_GBM,
#                                        num_divt = num_divt, dataset = "post_GBM_FILT_out",
#                                        node_calib = paste( outchecks_dir,
#                                                            "Only_calibnodes.csv",
#                                                            sep = "" ),
#                                        perc = perc, clean = clean,
#                                        out_dat = path_post_GBM, prior = FALSE )
# num_filt_chains  <- length( filt_chains )
# half_filt_chains <- as.integer( num_filt_chains/2 )
# if( length(1:half_filt_chains) == 1 ){
#   post_GBM_FILT_sum_half1 <- post_GBM_FILT_sum$mean[1:half_filt_chains,]
# }else{
#   post_GBM_FILT_sum_half1 <- apply( X = post_GBM_FILT_sum$mean[1:half_filt_chains,],
#                                      MARGIN = 2, FUN = mean )
# }
# if( length(c(half_filt_chains+1):num_filt_chains) == 1 ){
#   post_GBM_FILT_sum_half2 <- post_GBM_FILT_sum$mean[c(half_filt_chains+1):num_filt_chains,]
# }else{
#   post_GBM_FILT_sum_half2 <- apply( X = post_GBM_FILT_sum$mean[c(half_filt_chains+1):num_filt_chains,],
#                                      MARGIN = 2, FUN = mean )
# }
# 
# pdf( paste( outchecks_dir, "plots/Convergence_plot_postGBM_filt.pdf", sep = "" ), paper = "a4" )
# plot_convergence( name_dir_subt = "Posterior (filtered) - GBM", mean_divt1 = post_GBM_FILT_sum_half1,
#                   mean_divt2 = post_GBM_FILT_sum_half2, num_runs = num_filt_chains )
# dev.off()


#\\\\\\\\\\\\\\\\\#
# POSTERIOR - ILN #
#-----------------#
# 1. Obtain object with all chains and generate first convergence plot
posterior_ILN_sum  <- find_prob_MCMC( num_dirs = num_chains, delcol = delcol_post, 
                                      data_dir = path_post_ILN,
                                      num_divt = num_divt, dataset = "post_ILN_out",
                                      node_calib = paste( outchecks_dir,
                                                          "Only_calibnodes.csv",
                                                          sep = "" ),
                                      perc = perc, clean = clean,
                                      out_dat = path_post_ILN, prior = FALSE )

half_chains <- num_chains/2
posterior_ILN_sum_half1 <- apply( X = posterior_ILN_sum$mean[1:half_chains,], MARGIN = 2,
                                  FUN = mean )
posterior_ILN_sum_half2 <- apply( X = posterior_ILN_sum$mean[(half_chains+1):num_chains,],
                                  MARGIN = 2, FUN = mean )
pdf( paste( outchecks_dir, "plots/Convergence_plot_postILN.pdf", sep = "" ),
     paper = "a4" )
plot_convergence( name_dir_subt = "Posterior - ILN", mean_divt1 = posterior_ILN_sum_half1,
                  mean_divt2 = posterior_ILN_sum_half2, num_runs = num_chains )
dev.off()

#> CHECK: Looking for q97.5% and q2.5% issues across the chains. Basically, this 
#> function compares the divergence times estimated for each node across all 
#> chains. The difference between time estimates is then computed. Differences
#> larger than the threshold set by the user will be flagged, and hence the 
#> corresponding chain will be flagged as problematic.
#> 
#> If you are working with deep phylogenies, you may want to be more generous 
#> with the threshold (e.g., >0.4). Otherwise, a threshold lower than 0.4 can 
#> be too stringent and many chains will be flagged as "problematic". Bear in mind
#> that should not bee too vague (e.g., you would accept a different of +-threshold
#> in time estimates at the 97.5% and 2.5% quantile) as the larger the threshold 
#> the larger the differences between the estimates collected across chains for 
#> the same nodes.
#> For shallow datasets, you may try a threshold `th = 0.2` so you can better
#> refine which chains you keep. Once you run the ESS checks, you can make sure
#> that Rhat is not larger than 1.05 too.
#> If your analyses has returned flagged chains, please check the output file
#> `check_chains.txt` that will be generated before deciding whether a chain is
#> to be kept or deleted from your analaysis.
#> 
#> 1. Set threshold and run preliminary comparison anlaysis across MCMC runs
th <- 0.1
sum_quantiles <- check_quantiles( dat = posterior_ILN_sum, threshold = th,
                                  num_chains = num_chains, compare_all = TRUE )
#> 2. Find out which is the one that should be used as "main chain" against which
#> the rest should be compared (i.e., the one with less differences in q97.5% and q2.5%)
chain_ind <- which( sum_quantiles$sum_chains[,1] %in% min( sum_quantiles$sum_chains[,1] ) )
main_ch   <- sum_quantiles$sum_chains[chain_ind,2]
if( length( main_ch ) > 1 ){
  # Use by default one of the chains that ran first
  # i.e., starting from chain labelled "1"
  main_ch <- min( main_ch )
}
#> 3. Find out which chains should be removed by using this main chain
#>    In this case, the main chain is chain-32
check_quantiles( dat = posterior_ILN_sum, threshold = th, main_chain = main_ch,
                 num_chains = num_chains, compare_all = FALSE, outdir = path_post_ILN )
#>
#> Chains to be deleted: NONE
#> END CHECK

# 2. Get filtered summary object and generate convergence plots with filtered dataset
# ##> COMMENTED AS NO CHAINS HAD TO BE DELETED.
# filt_chains <- c( )
# post_ILN_FILT_sum  <- find_prob_MCMC( num_dirs = filt_chains, delcol = delcol_post, 
#                                        data_dir = path_post_ILN,
#                                        num_divt = num_divt, dataset = "post_ILN_FILT_out",
#                                        node_calib = paste( outchecks_dir,
#                                                            "Only_calibnodes.csv",
#                                                            sep = "" ),
#                                        perc = perc, clean = clean,
#                                        out_dat = path_post_ILN, prior = FALSE )
# num_filt_chains  <- length( filt_chains )
# half_filt_chains <- as.integer( num_filt_chains/2 )
# if( length(1:half_filt_chains) == 1 ){
#   post_ILN_FILT_sum_half1 <- post_ILN_FILT_sum$mean[1:half_filt_chains,]
# }else{
#   post_ILN_FILT_sum_half1 <- apply( X = post_ILN_FILT_sum$mean[1:half_filt_chains,],
#                                      MARGIN = 2, FUN = mean )
# }
# if( length(c(half_filt_chains+1):num_filt_chains) == 1 ){
#   post_ILN_FILT_sum_half2 <- post_ILN_FILT_sum$mean[c(half_filt_chains+1):num_filt_chains,]
# }else{
#   post_ILN_FILT_sum_half2 <- apply( X = post_ILN_FILT_sum$mean[c(half_filt_chains+1):num_filt_chains,],
#                                      MARGIN = 2, FUN = mean )
# }
# 
# pdf( paste( outchecks_dir, "plots/Convergence_plot_postILN_filt.pdf", sep = "" ), paper = "a4" )
# plot_convergence( name_dir_subt = "Posterior (filtered) - ILN", mean_divt1 = post_ILN_FILT_sum_half1,
#                   mean_divt2 = post_ILN_FILT_sum_half2, num_runs = num_filt_chains )
# dev.off()

#---------------#
# CALCULATE ESS #
#---------------#
#> ESS with RStan
# Each column is assumed to be an MCMC. Rows are iterations for parameter X
# Source explaining why it is preferable than the function in coda:
# https://nature.berkeley.edu/~pdevalpine/MCMC_comparisons/nimble_MCMC_comparisons.html

#\\\\\\\\\\\\\#
# PRIOR - CLK #
#-------------#
ESS_prior_CLK  <- sum_MCMC_ESS( x = prior_CLK_sum$all_mcmc, coda_fun = TRUE )
ESS_prior_CLK$tab
#      Tail-ESS Bulk-ESS  coda-ESS
# Med.    49494    49366  99756.04
# Min.    48780    49310  98554.67
# Max.    49790    49374 100000.00
min( ESS_prior_CLK$stats$Rhat )
#[1] 1.000011
max( ESS_prior_CLK$stats$Rhat )
#[1] 1.00008
dim( prior_CLK_sum$all_mcmc )
#[1] 100000      3


#\\\\\\\\\\\\\\\\#
# POSTERIOR- GBM #
#----------------#
ESS_post_GBM  <- sum_MCMC_ESS( x = posterior_GBM_sum$all_mcmc, coda_fun = TRUE )
ESS_post_GBM$tab
#      Tail-ESS Bulk-ESS  coda-ESS
# Med.    48314    50191 100000.00
# Min.    48126    49703  98946.45
# Max.    50138    50594 100000.00
min( ESS_post_GBM$stats$Rhat )
#[1] 0.9999802
max( ESS_post_GBM$stats$Rhat )
#[1] 0.999999
dim( posterior_GBM_sum$all_mcmc )
#[1] 100000      3


#\\\\\\\\\\\\\\\\#
# POSTERIOR - ILN #
#----------------#
ESS_post_ILN  <- sum_MCMC_ESS( x = posterior_ILN_sum$all_mcmc, coda_fun = TRUE )
ESS_post_ILN$tab
#      Tail-ESS Bulk-ESS  coda-ESS
# Med.    47944    49261   100000
# Min.    47643    48946   100000
# Max.    49876    50390   100000
min( ESS_post_ILN$stats$Rhat )
#[1] 0.9999814
max( ESS_post_ILN$stats$Rhat )
#[1] 1.00004
dim( posterior_ILN_sum$all_mcmc )
#[1] 100000      3
