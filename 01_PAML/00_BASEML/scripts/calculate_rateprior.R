#-------------------#
# CLEAN ENVIRONMENT #
#-------------------#
rm( list = ls( ) )

#-----------------------#
# SET WORKING DIRECTORY #
#-----------------------#
library( rstudioapi ) 
# Get the path to current open R script and find main dir "00_Gene_filtering"
path_to_file <- getActiveDocumentContext()$path
script_wd <- paste( dirname( path_to_file ), "/", sep = "" )
wd <- gsub( pattern = "/scripts", replacement = "", x = script_wd )
setwd( wd )

#--------------#
# LOAD OBJECTS #
#--------------#
# Load tree
raw_tt_paths <- list.files( path = "../../00_data_formatting/00_raw_data/trees/",
                            pattern = "tree_ML", full.names = TRUE )
raw_tt <- vector( mode = "list", length = length( raw_tt_paths ) )
for( i in 1:length(raw_tt) ){
  raw_tt[[ i ]] <- ape::read.tree( file = raw_tt_paths[i] )
}

#-------#
# TASKS #
#-------#
# 1. Get an estimate of the calibration set for the root to have the 
#    time for the speciation event at the root, what we will use 
#    to estimate the mean evolutionary rate later. If we had a soft-bound
#    calibration ir a hard-bound calibration, we would just get the mean of
#    the minimum and maximum ages and use an approximate mean root age.
#    We have a ST distribution, and so we can just use the first parameter
#    value (xi), which is an approximation of the "location" which, in this
#    case, would correspond to the mean root age
#    'ST(5.83,0.059,0.112,109.124)' --> xi = 5.83
root_age <- 5.83 # 5.83 (x100 Mya)

# 2. Find tree height. You can use the function `phytools::nodeHeights` to
#    calculateall the heights of the tree. Then, we can extract the maximum
#    height calculated, which will correspond to the length from the root to
#    the highest tip.
tree_height <- vector( mode = "numeric", length( raw_tt ) )
for( i in 1:length(raw_tt) ){
  tree_height[i] <- max( phytools::nodeHeights( raw_tt[[ i ]] ) )
}
# tree height: 0.3663503

# 3. Estimate mean rate based on the two different time units
#
#    tree_height = mean_rate x root_age --> mean_rate = tree_height / root_age
mean_rate <- tree_height / root_age # 0.06283882 subst/site/time unit

# If we want to know the mean rate in subst/site/year, we apply the time unit. We
# We should get the same estimate regardless the time unit used:
#
# Time unit 100 May (10^8y): 0.06283882 subst/site/1e+08 = 6.28e-10 subst/site/year

# 4. Now, we can build the gamma distribution given that we now have an estimate
#    for the mean rate. We will also use `alpha = 2` as we will start with a 
#    vague distribution. Nevertheless, if you were very sure about the mean 
#    rate, you could build a more constraint prior.
#
#    mean_Gamma = mean_rate = alpha / beta --> beta = alpha / mean_rate
alpha <- 2
beta  <- alpha/mean_rate # 31.82746 ~ 31.8

# We can plot this distribution
curve( dgamma( x, shape = 2, rate = beta ), from = 0, to = 3 )
legend( "topright", legend = c( "G(2,31.8) " ), 
        col = "black", lty = 1, box.lty = 2 )

# 5. Plot the gamma distributions
if ( ! dir.exists( "out_RData" ) ){
  dir.create( "out_RData" )
}
pdf( file = "out_RData/gamma_dist.pdf", paper = "a4" )
curve( dgamma( x, shape = 2, rate = beta), from = 0, to = 3, col = "black" )
legend( "topright", legend = c( "G(2,31.8) " ), 
        col = "black", lty = 1, box.lty = 2 )
dev.off()

