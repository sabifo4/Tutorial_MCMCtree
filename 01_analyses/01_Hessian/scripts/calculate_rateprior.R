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
raw_tt <- ape::read.tree( file = "../00_data/00_raw_data/tree_ML.tree" )

#-------#
# TASKS #
#-------#
# 1. Find tree height. You can use the function `phytools::nodeHeights` to calculate 
#    all the heights of the tree. Then, we can extract the maximum height calculated,
#    which will correspond to the length from the root to the highest tip.
tree_height <- max( phytools::nodeHeights( raw_tt ) ) #0.3663503

# 2. Get the mean of the calibration set for the root to have the 
#    time for the speciation event at the root, what we will use 
#    to estimate the mean evolutionary rate later. As we have an
#    ST distribution, we will get the first parameter as the mean
#    age. We will test two different time units: 1Myr and 100Myr.
root_age1   <- 583
root_age100 <- 5.83

# 3. Estimate mean rate based on the two different time units
#
#    tree_height = mean_rate x root_age --> mean_rate = tree_height / root_age
mean_rate1   <- tree_height / root_age1   # 0.0006283882 subst/site per time unit
mean_rate100 <- tree_height / root_age100 # 0.06283882 subst/site per time unit

# If we want to know the mean rate in subst/site/year, we apply the time unit. We
# We should get the same estimate regardless the time unit used:
#
# Time unit 1Myr (10^6y): 0.0006283882 subst/site/10^6 = 6.28e-10 subst/site/year
# Time unit 100Myr (10^8y): 0.06283882 subst/site/10^8 = 6.28e-10 subst/site/year

# 4. Now, we can build the gamma distribution given that we now have an estimate
#    for the mean rate. We will also use `alpha = 2` as we will start with a 
#    vague distribution. Nevertheless, if you were very sure about the mean 
#    rate, you could build a more constraint prior.
#
#    mean_Gamma = mean_rate = alpha / beta --> beta = alpha / mean_rate
alpha    <- 2
beta1    <- alpha/mean_rate1   # 3182.746 ~ 3182
beta100  <- alpha/mean_rate100 # 31.83 ~ 32

# Here, the time unit 100Myrs looks better than 1Myrs as the times and rates
# are neither too large nor too small. In theory, the results should not
# depend on the chosen time unit. Nevertheless, we will choose the results with
# time unit = 100Myr and hence use G(2,32) for a vague distribution.
curve( dgamma( x, shape = 2, rate = beta100 ), from = 0, to = 0.5, col = "black" )
legend( "topright", legend = c( "G(2,32) " ), 
        col = "black", lty = 1, box.lty = 2 )

# 5. If we were to use a more constraint distribution (i.e., we are quite sure 
#    about the estimated mean rate being this value), we can change the alpha
#    parameter to a higher value to do so.
alpha_const1 = 10
beta_const1 = alpha_const1 / mean_rate100
alpha_const2 = 50
beta_const2 = alpha_const2 / mean_rate100
alpha_const3 = 100
beta_const3 = alpha_const3 / mean_rate100
curve( dgamma( x, shape = alpha_const3, rate = beta_const3 ), from = 0, to = 0.2, col = "green" )
curve( dgamma( x, shape = alpha_const2, rate = beta_const2 ), col = "blue", add = TRUE )
curve( dgamma( x, shape = alpha_const1, rate = beta_const1 ), col = "purple", add = TRUE)
curve( dgamma( x, shape = alpha, rate = beta100 ), col = "black", add = TRUE )
legend( "topright", legend = c( "G(2,32)", "G(10,159)", "G(50,796)", "G(100,1591)" ), 
        col = c("black", "purple", "blue","green"),
        lty = 1, box.lty = 2 )

# The plot above shows how confident you are with the mean rate being around
# 0.06283882 subst/site per 100Myr. The narrower the distribution is, the more
# confident you are

# 6. Plot the gamma distributions
if ( ! dir.exists( "out_RData" ) ){
  dir.create( "out_RData" )
}
pdf( file = "out_RData/gamma_dists.pdf", paper = "a4" )
par( mfrow = c(1,2) )
curve( dgamma( x, alpha, rate = beta100 ), from = 0, to = 2, col = "black" )
legend( "topright", legend = c( "G(2,4) " ), 
        col = "black", lty = 1, box.lty = 2 )
curve( dgamma( x, shape = alpha_const3, rate = beta_const3 ), from = 0, to = 2, col = "green" )
curve( dgamma( x, shape = alpha_const2, rate = beta_const2 ), col = "blue", add = TRUE )
curve( dgamma( x, shape = alpha_const1, rate = beta_const1 ), col = "purple", add = TRUE)
curve( dgamma( x, shape = alpha, rate = beta100 ), col = "black", add = TRUE )
legend( "topright", legend = c( "G(2,4)", "G(10,22)", "G(50,108)", "G(100,216)" ), 
        col = c("black", "purple", "blue","green"),
        lty = 1, box.lty = 2 )
dev.off()

