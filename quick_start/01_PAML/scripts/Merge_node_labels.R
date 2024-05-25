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
# Read trees
## tt_list[[ j ]][[ 1 ]]: tree with node labels as of `MCMCtree`
## tt_list[[ j ]][[ 2 ]]: tree with node labels including `MCMCtree` notation
#                         and node name
# Put them in a list
dat <- c( "mtcdnapri" )
tt_list <- vector( mode = "list", length = 1 )
names( tt_list ) <- dat
for( i in 1:length( dat ) ){
  cat( "[[ Evaluating dataset ", dat[i], " ]]\n" )
  tt_list[[ i ]] <- vector( mode = "list", length = 2 )
  if( length( dat ) > 1 ){
    tt_list[[ i ]][[ 1 ]] <- ape::read.tree( file = paste( "../01_MCMCtree/00_prior/",
                                                           "node_tree_", i, ".tree",
                                                           sep = "" ) )  
  }else{
    tt_list[[ i ]][[ 1 ]] <- ape::read.tree( file = paste( "../01_MCMCtree/00_prior/",
                                                           "node_tree.tree",
                                                           sep = "" ) )
  }
  # Get the tree in a format that can be read by FigTree, where next to 
  # `MCMCtree` notation you will see the name you gave to such calibration
  tt_list[[ i ]][[ 2 ]] <- ape::read.tree( file = paste( "../../00_inp_data/tree_display/",
                                                         dat[i],
                                                          "_fordisplay_calib_MCMCtree.tree",
                                                          sep = "" ) )
 
}
# Start calibration information file for MCMC diagnostics with header
if( ! dir.exists( "../../00_inp_data/calibs/inp_calibs" ) ){
  dir.create( "../../00_inp_data/calibs/inp_calibs" ) 
}
for( i in 1:length(dat) ){
  write( x = "Calib;node;Prior",
         file = paste( "../../00_inp_data/calibs/inp_calibs/Calibnodes_",
                       dat[i], ".csv", sep = "" ) )
  write( x = "Calib;node;Prior",
         file = paste( "../../00_inp_data/calibs/inp_calibs/Calibnodes_",
                       dat[i], "_margVScalib.csv",
                       sep = "" ) )
}
# Now, populate the output file with the info for the rest of the nodes
for( j in 1:length(dat) ){
  # Start counter to see if that file has been created or not
  count_margVScalib <- 0
  cat( "\n\n[[ Evaluating tree for dataset ", dat[j], " ]]\n" )
  for( i in 1:length( tt_list[[ j ]][[ 2 ]]$node.label ) ){
    
    if( tt_list[[ j ]][[ 2 ]]$node.label[i] != "" ){
      tmp_name <- gsub( x = gsub( x = tt_list[[ j ]][[ 2 ]]$node.label[i],
                                  pattern = "..*\\)-",
                                  replacement = "" ), pattern = "'",
                        replacement =  "" )
      tmp_dist <- gsub( x = gsub( x = tt_list[[ j ]][[ 2 ]]$node.label[i],
                                  pattern = "\\)-..*",
                                  replacement = "\\)" ), pattern = "'",
                        replacement = "" )
      is_notdist <- grep( x = tmp_dist, pattern = "flag" )
      if( length( is_notdist ) > 0 ){
        write( x = paste( tmp_name, ";", tt_list[[ j ]][[ 1 ]]$node.label[i],
                          ";", tmp_dist, sep = "" ), 
               file = paste( "../../00_inp_data/calibs/inp_calibs/Calibnodes_",
                             dat[j], "_margVScalib.csv", sep = "" ),
               append = TRUE )
        count_margVScalib <- count_margVScalib + 1
      }
      cat( paste( tmp_name, ";", tt_list[[ j ]][[ 1 ]]$node.label[i], ";",
                  tmp_dist, "\n", sep = "" ) )
      write( x = paste( tmp_name, ";", tt_list[[ j ]][[ 1 ]]$node.label[i],
                        ";", tmp_dist, sep = "" ), 
             file = paste( "../../00_inp_data/calibs/inp_calibs/Calibnodes_",
                           dat[j], ".csv", sep = "" ), append = TRUE )
    }
    
  }
  # Remove second csv if not needed
  if( count_margVScalib == 0 ){
    unlink( x = paste( "../../00_inp_data/calibs/inp_calibs/Calibnodes_",
                       dat[j], "_margVScalib.csv", sep = "" ) )
  }
  
}


