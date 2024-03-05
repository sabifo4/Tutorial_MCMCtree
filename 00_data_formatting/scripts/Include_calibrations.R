#-------------------#
# CLEAN ENVIRONMENT #
#-------------------#
rm( list = ls( ) )

#-----------------------#
# SET WORKING DIRECTORY #
#-----------------------#
# Load package to help find the path to this source file 
library(rstudioapi) 
# Get the path of current open file
path_to_file <- getActiveDocumentContext()$path 
# Get working directory path
wd      <- paste( dirname( path_to_file ), "/", sep = "" )
wd.name <- dirname( path_to_file )
# Set wd. Note that the wd will be the same directory where this 
# script is saved.
setwd( wd )
# Load in-house R script with main function to add node age constraints
source( "../../src/Functions.R" )

#----------------------------#
# DEFINE GLOBAL VARS BY USER #
#----------------------------#
# Path to your input tree with calibrations. If in the same directory,
# you only need to write the name. If in a different directory, please
# type the absolute or relative path to this file. You need to include the
# flags within square brackets (e.g., [Mammalia]) and write them on the node
# that is to be calibrated.
#
# NOTE: Make always sure that there is at least one blank line at the 
# end of the this text file! Otherwise, you will get an error telling you that 
# there is an incomplete final line in these files.


# Path to your input text file that allows you to match the flags you have 
# used to label the nodes that are to be calibrated with the calibration you 
# want to use in `MCMCtree` format. If in the same directory,
# you only need to write the name. If in a different directory, please
# type the absolute or relative path to this file.
# The format of this file meets the following criteria:
#
#   - Header.
#   - One row per calibration.
#   - No spaces at all, semi-colon separated,
#   - There are 8 columns:
#       - Name of calibration (no spaces!).
#       - Name of tip 1 that leads to MRCA (no spaces!).
#       - Name of tip 2 that leads to MRCA (no spaces!).
#       - Minimum age.
#       - Percentage for the tail that will be allowed from minimum age.
#       - Maximum age.
#       - Percentage for the tail that will be allowed from maximum age.      
#       - If you want to specify directly the calibration in `MCMCtree` format,
#         then type this calibration in the 8th column and leave columns 5-7
#         blank (no spaces!).
# 
# E.g. 1: row in this text file to calibrate node "LBCA"
#
# ```
# # With columns 4-7
# LBCA;Dorm3;ProtG29;3.225;1e-4;4.520;1e-5;
# # With column 8 only
# LBCA;Dorm3;ProtG29;;;;;'B(3.225,4.520,1e-4,1e-5)'
# ```
#
# If in doubt, please follow the same format as used in the calibrations file 
# used for this analysis
path_textconv <- c( "../00_raw_data/calibs/calibrations.csv" )
all_calibs    <- vector( mode = "list", length(path_textconv) )
names( all_calibs ) <- c( "example" )
for( i in 1:length( all_calibs ) ){
  all_calibs[[ i ]] <- read.table( file = path_textconv[i],
                                   stringsAsFactors = FALSE, sep = ";",
                                   blank.lines.skip = TRUE, header = TRUE,
                                   colClasses = rep( "character", 8 ) )
}
# Path to trees
path_trees <- c( "../01_inp_data/tree_example_uncalib.tree" )
# Run `add_const` in-house pipeline to calibrate the tree
# 
# Arguments:
# tt            Phylo, object with the tree, previously generated
# calibrations  Matrix, element of the list generated with the calibration
#               info. E.g., `all_calibs[[1]]`
# out_name      Character, name of the dataset analysed
# out_dir_raw   Character, abs/rel path to the output directory for raw
#               trees
# out_dir_inp   Character, abs/rel path to th eoutput directory for inp data
#
for( i in 1:length( path_trees ) ){
  tt_ape <- ape::read.tree( file = path_trees[i] )
  test   <- add_node_const( tt = tt_ape, calibrations = all_calibs[[ i ]],
                            out_name = names( all_calibs )[i],
                            out_dir_raw = "../00_raw_data/trees/",
                            out_dir_inp = "../01_inp_data/")
  
}

