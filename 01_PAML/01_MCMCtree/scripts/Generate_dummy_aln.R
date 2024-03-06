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
# Get tree topology so that we can get the taxa names
tt <- ape::read.tree( file = "../../00_data_formatting/01_inp_data/tree_example_uncalib.tree" )

#-------#
# TASKS #
#-------#
# 1. Find number of taxa 
num_sp        <- length( tt$tip.label )
spnames       <- tt$tip.label
phylip_header <- paste( num_sp, "  1", sep = "" )

phylip_header_aln <- paste( num_sp, "  2\n", sep = "" )
spnames_2chars      <- paste( spnames, "     AT", sep = "" ) # Change this if you have AA!

# 2. Generate dummy aln
if( ! dir.exists( "dummy_aln/" ) ){
  dir.create( "dummy_aln/" )
}
write( x = phylip_header_aln, file = paste( "dummy_aln/dummy_aln.aln", sep = "" ) )
write( x = spnames_2chars, file = paste( "dummy_aln/dummy_aln.aln", sep = "" ),
       append = TRUE )

