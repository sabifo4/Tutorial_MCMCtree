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
script_wd <- paste( dirname( path_to_file ), "/", sep = "" )
wd <- gsub( pattern = "/scripts", replacement = "", x = script_wd )
# Set working directory
setwd( wd )

#--------------#
# LOAD OBJECTS #
#--------------#
# The tips in the uncalibrated tree will be used to generate the dummy
# alignments
tt <- ape::read.tree( file = "../00_inp_data/tree_display/mtcdnapri_uncalib.tree" )

#-------#
# TASKS #
#-------#
# 1. Find number of taxa 
num_sp        <- length( tt$tip.label )
spnames       <- tt$tip.label
phylip_header <- paste( num_sp, "  1", sep = "" )

phylip_header_aln <- paste( num_sp, "  2\n", sep = "" )
spnames_2chars      <- paste( spnames, "     AT", sep = "" )

# 2. Generate dummy aln
if( ! dir.exists( "dummy_aln/" ) ){
  dir.create( "dummy_aln/" )
}
num_parts <- 3
for( i in 1:num_parts ){
  if( i == 1 ){
    write( x = phylip_header_aln, file = paste( "dummy_aln/dummy_aln.aln",
                                                sep = "" ) )
    write( x = spnames_2chars, file = paste( "dummy_aln/dummy_aln.aln",
                                             sep = "" ),
           append = TRUE )
  }else{
    write( x = paste( "\n", phylip_header_aln, sep = "" ),
           file = paste( "dummy_aln/dummy_aln.aln", sep = "" ),
           append = TRUE )
    write( x = spnames_2chars, file = paste( "dummy_aln/dummy_aln.aln",
                                             sep = "" ),
           append = TRUE )
  }
  
}


