#-------------------#
# CLEAN ENVIRONMENT #
#-------------------#
rm( list = ls( ) )

#----------------#
# LOAD LIBRARIES #
#----------------#
library( ape )

#-----------------------#
# SET WORKING DIRECTORY #
#-----------------------#
library( rstudioapi ) 
# Get the path to current open R script and find main dir "00_Gene_filtering"
path_to_file <- getActiveDocumentContext()$path
script_wd <- paste( dirname( path_to_file ), "/", sep = "" )
main_dir  <- gsub( pattern = "/scripts", replacement = "", x = script_wd )
wd        <- paste( main_dir, "00_raw_data/", sep = "" )
setwd( wd )

#----------#
# GET DATA #
#----------#
# Load gene trees
trees_dir   <- paste( wd, "trees/", sep = "" )
fname_trees <- list.files( path = trees_dir, pattern = "tree" )
all_trees   <- vector( mode = "list", length( fname_trees ) )

# We know that the calibration that we will use on the root of the phylogeny
# will be a soft bound: `B(0.94,1.2537698)`. We will convert this calibration so
# time unit is 100 Myr given that the we will obtain smaller values for the beta
# parameter of the gamma distribution for the prior on the rates, which will be
# more convenient when subsequently running `MCMCtree`. We will take the mean of
# the min and max age used to set the soft bound as the mean root age. This will
# be the estimated time of divergence at the root of the phylogeny, which we will
# use, together with the tree height of each gene tree, to estimate the mean
# evolutionary rate of each gene tree
mean_root_age <- mean( c(0.94, 1.2537698) ) # 1.096885 (100 Myr)

# Generate object to keep track of sum data
sum_out  <- as.data.frame( matrix( 0, nrow = length( all_trees ), ncol = 6 ) )

# The order of the genes corresponds to the order of the alignments, 
# so we can match their names
colnames( sum_out ) <- c( "Tree.length", "%Largest.rel.blength",
                          "Tip.name", "Delete?", "Gene_name", "Evol_rate" )

# Create dirs for out results 
if( ! dir.exists( paste( wd, "out_RData", sep = "" ) ) ){
  dir.create(  paste( wd, "out_RData", sep = "" ) )
}
if( ! dir.exists( paste( wd, "out_logs", sep = "" ) ) ){
  dir.create(  paste( wd, "out_logs", sep = "" ) )
}

#---------------------------------#
# RUN RELATIVE BRANCH LENGTH TEST #
#---------------------------------#
# Start search for outliers running the relative branch length test
counter <- 0
for( i in 1:length( fname_trees ) ){
  
  all_trees[[ i ]]       <- ape::read.tree( file = paste( trees_dir,
                                                          fname_trees[i],
                                                          sep = "" ) )
  sum_out[,5]            <- gsub( x = fname_trees, pattern = ".tree",
                                  replacement = "" ) 
  tmp_theight            <- max( phytools::nodeHeights( all_trees[[ i ]] ) )
  tmp_meanrate           <- tmp_theight / mean_root_age
  sum_out[i,6]           <- tmp_meanrate
  names( all_trees )[i]  <- paste( "gene_", i, sep = "" )
  rownames( sum_out )[i] <- paste( "gene_", i, sep = "" )
  sum_out[i,1]           <- sum( all_trees[[ i ]]$edge.length )

  # Divide all branches into the tree length so we can obtain the relative 
  # branch lengths, i.e., r_{ij} = b{ij} / SUM_j (b_{ij}); as detailed in
  # dos Reis et al. 2012.
  # Then, take the largest relative branch length (in %) for each gene and store 
  # it in second  column in data.frame `sum_out`.
  sum_out[i,2] <- max( all_trees[[i]]$edge.length/sum_out[i,1] )*100
  print( max( all_trees[[i]]$edge.length/sum_out[i,1] )*100 )
  
  # If the largest relative branch length is equal to or larger than 60%...
  if ( round( sum_out[i,2] ) >= 60 ){
    cat( "Gene ", rownames( sum_out)[i], 
         "- has a branch length longer than 60% of the total tree length\n" )
    write( paste( "Gene ", rownames( sum_out )[i],
                  "- has a branch length longer than 60% of the total tree length: ",
                  sum_out[i,2], "%\n", sep = ""),
           file = paste( wd, "out_logs/log_genes_blength_longer_60pc.txt", sep = "" ),
           append = T )
    counter <- counter + 1
    sum_out[i,5] <- c("Y")
    
    # 1. Find which position has largest branch length
    ind1   <- which( round( all_trees[[i]]$edge.length/sum_out[i,1] ) >= 0.6 )
    # 2. In case more than one branch length was larger than 0.6 and stored as
    #   `ind1`, find the largest and save it as `ind1_1`
    ind1_1 <- max( c(all_trees[[i]]$edge.length/sum_out[i,1])[ind1] )
    # 3. Find the position of this largest branch which index is `ind1_1`, just
    #    in case there were more than one saved in step 1
    ind1_2 <- which( all_trees[[i]]$edge.length/sum_out[i,1] == ind1_1 )
    # 4. Find the node position of the individual which branch length is `ind1_2`
    ind2   <- all_trees[[i]]$edge[ind1_2,2]
    sum_out[i,3] <- all_trees[[i]]$tip.label[ind2]
    
    # 5. Now we want to find the corresponding taxon name. But if
    #    ind2 is larger than the amount of taxa available in
    #    this gene tree, then increase ind1_2 by 1 until
    #    it finds the correct taxon name.
    if( ind2 >= length( all_trees[[i]]$tip.label ) ){
      while( ind2 >= length( all_trees[[i]]$tip.label ) ){
        ind1_2 <- ind1_2 + 1
        ind2   <- all_trees[[i]]$edge[ind1_2,2]
      }
      sum_out[i,3] <- paste( all_trees[[i]]$tip.label[ind2], ",",
                              all_trees[[i]]$tip.label[ind2+1], sep = "" )
    }else{
      sum_out[i,3] <- all_trees[[i]]$tip.label[ind2]
    }
  }else{
    sum_out[i,3] <- "NULL"
    sum_out[i,4] <- c("N")
  }
  
  # Append that the analysis has finished in the log file!
  if( i == length( all_trees ) ){
    cat("\nA total of ", counter, "genes had a branch length longer than 60% of
        the total tree length\n")
    write( paste( "\nA total of ", counter, " genes had a branch length longer than 60% of
                  the total tree length\n", sep = "" ),
           file = paste( wd, "out_logs/log_genes_blength_longer_60pc.txt", sep = "" ),
           append = T )
  }
  
}

# Save RData for later
save( sum_out, file = "out_RData/sum_out.RData" )
save( all_trees, file = "out_RData/all_trees.RData" )
# [Run only if you have already generated the `sum_out.RData`]
# The command below loads the RData file
#load( "out_RData/sum_out.RData" )
#load( "out_RData/all_trees.RData" )

# Plot results
pdf( file = "out_RData/check_relblVStreelength.pdf", paper = "a4" )
plot( sum_out[,1], sum_out[,2],
      main = "Quality control check of relative branch lengths",
      xlab = "Tree length", ylab = "Relative branch lengths (%)" )
dev.off()

pdf( file = "out_RData/check_relblVStreelength_log.pdf", paper = "a4" )
plot( log(sum_out[,1]), sum_out[,2],
      main = "Quality control check of relative branch lengths",
      xlab = "Total tree length (log-scale)",
      ylab = "Relative branch lengths (%)" )
dev.off()

# There seem to be two outliers, so we may want to remove them.
# So you can see how different they are to the rest, however, we will keep 
# these taxa for now and continue the analyses as if we had not yet
# found them

# Save csv file with results
write.table( sum_out[,c(1:2,5:6)],
             file = "out_RData/sum_out.csv",
             quote = F, sep = "," )

#-----------------------------------------#
# ORDER GENES FROM SLOW- TO FAST-EVOLVING #
#-----------------------------------------#
# Now, we can order the trees from slowest- to fastest-evolving
sum_out_s2f  <- sum_out[order( round( as.numeric( sum_out$Evol_rate ), 3) ),]
# Slowest and fastest gene in each data subset
round( as.numeric( min( sum_out_s2f$Evol_rate ) ), 3 )
round( as.numeric( max( sum_out_s2f$Evol_rate ) ), 3 )
# (0.125, 15.656)
# Now, you can see that `gene_10` and `gene_5` seem to be the two outliers you
# had already detected with the visual plot. There cannot be enough checks!
# Let's delete these two genes :)
sum_out_s2f <- sum_out_s2f[-c(14,15),]
# Plot filtered results
pdf( file = "out_RData/check_relblVStreelength_filt.pdf", paper = "a4" )
plot( sum_out_s2f[,1], sum_out_s2f[,2],
      main = "Quality control check of relative branch lengths",
      xlab = "Tree length", ylab = "Relative branch lengths (%)" )
dev.off()

pdf( file = "out_RData/check_relblVStreelength_log_filt.pdf", paper = "a4" )
plot( log(sum_out_s2f[,1]), sum_out_s2f[,2],
      main = "Quality control check of relative branch lengths",
      xlab = "Total tree length (log-scale)",
      ylab = "Relative branch lengths (%)" )
dev.off()
# Write tables with the matrices
write.table( sum_out_s2f[c(1,2,5)],
             file = "out_RData/filtered_s2f.csv",
             quote = F, sep = "," )

#---------------------------------------------------------#
# Copy gene alignments and gene trees in a file structure # 
# that ranks them from slow- to fast-evolving             #
#---------------------------------------------------------#
# Generate the data
# Check if directory where to copy files exist. Otherwise, create it.
if ( ! dir.exists( "../01_ordered_genes" ) ){
  dir.create( "../01_ordered_genes" )
}
# Copy genes and trees in the corresponding directories
name_aln <- gsub( x = list.files( path = paste( wd, "aln", sep = "" ),
                                  pattern = ".fasta" ),
                  pattern = ".fasta", replacement = "" )
counter <- 0
for( i in 1:length(rownames( sum_out_s2f )) ){
  # Start counter
  counter <- counter + 1
  # Check if directory where to copy files exist. Otherwise, create it.
  if ( ! dir.exists( paste( "../01_ordered_genes/", counter, "/", sep = "") ) ){
    dir.create( paste( "../01_ordered_genes/", counter, "/", sep = "" ) )
  }
  # Rename tags, remove duplicates, and copy resulting alignments
  tmp_g    <- sum_out_s2f[i,5]
  tmp_ind  <- which( name_aln %in% tmp_g ) 
  tmp_aln  <- paste( name_aln[tmp_ind], ".fasta", sep = "" )
  file.copy( paste( wd, "aln/", tmp_aln, sep = "" ),
             paste( "../01_ordered_genes/", counter, "/", sep = ""),
             recursive = TRUE )
  tmp_g_ind <- as.numeric( gsub( pattern = "gene_", replacement = "", 
                                 x = rownames( sum_out_s2f )[i] ) )
  cat( "Gene in position ", tmp_g_ind, "in the main tree file is ranked in position ", 
       counter, " with regards to the evolutionary rate. This gene is: ", tmp_g, "\n" )
  write( paste( "Gene in position ", tmp_g_ind, " in the main tree file is ranked in position ", 
                counter, " with regards to the evolutionary rate. This gene is: ", tmp_g, "\n",
                sep = "" ),
         file = "out_logs/log_R_copy_ordered_genes_s2f.txt",
         append = TRUE, sep = "\n" )
  # Copy gene tree
  ape::write.tree( all_trees[[tmp_g_ind]],
                   file = paste( "../01_ordered_genes/", counter, "/", tmp_g, ".tree", sep = "" ) )
}

  