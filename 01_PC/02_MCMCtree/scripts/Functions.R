#--------------------#
# DEFINE GLOBAL FUNS #
#--------------------#

# Function to set the main working directories
# No arguments needed
#
set_homedir <- function( ){
  # Get the path to current open R script and find main dir "scripts"
  path_to_file <- getActiveDocumentContext()$path
  wd <- paste( dirname( path_to_file ), "/", sep = "" )
  wd <- gsub( pattern = "scripts/", replacement = "ESS_and_chains_convergence/", x = wd )
  if( ! dir.exists( wd ) ){
    dir.create( wd )
  }
  if( ! dir.exists( paste( wd, "plots", sep = "" ) ) ){
    dir.create( paste( wd, "plots", sep = "" )  )
  }
  wd2 <- gsub( pattern = "ESS_and_chains_convergence/", replacement = "", x = wd )
  return( list( home = wd2, ESS = wd ) )
}


# Function to convert the matrix with the sampled values
# for each parameter in the MCMC into a 3D array readable 
# by `rstan::monitor`.
# 
# Parameters:
# x  Data frame, there are many rows as iterations and as many 
#    columns as parameters that had to be inferred during the 
#    MCMC
#
mat2arr <- function( x )
{
  
  # Get rid of "Gen" column:
  ## Not needed anymore
  #x <- x[,-c(1)]
  
  # Get array format
  arr <- array( 0, dim = c( dim( x )[1], 1, dim( x )[2] ),
                dimnames = list( paste( "Iter", 1:dim( x )[1], sep = "" ),
                                 "1 chain",
                                 #paste( "Param", 1:dim( x )[2], sep = "" ) )
                                 colnames( x ) )
  )
  
  for( p in 1:dim( x )[2] ){
    arr[,1,p] <- x[,p]
  }
  
  # Return object
  return( arr )
  
}


# Function to run `rstan::monitor` and export the 
# median, min, and max values for the tail-ESS and the 
# bulk-ESS
#
# Parameters:
# x        Data frame, there are many rows as iterations and as many 
#          columns as parameters that had to be inferred during the 
#          MCMC
#
# coda_fun Boolean, TRUE if you want to also compute the ESS
#          with the coda::effectiveSize function
sum_MCMC_ESS <- function( x, coda_fun = FALSE )
{
  
  # Define empty vector to report stats
  if( coda_fun == FALSE ){
    out_stats <- matrix( 0, nrow = 3, ncol = 2 )
    colnames( out_stats ) <- c( "Tail-ESS", "Bulk-ESS" )
  } else{ 
    out_stats <- matrix( 0, nrow = 3, ncol = 3 )
    colnames( out_stats ) <- c( "Tail-ESS", "Bulk-ESS", "coda-ESS" )
  }
  rownames( out_stats ) <- c( "Med.", "Min.", "Max." )
  
  # Compute stats
  stats_mcmc     <- rstan::monitor( sims = mat2arr( x ) )
  out_stats[1,1] <- median( stats_mcmc$Tail_ESS )
  out_stats[2,1] <- min( stats_mcmc$Tail_ESS )
  out_stats[3,1] <- max( stats_mcmc$Tail_ESS )
  out_stats[1,2] <- median( stats_mcmc$Bulk_ESS )
  out_stats[2,2] <- min( stats_mcmc$Bulk_ESS )
  out_stats[3,2] <- max( stats_mcmc$Bulk_ESS )
  
  if( coda_fun == TRUE ){
    ESS_coda       <- coda::effectiveSize( x = x )
    out_stats[1,3] <- median( ESS_coda )
    out_stats[2,3] <- min( ESS_coda )
    out_stats[3,3] <- max( ESS_coda )
  }
  
  # Return stats 
  return( list( tab = out_stats, stats = stats_mcmc, stats_CODA = ESS_coda ) )
  
}


# Function to load the `mcmc.txt` files.
# Used to get mean estimates for convergence plots.
# 
# Parameters:
# mcmc    Character, path to mcmc.txt file
# delcol  Numeric, number of columns that are to be deleted 
#         as they do not contain divtimes. You can read the header
#         of the `mcmc.txt` files and count the number of `mu*` and
#         `sigma2*` elements. Do not count the `lnL` value as this is
#         automatically deleted in the calculations below. 
#         Assuming an MCMC run under a relaxed-clock model with no 
#         partitions, we would see `mu` and `sigma2` columns. Therefore,
#         `delcol = 2`. The default value is 2, but you will need to change
#         it if you run `MCMCtree` with different settings.
# perc    Numeric. Percentile to calculate the quantiles. Default: 0.975.
# clean   Boolean. FALSE if a "clean" file removing incomplete lines in the
#         `mcmc.txt` has not been created. TRUE otherwise.
# prior   Boolean. TRUE if `MCMCtree` has been sampling from the prior, FALSE
#         otherwise
load_dat <- function( mcmc, delcol = 3, perc = 0.975, clean = FALSE,
                      prior = FALSE )
{
  
  # 1. Load files and get parameters
  cat( "Load combined mcmc.txt from path \n", mcmc, "\n... ...\n" )
  if( clean == TRUE){
    run <- read.table( mcmc, header = TRUE, sep = "\t" )
  }else{
    # Use arg `fill` so empty fields in last line do not cause an error
    # when reading the file, then delete last line
    run_tmp <- read.table( mcmc, header = TRUE, fill = TRUE, sep = "\t" )
    run     <- matrix( 0, nrow = c(dim( run_tmp )[1]-1), ncol = dim( run_tmp )[2] )
    run     <- run_tmp[1:c(dim( run_tmp )[1]-1),]
  }
  
  # 2. Summarise parameters
  cat( "Generate objects with summarised estimated divergence times... ...\n")
  dim.r   <- dim(run)[2]
  # divtimes <- run[,-c(1, dim(run)[2])]
  # If `MCMCtree` has sampled from the prior, then the `lnL` column will not
  # appear and we need to subtract `1` from the `delcol`
  if( prior == TRUE ){
    delcol <- delcol - 1
    divtimes <- run[,-c( 1, ( dim.r-delcol ):dim.r )]
  }else{
    divtimes <- run[,-c( 1, ( dim.r-delcol ):dim.r )]
  }
  
  mean_est_divt  <- apply( X = divtimes, MARGIN = 2, FUN = mean )
  quant_est_divt <- apply( X = divtimes, MARGIN = 2, FUN = quantile, probs = c( 1-perc,perc ) )
  quant_est_divt <- t( quant_est_divt )
  test_alleq     <- all.equal( quant_est_divt[1,], quantile( divtimes[,1], probs = c( 1-perc,perc ) ) )
  if( test_alleq != TRUE ){
    stop( "There was an issue calculating quantiles!" )
  }
  
  # 3. Return object 
  cat( "\nTasks done! Return objects\n\n")
  return( list( divt = divtimes, mean_divt = mean_est_divt, quant_divt = quant_est_divt ) )
  
}


# Function to find problematic runs.
# It prints out the mean div time, qlow, and qup 
# for each node and for each run
# 
# Parameters:
# num_dirs       Numeric. Default is 32 as the number of chains run in a first instance.
# delcol         Numeric, number of columns that are to be deleted as they do not contain
#                divtimes. You can read the header of the `mcmc.txt` files and count the
#                number of `mu*` and `sigma2*` elements. Do not count the `lnL` value as
#                this is automatically ignored in the function `load_dat`. Assuming an MCMC
#                run under a relaxed-clock model with no partitions, we would see `mu` and
#                `sigma2` columns. Therefore, `delcol = 2`. The default value is 2, but you
#                will need to change it if you run `MCMCtree` with different settings.
# data_dir       Character. Path to the directory where the analyses for this datasubset ran 
#                and the `mcmc.txt` and `FigTree` files can be found.
# num_divt       Numeric. Number of columns of the `mcmc.txt` that correspond to the samples 
#                collected for the times.
# node_calib     Character. CSV file with two columns, where the first column has the names 
#                of the calibrations and in the second the corresponding node.
# perc           Numeric. Percentile to calculate the quantiles. Default: 0.975.
# clean          Boolean. FALSE if a "clean" file removing incomplete lines in the `mcmc.txt`
#                has not been created. TRUE otherwise.
# prior          Boolean. TRUE if `MCMCtree` has been sampling from the prior, FALSE
#                otherwise.
# dataset        Character, name of the data subset that will be used to generate the name 
#                of the output directory.
# out_dat        Character, path to the directory to save output data. Do not add the last "/".
find_prob_MCMC <- function ( num_dirs = 32, delcol = 3,
                             data_dir, num_divt, node_calib,
                             perc = 0.975, clean = FALSE, prior = FALSE,
                             out_dat, dataset )
{
  
  # 0. Allow for numbers not using exponentials if too low !
  options( scipen = 999 )
  
  # 1. Create global vars
  if( length( num_dirs ) == 1 ){
    run_dirs   <- num_dirs           # If all dirs are there (i.e., only one num)
    seq_dirs   <- c(1:run_dirs)
  }else{
    run_dirs   <- length( num_dirs ) # If some dirs have been deleted (i.e., vector)
    seq_dirs   <- num_dirs
  }
  total_runs         <- run_dirs
  subtree_list       <- vector( mode = "list", length = total_runs )
  subtree_meandivt   <- matrix( 0, nrow = total_runs, ncol = num_divt )
  subtree_qup        <- matrix( 0, nrow = total_runs, ncol = num_divt )
  subtree_qlow       <- matrix( 0, nrow = total_runs, ncol = num_divt )
  names( subtree_list )        <- rownames( subtree_meandivt ) <-
    rownames( subtree_qup ) <- rownames( subtree_qlow ) <- paste( "run", 1:total_runs, sep = "" )
  
  # 2. Get summarised data
  count <- 0
  count_dim <- vector( "numeric", total_runs )
  #for( i in 1:run_dirs ){
  for( i in seq_dirs ){
    count <- count + 1
    cat( "Loading run ", i, "... ...\n" )
    if( clean == FALSE ){
      subtree_list[[ count ]]  <- load_dat( mcmc = paste( data_dir, i, "/mcmc.txt", sep = "" ),
                                            delcol = delcol, perc = perc, prior = prior )
    }else if( clean == TRUE ){
      subtree_list[[ count ]]  <- load_dat( mcmc = paste( data_dir, i, "/mcmc_clean.txt", sep = "" ),
                                            delcol = delcol, perc = perc, prior = prior )
    }
    
    count_dim[count]         <- dim( subtree_list[[ count ]]$divt )[1]
    subtree_meandivt[count,] <- subtree_list[[ count ]]$mean_divt
    subtree_qlow[count,]     <- subtree_list[[ count ]]$quant_divt[,1]
    subtree_qup[count,]      <- subtree_list[[ count ]]$quant_divt[,2]
  }
  colnames( subtree_meandivt ) <- colnames( subtree_qup ) <-
    colnames( subtree_qlow ) <- rownames( subtree_list[[ 1 ]]$quant_divt )
  
  # 3. Get mean data across runs. Apparently it is more accurate 
  #    than actually using meandivt per run and it matches metrics 
  #    calculated by MCMCtree :|
  mcmc_all <- matrix( 0, nrow = sum( count_dim ), ncol = num_divt )
  colnames( mcmc_all ) <- rownames( subtree_list[[ 1 ]]$quant_divt )
  start <- 1
  stop  <- 0
  cat( "Calculating mean and quantiles for all samples... ...\n\n")
  for( i in 1:total_runs ){
    # cat( "Start: ", start, "\n")
    stop <- stop + count_dim[i]
    # cat( "Stop: ", stop, "\n" )
    mcmc_all[c(start:stop),] <- matrix( unlist( subtree_list[[ i ]]$divt ), nrow = num_divt, byrow = FALSE ) 
    start <- stop + 1
  }
  
  mean_divt    <- apply( X = mcmc_all, MARGIN = 2, FUN = mean )
  mean_quants  <- apply( X = mcmc_all, MARGIN = 2, FUN = quantile, probs = c( 1-perc,perc ) )
  mean_est         <- matrix( 0, nrow = length( mean_divt ), ncol = 3 )
  mean_est_priors  <- matrix( 0, nrow = length( mean_divt ), ncol = 4 )
  mean_est[,1]     <- mean_est_priors[,1] <- round( mean_divt*100, digits = 3 )
  mean_est[,2]     <- mean_est_priors[,2] <- round( as.numeric( format( mean_quants[1,], scientific = FALSE ) )*100, digits = 2 )
  mean_est[,3]     <- mean_est_priors[,3] <- round( as.numeric( format( mean_quants[2,], scientific = FALSE ) )*100, digits = 2 )
  colnames( mean_est )        <- c( "Mean_time", "Mean_qlow", "Mean_qup" )
  colnames( mean_est_priors ) <- c( "Mean_time", "Mean_qlow", "Mean_qup", "Priors" )
  test_names <- all.equal( names( mean_divt ), colnames( mean_quants ) )
  if( test_names != TRUE ){
    stop( "Issue with names for mean divt, mean qup, and mean q low!" )
  }
  rownames( mean_est ) <- rownames( mean_est_priors ) <- names( mean_divt )
  
  # 4. Get matching nodes with calib names
  match_csv <- read.table( node_calib, header = TRUE, sep = ";", stringsAsFactors = FALSE )
  ind_match <- which( as.numeric( gsub( pattern = "t_n", replacement = "", x = rownames( mean_est_priors ) ) )
                      %in% match_csv[,2] )
  for( i in ind_match ){
    tmp_ind <- which( paste( "t_n", match_csv[,2], sep = "" ) == rownames( mean_est_priors )[i] )
    rownames( mean_est_priors )[i] <- paste( "t_n", match_csv[tmp_ind,2], "_", match_csv[tmp_ind,1], sep = "" )
    mean_est_priors[i,4] <- match_csv[tmp_ind,3]
  }
  
  # 5. Write separate output files for mean times, mean qup, and mean qlow for each 
  #    run
  if ( ! dir.exists( paste( out_dat, "/", dataset, sep = "" ) ) ){
    dir.create( paste( out_dat, "/", dataset, sep = "" ) )
  }
  write.table( x = subtree_meandivt, file = paste( out_dat, "/", dataset, "/mean_divt.tsv", sep = "" ),
               sep = "\t", quote = FALSE )
  write.table( x = subtree_qup, file = paste( out_dat, "/", dataset, "/mean_qup.tsv", sep = "" ), sep = "\t",
               quote = FALSE )
  write.table( x = subtree_qlow, file = paste( out_dat, "/", dataset, "/mean_qlow.tsv", sep = "" ), sep = "\t",
               quote = FALSE )
  write.table( x = mean_est_priors,
               file = paste( out_dat, "/", dataset, "/", dataset, "_all_mean_est.tsv", sep = "" ),
               sep = "\t", quote = FALSE )
  
  cat( "Output files available! Check returned list with all objects generated too :) \n\n")
  
  return( list( tt_all = subtree_list, mean = subtree_meandivt, qup = subtree_qup, qdown = subtree_qlow,
                all_mean_est = mean_est, all_mcmc = mcmc_all ) )
  
}


# Plot convergence plot 
# Function to plot a convergence plot 
# Arguments:
#
# name_dir_subt  Character. Name of the directory where the anlayses for this subtree ran.
# mean_divt1     Character. Name of the object that contains the min divtimes for the
#                first half of runs.
# mean_divt2     Character. Name of the object that contains the min divtimes for the
#                second half of runs.
# num_runs       Integer. Number of runs.
plot_convergence <- function ( name_dir_subt, mean_divt1, mean_divt2, num_runs )
{
  half <- round( num_runs/2, digits = 0 )
  tmp <- bquote( paste(  hat( italic(t) ), " | MCMC - run 1-", .(half), sep = "" ) )
  tmp2 <- bquote( paste( hat( italic(t) ), " | MCMC - run ", .(half+1), "-", .(num_runs), sep = "" ) )
  plot( x = mean_divt1, y = mean_divt2,
        main = paste( "Convergence plot - ", name_dir_subt, sep = "" ),
        xlab = tmp, ylab = tmp2 )
  abline( lm( mean_divt2~0 + mean_divt1 ),
          col="red", lty = 2 )
}


# Check quantiles in chains
# Function to check the quantiles calculated based on the estimated divergence times in 
# each chain. It can use either a specific chain against which the rest of the chains are
# compared to or do an "all-against-all" comparison
#
# Arguments
# 
# dat          List, object created with the functions above to store the summary of the
#              chains ran with MCMC. It contains the q2.5% and q97.5% as the third and 
#              fourth entries, which are used in this function.
# threshold    Numeric, threshold over which chains will be evaluated as "problematic".
#              Default: 0.2. NOTE: This might be too restrictive sometimes, so maybe 
#              0.3 could be used too.
# main_chain   Integer, number of the main_chain. Required when a specific chain is to be
#              used against which the rest will be compared.
# num_chains   Integer, number of chains to be analysed for each dataset. Default: 32.
# compare_all  Boolean, TRUE if an all-against-all comparison is to be made. FALSE otherwise.
# outdir       Character, path to where the output file with flagged chains should be output.
check_quantiles <- function( dat, threshold = 0.2, main_chain = 1, num_chains = 32,
                             compare_all = TRUE, outdir )
{
  
  if( compare_all == TRUE ){
    eval_compare <- vector( mode = "list", length = c(num_chains-1) )
    eval_mat <- matrix( 0, ncol = 2, nrow = c(num_chains-1) )
    colnames( eval_mat )   <- c( "Prob_chain", "Prob_node/s" )
    names( eval_compare )  <- paste( "main_chain_", 1:c(num_chains-1), sep = "" )
    sum_chains             <- matrix( 0, nrow = num_chains, ncol = 2 )
    colnames( sum_chains ) <- c( "Num_prob_chain/s", "Main_chain" )
    for ( k in 1:num_chains ){
      #cat( "Main chain compared is: ", k, "\n" )
      rownames( eval_mat ) <- rep( paste( "main_chain_", k, sep = "" ), c(num_chains-1) )
      tmp_ch       <- c( 1:num_chains ) 
      ind_ch       <- which( tmp_ch %in% k )
      eval_ch      <- tmp_ch[-ind_ch]
      count        <- 0
      count_row    <- 0
      for ( j in eval_ch ){
        count_row <- count_row + 1
        tmp_prob_chains <- rep( 0, 2 )
        tmp_prob_nodes  <- rep( NA, 2 )
        tmp_qup         <- dat[[3]][k,] - dat[[3]][j,]
        tmp_ind_qup     <- which( abs( tmp_qup ) > threshold )
        if( length( tmp_ind_qup ) > 0 ){
          #cat( "q97.5%: Check the following nodes in chain ", j, ":", names( tmp_qup[tmp_ind_qup] ), "\n" )
          #cat(  "Difference: ", tmp_qup[tmp_ind_qup], "\n")
          count <- count + 1
          tmp_prob_chains[1] <- j
          tmp_prob_nodes[1] <- paste0( names( tmp_qup[tmp_ind_qup] ), collapse = "-" )
        }
        tmp_qdown    <- dat[[4]][k,] - dat[[4]][j,]
        tmp_ind_qdown<- which( abs( tmp_qdown ) > threshold )
        if( length( tmp_ind_qdown ) > 0 ){
          #cat( "q2.5%: Check the following nodes in chain ", j, ":", names( tmp_qdown[tmp_ind_qdown] ), "\n" )
          #cat(  "Difference: ", tmp_qdown[tmp_ind_qdown], "\n")
          count <- count + 1
          tmp_prob_chains[2] <- j
          tmp_prob_nodes[2]  <- paste0( names( tmp_qdown[tmp_ind_qdown] ), collapse = "-" )
        }
        # If there have been no problematic chains found, then the first condition is TRUE
        if( all.equal( tmp_prob_chains, c( 0, 0) ) == TRUE ){
          eval_mat[count_row,1]   <- NA
          eval_mat[count_row,2]   <- NA
        }else{
          # Otherwise, just find the problematic chain and nodes and append them to the
          # vector
          tmp_prob_chains <- unique( tmp_prob_chains )
          ind_rm          <- which( tmp_prob_chains == 0 )
          ind_rm2         <- which( tmp_prob_nodes %in% NA )
          if( length( ind_rm ) > 0 ){
            tmp_prob_chains <- tmp_prob_chains[-ind_rm]
          }
          if( length( ind_rm2 ) > 0 ){
            tmp_prob_nodes  <- tmp_prob_nodes[-ind_rm2]
          }
          tmp_prob_nodes  <- paste0( unique( tmp_prob_nodes ), collapse = "-" )
          eval_mat[count_row,1]   <- tmp_prob_chains
          eval_mat[count_row,2]   <- tmp_prob_nodes
        }
        
      }
      
      # Check how many chains 
      main_ch <- c( max(eval_ch), k )
      tmp_all <- eval_ch
      rm_ind  <- which( eval_mat[,1] %in% NA )
      if( length( rm_ind ) > 0 ){
        cat( "Main chain ", k, "-- number of problematic chains : ", length( tmp_all[-rm_ind] ), "\n" )
        if ( length( tmp_all[-rm_ind] ) < main_ch[1] ){
          main_ch[1] <- c( length( tmp_all[-rm_ind] ) )
        }
      }else{
        cat( "Chain ", k, "is too problematic when used as main! \n" )
        main_ch[1] <- num_chains-1
      }
      
      sum_chains[k,]      <- main_ch
      eval_compare[[ k ]] <- eval_mat
    }
    
    # Return object
    return( list( eval_obj = eval_compare, sum_chains = sum_chains ) )
    
  }else{
    tmp_ch  <- c( 1:num_chains ) 
    ind_ch  <- which( tmp_ch %in% main_chain )
    eval_ch <- tmp_ch[-ind_ch]
    for ( j in eval_ch ){
      tmp_qup     <- dat[[3]][ind_ch,] - dat[[3]][j,]
      tmp_ind_qup <- which( abs( tmp_qup ) > threshold )
      if( length( tmp_ind_qup ) > 0 ){
        cat( "q97.5%: You will have to check some nodes in chain ", j, "\n" )
        #cat(  "Difference: ", tmp_qup[tmp_ind_qup], "\n")
        write( x = paste( "q97.5%: Check the following nodes in chain ", j, "\nNodes: ",
                          paste0(names( tmp_qup[tmp_ind_qup] ), collapse = " | " ), "\n", 
                          "Difference: ", paste0(tmp_qup[tmp_ind_qup], collapse = " | "),
                          "\n", sep = "" ),
               file = paste( outdir, "check_chains.txt", sep = "" ), append = TRUE )
      }
      tmp_qdown    <- dat[[4]][ind_ch,] - dat[[4]][j,]
      tmp_ind_qdown<- which( abs( tmp_qdown ) > threshold )
      if( length( tmp_ind_qdown ) > 0 ){
        cat( "q2.5%: You will have to check some nodes in chain ", j, "\n" )
        #cat(  "Difference: ", tmp_qdown[tmp_ind_qdown], "\n")
        write( x = paste( "q2.5%: Check the following nodes in chain ", j, "\nNodes: ",
                          paste0(names( tmp_qdown[tmp_ind_qdown] ), collapse = " | "), "\n",
                          "Difference: ", paste0(tmp_qdown[tmp_ind_qdown], collapse = " | "),
                          "\n", sep = "" ),
               file = paste( outdir, "check_chains.txt", sep = "" ),
               append = TRUE )
      }
    }
    
    # Informative message
    cat( "\n>> If any conflicts have been found, you will find file",
         paste( outdir, "check_chains.txt", sep = "" ), ".\n",
         "Please open this file to check conflictive nodes if created.\n")
  }
  
} 