#--------------------#
# DEFINE GLOBAL FUNS #
#--------------------#

# Load packages required in functions
library( sn )
library( rstan )
library( stringr )

# Function to set the main working directories
# No arguments needed
#
set_homedir <- function( ){
  # Get the path to current open R script and find main dir "scripts"
  path_to_file <- getActiveDocumentContext()$path
  wd <- paste( dirname( path_to_file ), "/", sep = "" )
  #wd <- gsub( pattern = "scripts/", replacement = "ESS_and_chains_convergence/", x = wd )
  wd <- gsub( pattern = "scripts/", replacement = "plots/", x = wd )
  if( ! dir.exists( wd ) ){
    dir.create( wd )
  }
  if( ! dir.exists( paste( wd, "ESS_and_chains_convergence", sep = "" ) ) ){
    dir.create( paste( wd, "ESS_and_chains_convergence", sep = "" )  )
  }
  wd2 <- gsub( pattern = "plots/", replacement = "", x = wd )
  # if( ! dir.exists( paste( wd, "plots", sep = "" ) ) ){
  #   dir.create( paste( wd, "plots", sep = "" )  )
  # }
  # wd2 <- gsub( pattern = "ESS_and_chains_convergence/", replacement = "", x = wd )
  return( list( home = wd2, ESS = wd ) )
}

# Function to load the `mcmc.txt` files.
# Used to get mean estimates for convergence plots.
# 
# Parameters:
# mcmc         Character, path to `mcmc.txt` file.
# delcol       Numeric, number of columns that are to be deleted 
#              as they do not contain divtimes. You can read the header
#              of the `mcmc.txt` files and count the number of `mu*` and
#             `sigma2*` elements. Do not count the `lnL` value as this is
#              automatically deleted in the calculations below. 
#              Assuming an MCMC run under a relaxed-clock model with no 
#              partitions, we would see `mu` and `sigma2` columns. Therefore,
#              `delcol = 2`. The default value is 2, but you will need to change
#              it if you run `MCMCtree` with different settings.
# perc         Numeric. Percentile to calculate the quantiles. Default: 0.975.
# def_samples  Numeric. Number of samples that the user defined through the
#              `MCMCtree` option `nsample. 
# prior        Boolean. TRUE if `MCMCtree` has been sampling from the prior, FALSE
#              otherwise
load_dat <- function( mcmc, delcol = 3, perc = 0.975, def_samples = 20000,
                      prior = FALSE )
{
  
  # 1. Load files and get parameters. Load everything with option `fill=TRUE`
  # in case a line is not complete
  cat( "Load \"mcmc.txt\" from path \n", mcmc, "\n... ...\n" )
  run_tmp  <- read.table( mcmc, header = TRUE, fill = TRUE, sep = "\t" )
  test_len <- def_samples+1
  #print( test_len )
  #print( dim(run_tmp))
  if( test_len != dim(run_tmp)[1] ){
    cat( " [[ The last line of the \"mcmc.txt\" is incomplete and will be removed ]] \n\n" )
    run <- matrix( 0, nrow = c(dim( run_tmp )[1]-1), ncol = dim( run_tmp )[2] )
    run <- run_tmp[1:c(dim( run_tmp )[1]-1),]
  }else{
    cat( " [[ MCMCtree collected the same amount of samples you specified ]] \n",
         " [[       All lines will be kept in the \"mcmc.txt\"            ]] \n\n" )
    run <- run_tmp
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
  # Calculate mean and quantiles
  mean_est_divt  <- apply( X = divtimes, MARGIN = 2, FUN = mean )
  quant_est_divt <- apply( X = divtimes, MARGIN = 2, FUN = quantile, probs = c( 1-perc,perc ) )
  quant_est_divt <- t( quant_est_divt )
  test_alleq     <- all.equal( quant_est_divt[1,], quantile( divtimes[,1], probs = c( 1-perc,perc ) ) )
  if( test_alleq != TRUE ){
    stop( "There was an issue calculating quantiles!" )
  }
  
  # 3. Return object 
  cat( "\nTasks done! Returning objects\n\n")
  return( list( divt = divtimes, mean_divt = mean_est_divt, quant_divt = quant_est_divt ) )
  
}


# Function to read the csv file/s (semicolon separated) with the information
# about the calibrations used for the input tree in `MCMCtree`.
#
# Arguments>
#
# f_dir       Character, (absolute or relative) path to the directory where the
#             file/s are saved. Please type a "/" at the end of the path.
# f_names     Character, vector with the file name/s. If more than one file, please
#             separate the names with comas as in a character vector.
# dat         Character, global vector created at the beginning of he script 
#             with the name that users have given to each analysis. The length of 
#             this vector is equivalent to the number of `MCMCtree` jobs that have
#             been run for each hypothesis (e.g., different tree hypotheses, 
#             different calibration hypotheses, etc.). The length of this vector 
#             must be equal to the length of the character vector passed to 
#             argument `f_names`.
# head_avail  Boolean, TRUE if the header is available in the input files.
#             FALSE otherwise.
#
read_calib_f <- function( main_dir, f_names, dat, head_avail = TRUE )
{
  
  # Check length(dat) == length(f_names), otherwise abort
  if( length( dat ) != length( f_names ) ){
    stop( paste( "\nThe length of character vector passed to argument \"dat\" must be",
                 " the same length as of the character vector passed to argument \"f_names\"\n",
                 sep = "" ) ) 
  }
  # Generated list vector to save files
  calib_nodes          <- vector( mode = "list", length = length( dat ) )
  names( calib_nodes ) <- dat
  # Now, load the file/s on the list object `calib_nodes`
  for( i in 1:length(dat) ){
    
    if( head_avail == TRUE ){
      calib_nodes[[ i ]] <- read.table( file = paste( main_dir, f_names[i], sep = "" ),
                                        header = TRUE, sep = ";",
                                        stringsAsFactors = FALSE )
    }else{
      calib_nodes[[ i ]] <- read.table( file = paste( main_dir, f_names[i], sep = "" ),
                                        header = FALSE, sep = ";",
                                        stringsAsFactors = FALSE )
      colnames( calib_nodes ) <- c( "Calib", "node", "Prior" )
    }
    
  }
  
  # Return final object
  return( calib_nodes )
  
}

# Function to summarise runs.
# It prints out the mean div time, qlow, and qup 
# for each node and for each run
# 
# Parameters:
# num_dirs       Numeric. Default is 36 as the number of chains run in a first
#                instance.
# delcol         Numeric, number of columns that are to be deleted as they do
#                not contain divtimes. You can read the header of the `mcmc.txt`
#                files and count the number of `mu*` and `sigma2*` elements. Do
#                not count the `lnL` value as this is automatically ignored in
#                the function `load_dat`. Assuming an MCMC run under a
#                relaxed-clock model with no partitions, we would see `mu` and
#                `sigma2` columns. Therefore, `delcol = 2`. The default value is
#                2, but you will need to change it if you run `MCMCtree` with
#                different settings.
# data_dir       Character. Path to the directory where the analyses for this
#                dataset ran and the `mcmc.txt` file can be found.
# num_divt       Numeric. Number of columns of the `mcmc.txt` that correspond
#                to the samples collected for the times.
# node_calib     Character. CSV file with two columns, where the first column
#                has the names of the calibrations and in the second the
#                corresponding node.
# perc           Numeric. Percentile to calculate the quantiles. Default: 0.975.
# def_samples    Numeric. Number of samples that the user defined through the
#                `MCMCtree` option `nsample. 
# prior          Boolean. TRUE if `MCMCtree` has been sampling from the prior,
#                FALSE otherwise.
# dataset        Character, name of the dataset that will be used to  
#                generate the name of the output directory.
# out_dat        Character, path to the directory to save output data. Do not
#                add the last "/".
# time_unit      Integer, value that you need to multiply the estimated divtimes
#                with
sum_MCMC <- function ( num_dirs = 36, delcol = 3,
                       data_dir, num_divt, node_calib,
                       perc = 0.975, def_samples = 20000, prior = FALSE,
                       out_dat, dataset, time_unit = 100 )
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
  for( i in seq_dirs ){
    count <- count + 1
    cat( "Loading run ", i, "... ...\n" )
    subtree_list[[ count ]]  <- load_dat( mcmc = paste( data_dir, i, "/mcmc.txt", sep = "" ), 
                                          delcol = delcol, perc = perc, def_samples = 20000,
                                          prior = prior )
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
  # Generate empty array that fits Rstan format
  # dim = iterations/samples, num chains, num parameters
  # NOTE: I am adding 1 to `def_samples` value just in case all samples have 
  # been collected!
  arr_all <- array( 0, dim = c( def_samples+1, run_dirs, num_divt ),
                    dimnames = list( samples = paste( "samp", 1:c(def_samples+1), sep = "" ),
                                     chains = paste( "chain", 1:run_dirs, sep = "" ),
                                     parameters = paste( "t_n",c(c(num_divt+2):c(num_divt+num_divt+1)), sep = "" )
                    ) )
  for( i in 1:total_runs ){
    #cat( "Start: ", start, "\n")
    stop <- stop + count_dim[i]
    #cat( "Stop: ", stop, "\n" )
    mcmc_all[c(start:stop),] <- matrix( unlist( subtree_list[[ i ]]$divt ), nrow = num_divt, byrow = FALSE ) 
    tmp_rows <- dim( mcmc_all[c(start:stop),] )[1]
    tmp_cols <- dim( mcmc_all[c(start:stop),] )[2]
    # Now, add samples to each row in the array generated
    for( p in 1:tmp_cols ){
      arr_all[1:tmp_rows,i,p] <- mcmc_all[c(start:stop),][,p]
    }
    start <- stop + 1
  }
  
  mean_divt    <- apply( X = mcmc_all, MARGIN = 2, FUN = mean )
  mean_quants  <- apply( X = mcmc_all, MARGIN = 2, FUN = quantile, probs = c( 1-perc,perc ) )
  mean_est         <- matrix( 0, nrow = length( mean_divt ), ncol = 3 )
  mean_est_priors  <- matrix( 0, nrow = length( mean_divt ), ncol = 4 )
  mean_est[,1]     <- mean_est_priors[,1] <- round( mean_divt*time_unit, digits = 3 )
  mean_est[,2]     <- mean_est_priors[,2] <- round( as.numeric( format( mean_quants[1,], scientific = FALSE ) )*time_unit, digits = 2 )
  mean_est[,3]     <- mean_est_priors[,3] <- round( as.numeric( format( mean_quants[2,], scientific = FALSE ) )*time_unit, digits = 2 )
  colnames( mean_est )        <- c( "Mean_time", "Mean_qlow", "Mean_qup" )
  colnames( mean_est_priors ) <- c( "Mean_time", "Mean_qlow", "Mean_qup", "Priors" )
  test_names <- all.equal( names( mean_divt ), colnames( mean_quants ) )
  if( test_names != TRUE ){
    stop( "Issue with names for mean divt, mean qup, and mean q low!" )
  }
  rownames( mean_est ) <- rownames( mean_est_priors ) <- names( mean_divt )
  
  # 4. Get matching nodes with calib names. Updated 231030: `node_calib` is 
  # the csv file formatted into a matrix, not the list!
  #match_csv <- node_calib[[1]]
  match_csv <- node_calib
  ind_match <- which( as.numeric( gsub( pattern = "t_n", replacement = "",
                                        x = rownames( mean_est_priors ) ) )
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
                all_mean_est = mean_est, all_mcmc = mcmc_all, arr4stan = arr_all,
                samp_per_chain = count_dim ) )
  
}


# Plot convergence plot 
# Function to plot a convergence plot 
# Arguments:
#
# name_dir_dat  Character. Name of the directory where the anlayses for this subtree ran.
# mean_divt1    Character. Name of the object that contains the min divtimes for the
#               first half of runs.
# mean_divt2    Character. Name of the object that contains the min divtimes for the
#               second half of runs.
# num_runs      Integer. Number of runs.
plot_convergence <- function ( name_dir_dat, mean_divt1, mean_divt2, num_runs )
  {
  half <- round( num_runs/2, digits = 0 )
  tmp <- bquote( paste(  hat( italic(t) ), " | MCMC - run 1-", .(half), sep = "" ) )
  tmp2 <- bquote( paste( hat( italic(t) ), " | MCMC - run ", .(half+1), "-", .(num_runs), sep = "" ) )
  plot( x = mean_divt1, y = mean_divt2,
        main = paste( "Convergence plot - ", name_dir_dat, sep = "" ),
        xlab = tmp, ylab = tmp2 )
  abline( lm( mean_divt2~0 + mean_divt1 ),
          col="red", lty = 2 )
}

# Plot traces 
# Function to generate a trace plot
# Arguments:
#
# name_dir_dat  Character. Name of the directory where the analyses ran.
# sum_dat       List. Object generated with function `sum_MCMC`.
# divt          Numeric. Node labels for which you want to get a trace plot.
# n_chains      Integer. Number of chains for which you want to get the trace 
#               plot.
plot_traces <- function ( name_dir_dat, sum_dat, divt, n_chains, out_dir )
{
  if( length( n_chains ) == 1 ){
    num_chains <- c(1:n_chains)
  }else{
    num_chains <- n_chains
  }
  for( i in num_chains ){
    if( ! dir.exists( paste( out_dir, i, "/traceplots/", sep = "" ) ) ){
      dir.create( paste( out_dir, i, "/traceplots/", sep = "" ) )
    }
    # sum_dat[[1]] --> accesses object `tt_all`
    # sum_dat[[1]][[i]] --> accesses object `tt_all` for chain `i`
    # sum_dat[[1]][[i]][[1]] --> accesses object `tt_all` for chain `i` and
    # gets the `divt`!
    # sum_dat[[1]][[i]][[1]][,tmp_node2] --> accesses object `tt_all` for chain `i`,
    # gets the `divt` only for node `j`!
    for( j in divt ){
      tmp <- "Iterations"
      tmp2 <- bquote( paste( hat( italic(t) ), " - node ", .(j),
                             " (time unit = 100 Myr)", sep = "" ) )
      tmp_node <- gsub( pattern = "t_n", replacement = "",
                        x = colnames(sum_dat[[1]][[i]][[1]] ) )
      tmp_node2 <- which( tmp_node %in% j )
      pdf( paste( out_dir, i, "/traceplots/traceplot_", name_dir_dat, "_n", j, 
                  "_chain", i, ".pdf", sep = "" ),
           paper = "a4" )
      plot( x = 1:length(sum_dat[[1]][[i]][[1]][,tmp_node2]),
            y = sum_dat[[1]][[i]][[1]][,tmp_node2], type = "l",
            main = paste( "Trace plot - ", name_dir_dat, sep = "" ),
            xlab = tmp, ylab = tmp2 )
      dev.off()
    }
  }
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
    # Create file, will be deleted if no chains to be checked!
    write( x = c( "Checking problematic nodes\n"),
           file = paste( outdir, "check_chains.txt", sep = "" ) )
    visited <- 0
    for ( j in eval_ch ){
      tmp_qup     <- dat[[3]][ind_ch,] - dat[[3]][j,]
      tmp_ind_qup <- which( abs( tmp_qup ) > threshold )
      if( length( tmp_ind_qup ) > 0 ){
        visited <- visited + 1
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
        visited <- visited + 1
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
    
    # Check whether the out file is to be kept or deleted
    if( visited == 0 ){
      unlink( paste( outdir, "check_chains.txt", sep = "" ) )
    }
    
    # Informative message
    cat( "\n>> If any conflicts have been found, the following file should have been created:\n",
         paste( outdir, "check_chains.txt", sep = "" ),"\n",
         "If so, please open this file to check conflictive nodes.\n")
  }
  
} 


# Function to run `rstan::monitor` and export the 
# median, min, and max values for the tail-ESS and the 
# bulk-ESS
#
# Parameters:
# x              Array, object generated by function `sum_MCMC`, part of the 
#                list that said function returns. This array is already in the 
#                format required for the stand option.
#
# samp_per_chain Numeric, object generated by function `sum_MCMC`, part of the  
#                list that said function returns. This vector contains the 
#                number of samples collected for each of the chains.
sum_MCMC_ESS <- function( x, samp_per_chain )
{
  
  # Define empty vector to report stats
  out_stats <- matrix( 0, nrow = 3, ncol = 2 )
  colnames( out_stats ) <- c( "Tail-ESS", "Bulk-ESS" )
  rownames( out_stats ) <- c( "Med.", "Min.", "Max." )
  
  # Compute stats
  min_samp       <- min( samp_per_chain )
  tmp_arr        <- x[c(1:min_samp),,]
  stats_mcmc     <- rstan::monitor( sims = tmp_arr )
  out_stats[1,1] <- median( stats_mcmc$Tail_ESS )
  out_stats[2,1] <- min( stats_mcmc$Tail_ESS )
  out_stats[3,1] <- max( stats_mcmc$Tail_ESS )
  out_stats[1,2] <- median( stats_mcmc$Bulk_ESS )
  out_stats[2,2] <- min( stats_mcmc$Bulk_ESS )
  out_stats[3,2] <- max( stats_mcmc$Bulk_ESS )
  
  # Return stats 
  return( list( tab = out_stats, stats = stats_mcmc ) )
  
}


# Function to plot the distribution for soft bounds that is implemented
# in MCMCtree following Eq. 17 in Yang and Rannala 2005.
# Paper: https://academic.oup.com/mbe/article/23/1/212/1193630
#
# Parameters:
#
# t    Numeric, range of time values for which the distribution wants to be plotted.
# tL   Numeric, minimum age the user has specified for the soft-bound distribution.
# tU   Numeric, maximum age the user has specified for the soft-bound distribution.
softbounds <- function( t, tL, tU, pL = 0.025, pU = 0.025 )
{
  
  # If pL = pU = 0.025, then pmax = 0.95
  pmax  <- 1 - pL - pU
  # For minimum bound, power decay area
  theta1 <- pmax*tL/(pL*(tU-tL))
  # For maximum bound, exponential decay area
  theta2 <- pmax/(pU*(tU-tL))
  # Save all values in this object
  all    <- rep( 0, length(t) )
  
  count <- 0
  for ( i in t ){
    count <- count + 1
    if ( i > 0 && i <= tL ){
      # Power decay on the left side with probability `pL`
      c1 <- pL * (theta1/tL) * ( i/tL)^(theta1-1)
      all[count] <- c1
    }
    else if ( i > tL && i <=  tU ){
      # Uniform distribution with probability `pmax`
      c2 <- pmax/(tU-tL)
      all[count] <- c2
    }
    else if ( i > tU ){
      # Exponential decay on the right side with probability `pU`
      c3 <- pU*theta2*exp( -theta2*(i-tU) ) # for the maximum tail
      all[count] <- c3
    }
  }
  
  return( all )
}


# Function to plot the distribution for lower bounds that is implemented
# in MCMCtree following Eq. 25 in Inoue et al. 2010.
# Paper: https://academic.oup.com/sysbio/article/59/1/74/1725802?login=true
#
# Parameters:
#
# t    Numeric, range of time values for which the distribution wants to be plotted.
# tL   Numeric, minimum age the user has specified for the lower bound.
# p    Numeric, offset proportion. Default: 0.1.
# c    Numeric, scale parameter. Default: 1.
# pL   Numeric, left tail probability. Default: 0.025.
lowerbounds <- function( t, tL, p = 0.1, c = 1, pL = 0.025 )
{
  
  # If pL = 0.025, then pmax = 0.975
  pmax <- 1 - pL
  # Normalising constant due to truncation of the Cauchy distribution
  A      <- 1/2 + 1/pi * atan(p/c)
  # Theta that allows density at tL to be continuous
  theta  <- pmax/pL * ( 1 / ( pi * A * c * ( 1 + (p/c)^2 ) ) )
  # Vector to save output vals of t
  all    <- rep( 0, length(t) )
  
  count <- 0
  for ( i in t ){
    count <- count + 1
    if ( i > 0 && i <= tL ){
      c1 <- pL * (theta/tL) * (i/tL)^(theta-1)
      all[count] <- c1
    }
    else if ( i > tL ){
      d1_c2 <- A * pi * c * tL 
      d2_c2 <- ( ( i - tL * (1+p) ) / ( c * tL ) )^2 
      c2    <- pmax * ( 1 / ( d1_c2 * ( 1 + d2_c2 ) ) )
      all[count] <- c2
    }
  }
  
  return( all )
}


# Function to plot the distribution for upper bounds that is implemented
# in MCMCtree following Eq. 16 in Yang and Rannala 2005.
# Paper: https://academic.oup.com/mbe/article/23/1/212/1193630
#
# Parameters:
#
# t    Numeric, range of time values for which the distribution wants to be plotted.
# tU   Numeric, maximum age the user has specified for the upper bound.
# pU   Numeric, right tail probability. Default: 0.025.
upperbounds <- function( t, tU, pU = 0.025 )
{
  
  # If pL = 0.025, then pmax = 0.975
  pmax <- 1 - pU
  # Theta that allows density at tL to be continuous
  theta  <- pmax/(pU*tU)
  # Vector to save output vals of t
  all    <- rep( 0, length(t) )
  
  count <- 0
  for ( i in t ){
    count <- count + 1
    if ( i < tU ){
      c1 <- pmax/tU
      all[count] <- c1
    }
    else if ( i >= tU ){
      c2    <- pU * theta * exp(-theta*(i-tU))
      all[count] <- c2
    }
  }
  
  return( all )
}



# Function to read the csv file/s (semicolon separated) with the information
# about the calibrations used for the input tree in `MCMCtree`.
#
# Arguments>
#
# f_dir       Character, (absolute or relative) path to the directory where the
#             file/s are saved. Please type a "/" at the end of the path.
# f_names     Character, vector with the file name/s. If more than one file, please
#             separate the names with comas as in a character vector.
# dat         Character, global vector created at the beginning of the script 
#             with the name that users have given to each analysis. The length of 
#             this vector is equivalent to the number of `MCMCtree` jobs that have
#             been run for each hypothesis (e.g., different tree hypotheses, 
#             different calibration hypotheses, etc.). The length of this vector 
#             must be equal to the length of the character vector passed to 
#             argument `f_names`.
# head_avail  Boolean, TRUE if the header is available in the input files.
#             FALSE otherwise.
#
read_calib_f <- function( main_dir, f_names, dat, head_avail = TRUE )
{
  
  # Check length(dat) == length(f_names), otherwise abort
  if( length( dat ) != length( f_names ) ){
    stop( paste( "\nThe length of character vector passed to argument \"dat\" must be",
                 " the same length as of the character vector passed to argument \"f_names\"\n",
                 sep = "" ) ) 
  }
  # Generated list vector to save files
  calib_nodes          <- vector( mode = "list", length = length( dat ) )
  names( calib_nodes ) <- dat
  # Now, load the file/s on the list object `calib_nodes`
  for( i in 1:length(dat) ){
    
    if( head_avail == TRUE ){
      calib_nodes[[ i ]] <- read.table( file = paste( main_dir, f_names[i], sep = "" ),
                                        header = TRUE, sep = ";",
                                        stringsAsFactors = FALSE )
    }else{
      calib_nodes[[ i ]] <- read.table( file = paste( main_dir, f_names[i], sep = "" ),
                                        header = FALSE, sep = ";",
                                        stringsAsFactors = FALSE )
      colnames( calib_nodes ) <- c( "Calib", "node", "Prior" )
    }
    
  }
  
  # Return final object
  return( calib_nodes )
  
}

# Function to plot the calibration density VS the marginal density used by `MCMCtree`
#
# Parameters
#
# calibs     Data frame, object in which the csv file with the calibrations info is saved.
# divt_list  List, object in which the divtimes from `mcmc.txt` have been saved.
# dat        Character, vector with the name identifying the calibration file/s.
#            If more than one, just add as many labels as files you used.
# main_wd    Character, path to the main working directory.
# ind        Boolean, TRUE if individual plots for each calibrated node are to be created.
#            False otherwise.
# clock      Character, type of relaxed-clock model used when running `MCMCtree`.
# out        Character, vector with the name of the dataset/s you are 
#            evaluating.
plot_check_calibnodes <- function( calibs, divt_list, dat, main_wd,
                                   ind = TRUE, clock, n_col, n_row, out )
  
{
  
  # 0. Prepare output file
  if( ! dir.exists( paste( main_wd, "plots", sep = "" ) ) ){
    dir.create( paste( main_wd, "plots", sep = "" ) )
  }
  if( ! dir.exists( paste( main_wd, "plots/calVSmarg/", dat,  sep = "" ) ) ){
    dir.create( paste( main_wd, "plots/calVSmarg/", dat, sep = "" ) )
  }
  
  # In case users have not passed a string of length 3, just get the first three
  # letters
  tmp_dat_dir <- substr( x = dat, start = 1, stop = 3 )
  png( filename = paste( main_wd, "plots/calVSmarg/", dat, "/", out, "_",
                         tmp_dat_dir, "_", clock,".jpg", sep = "" ),
         width = 1024, height = 768 )
  
  par( mfrow = c(n_row,n_col), mai=c(0.1,0.1,0.1,0.1))
  for( i in 1:length( calibs[,2] ) ){
    
    # 1. Find matching node in calib nodes, get min and max vals, and then get values 
    #    of y for times from 0 to 6 (tU = 100My) for the soft bound distribution
    tmp_node <- which( colnames( divt_list[[ 1 ]] ) %in% 
                         paste( "t_n", calibs[i,2], sep = "" ) )
    # Check whether the calibration used was skew-T or soft-bound
    is_ST   <- grep( pattern = "ST", x = calibs[i,3])
    is_B    <- grep( pattern = "B", x = calibs[i,3])
    is_L    <- grep( pattern = "L", x = calibs[i,3])
    is_U    <- grep( pattern = "U", x = calibs[i,3])
    is_SN   <- grep( pattern = "SN", x = calibs[i,3])
    if( length( is_B ) == 1 ){
      cat( paste( "Node ", "t_n", calibs[i,2], " is calibrated with a soft-bound: ",
                  calibs[i,3], "\n", sep = "" ) )
      # NOTE: That's a nice wrapper, but users with R < 4.1 do not have `str_split_.
      # I have decided to then use `str_split`
      # tmp_allB <- as.numeric( stringr::str_split_1( string = gsub( pattern = "B\\(|\\)", replacement = "",
      #                                                              x = calibs[i,3] ),
      #                                               pattern = "," ) )
      tmp_allB <- as.numeric( stringr::str_split( string = gsub( pattern = "B\\(|\\)", replacement = "",
                                                                 x = calibs[i,3] ),
                                                  pattern = "," )[[1]] )
      #print(tmp_allB)
      if( length( tmp_allB ) == 2 ){
        tmp_min <- tmp_allB[1]
        tmp_max <- tmp_allB[2]
        tmp_pL  <- 0.025
        tmp_pU  <- 0.025
      }else if( length( tmp_allB ) == 4 ){
        tmp_min <- tmp_allB[1]
        tmp_max <- tmp_allB[2]
        tmp_pL  <- tmp_allB[3]
        tmp_pU  <- tmp_allB[4]
      }
      tmp_x   <- seq( tmp_min-0.5, tmp_max+0.5, 0.001 )
      tmp_y   <- softbounds( t = tmp_x, tL = tmp_min, tU = tmp_max, pL = tmp_pL, pU = tmp_pU )
      tmp_labcab <- substr( x = calibs[i,1], start = 1, stop = 12 )
    }
    if( length( is_ST ) == 1 ){
      cat( paste( "Node ", "t_n", calibs[i,2], " is calibrated with a skew-t dist: ",
                  calibs[i,3], "\n", sep = "" ) )
      tmp_ST_pars <- as.numeric( stringr::str_split( string = gsub( pattern = "ST\\(|\\)", replacement = "",
                                                                    x = calibs[i,3] ),
                                                     pattern = "," )[[1]] )
      tmp_xi    <- tmp_ST_pars[1] # location (mean age)
      tmp_omega <- tmp_ST_pars[2] # scale 
      tmp_alpha <- tmp_ST_pars[3] # shape
      tmp_nu    <- tmp_ST_pars[4] # df
      
      ext_len     <- 0.1 # Adjust according to the first time you plot -- it will help to
      # properly adjust automatically xlim and ylim
      tmp_x    <- round( max( density( sn::rst( n = 1000, dp = tmp_ST_pars ), adj = 1 )$x ) + 0.5, digits = 1 )
      tmp_y    <- round( max( density( sn::rst( n = 1000, dp = tmp_ST_pars ), adj = 1 )$y ), digits = 1 )
      tmp_labcab <- substr( x = calibs[i,1], start = 1, stop = 12 )
    }
    if( length( is_SN ) == 1 ){
      cat( paste( "Node ", "t_n", calibs[i,2], " is calibrated with a skew-normal dist: ",
                  calibs[i,3], "\n",  sep = "" ) )
      tmp_SN_pars <- as.numeric( stringr::str_split( string = gsub( pattern = "SN\\(|\\)", replacement = "",
                                                                    x = calibs[i,3] ),
                                                     pattern = "," )[[1]] )
      tmp_xi    <- tmp_SN_pars[1] # location (mean age)
      tmp_omega <- tmp_SN_pars[2] # scale 
      tmp_alpha <- tmp_SN_pars[3] # shape
      
      ext_len     <- 0.1 # Adjust according to the first time you plot -- it will help to
      # properly adjust automatically xlim and ylim
      tmp_x    <- round( max( density( sn::rsn( n = 1000, dp = tmp_SN_pars ), adj = 1 )$x ) + 0.5, digits = 1 )
      tmp_y    <- round( max( density( sn::rsn( n = 1000, dp = tmp_SN_pars ), adj = 1 )$y ), digits = 1 )
      tmp_labcab <- substr( x = calibs[i,1], start = 1, stop = 12 )
    }
    if( length( is_L ) == 1 ){
      cat( paste( "Node ", "t_n", calibs[i,2], " is calibrated with a lower bound: ",
                  calibs[i,3], "\n",  sep = "" ) )
      tmp_L <- as.numeric( stringr::str_split( string = gsub( pattern = "L\\(|\\)", replacement = "",
                                                              x = calibs[i,3] ),
                                               pattern = "," )[[1]] )
      if( length( tmp_L ) == 1 ){
        tmp_min <- tmp_L[1]
        tmp_p  <- 0.1
        tmp_c  <- 1
        tmp_pL <- 0.025
      }else if( length( tmp_L ) == 4 ){
        tmp_min <- tmp_L[1]
        tmp_p <- tmp_L[2]
        tmp_c  <- tmp_L[3]
        tmp_pL  <- tmp_L[4]
      }
      tmp_x   <- seq( tmp_min-1, tmp_min+2, 0.001 )
      tmp_y   <- lowerbounds( t=tmp_x, tL=tmp_min, p = tmp_p, c = tmp_c, pL = tmp_pL )
      tmp_labcab <- substr( x = calibs[i,1], start = 1, stop = 12 )
    }
    if( length( is_U ) == 1 ){
      cat( paste( "Node ", "t_n", calibs[i,2], " is calibrated with an upper bound: ",
                  calibs[i,3], "\n",  sep = "" ) )
      tmp_U <- as.numeric( stringr::str_split( string = gsub( pattern = "U\\(|\\)", replacement = "",
                                                              x = calibs[i,3] ),
                                               pattern = "," )[[1]] )
      if( length( tmp_U ) == 1 ){
        tmp_max <- tmp_U[1]
        tmp_pU  <- 0.025
      }else if( length( tmp_U ) == 2 ){
        tmp_max <- tmp_U[1]
        tmp_pU  <- tmp_U[2]
      }
      tmp_x   <- seq( tmp_max-2, tmp_max+1, 0.001)
      tmp_y   <- upperbounds( t = tmp_x, tU = tmp_max, pU = tmp_pU )
      tmp_labcab <- substr( x = calibs[i,1], start = 1, stop = 12 )
    }
    # 2. Plot calibration density VS marginal density used by MCMCtree
    
    # 2.1. Find limit axis
    max_x_chain <- round( max( density( divt_list[[ 1 ]][,tmp_node] )$x ) + 0.5 )
    min_x_chain <- round( min( density( divt_list[[ 1 ]][,tmp_node] )$x ) - 0.5 )
    x_lim       <- c( min( min_x_chain, min( tmp_x ) ),
                      max( max_x_chain, max( tmp_x ) ) )
    max_y_chain <- round( max( density( divt_list[[ 1 ]][,tmp_node] )$y ) + 0.5 )
    min_y_chain <- round( min( density( divt_list[[ 1 ]][,tmp_node] )$y ) - 0.5 )
    y_lim       <- c( min( min_y_chain, min( tmp_y ) ),
                      max( max_y_chain, max( tmp_y ) ) )
    
    # 2.2. Plot
    plot( density( divt_list[[ 1 ]][,tmp_node], adj = 1 ),
          xlim = c( x_lim[1], x_lim[2] ), ylim = c( y_lim[1], y_lim[2] ),
          main=NULL, xlab = '', ylab = '', cex.axis = 0.7, mgp = c(2.5,0,0),
          #xaxt = "n", yaxt = "n", 
          tck = 0.01, col.main = "white" )
    if( length( is_B ) == 1 | length( is_U ) == 1 | length( is_L ) == 1 ){
      lines( x = tmp_x, y = tmp_y, adj = 1 , col = "blue" )
    }
    if( length( is_ST ) == 1 ){
      curve( dst( x, xi = tmp_ST_pars[1], omega = tmp_ST_pars[2],
                  alpha = tmp_ST_pars[3], nu = tmp_ST_pars[4] ),
             from = 0, to = x_lim[2],
             n = 1e4, add = TRUE, col = "blue" )
    }
    if( length( is_SN ) == 1 ){
      curve( dst( x, xi = tmp_SN_pars[1], omega = tmp_SN_pars[2],
                  alpha = tmp_SN_pars[3] ),
             from = 0, to = x_lim[2],
             n = 1e4, add = TRUE, col = "blue" )
    }
    title( main = colnames( divt_list[[ 1 ]] )[tmp_node], 
           line = -2.7, sub = tmp_labcab, cex.main = 1.3, cex.sub = 0.7, adj = 0.9 )
    if ( ind == TRUE ){
      if( ! dir.exists( paste( main_wd, "plots/calVSmarg/", dat, "/ind", sep = "" ) ) ){
        dir.create( paste( main_wd, "plots/calVSmarg/", dat, "/ind", sep = "" ) )
      }
      png( filename = paste( main_wd, "plots/calVSmarg/",  dat, "/ind/CheckUSPvsEP_tn", 
                             calibs[i,2], "_", calibs[i,1], ".png", sep = "" ),
           width = 1024, height = 768 )
      
      plot( density( divt_list[[ 1 ]][,tmp_node], adj = 1 ),
            xlim = c( x_lim[1], x_lim[2] ), ylim = c( y_lim[1], y_lim[2] ),
            main=NULL, xlab = '', ylab = '', cex.axis = 0.7, mgp = c(2.5,0,0),
            #xaxt = "n", yaxt = "n", 
            tck = 0.01, col.main = "white" )
      if( length( is_B ) == 1 | length( is_U ) == 1 | length( is_L ) == 1 ){
        lines( x = tmp_x, y = tmp_y, adj = 1 , col = "blue" )
      }
      if( length( is_ST ) == 1 ){
        curve( dst( x, xi = tmp_ST_pars[1], omega = tmp_ST_pars[2],
                    alpha = tmp_ST_pars[3], nu = tmp_ST_pars[4] ),
               from = 0, to = x_lim[2],
               n = 1e4, add = TRUE, col = "blue" )
      }
      if( length( is_SN ) == 1 ){
        curve( dst( x, xi = tmp_SN_pars[1], omega = tmp_SN_pars[2],
                    alpha = tmp_SN_pars[3] ),
               from = 0, to = x_lim[2],
               n = 1e4, add = TRUE, col = "blue" )
      }
      title( main = colnames( divt_list[[ 1 ]] )[tmp_node], 
             line = -20, sub = tmp_labcab, cex.main = 3, cex.sub = 2, adj = 0.9 )
      info.legend <- c( "Calibration density", "Marginal density" )
      col.legend  <- c( "blue", "black" )
      legend( "topright", legend = info.legend, col = col.legend,
              lty = 1, bty = 'n', cex = 2 )
      dev.off()
    }
  }
  
  plot( x = 1, y = 1, col = "white", xaxt = "n", yaxt = "n", xlab = '', ylab = '' )
  info.legend <- c( "Calibration density",
                    "Marginal density" )
  col.legend  <- c( "blue", "black" )
  #coords.plot <- locator()
  legend( "top", legend = info.legend, col = col.legend,
          lty = 1, bty = 'n', cex = 1.3 )
  # Switch off main graphics device
  dev.off()
  
}

# Function to fit ST distributions
# 
# Arguments
#
# divtimes      List, object generated with `load_data` function.
# out_log       Boolean, TRUE if you want a log file to be output, FALSE otherwise.
# out_plots     Boolean, TRUE if you want plots of the fitted ST distributions to
#               the posterior distributions.
# out_path      Character, path to the main directory where output files will be
#               written out if enabled. It needs to have a "/" at the end of the
#               path.
# name_out      Character, name of the directory where you want to output your 
#               files,
# ST_empty_obj  List, empty list created before running the pipeline.
# ST_empty_dist List, empty list created before running the pipeline.
#
fitST <- function( divtimes, out_log = FALSE, out_plots = FALSE,
                   out_path = FALSE, name_out = FALSE,
                   ST_empty_obj, ST_empty_dist)
{
  
  # 1. Create directories for output files if required by user
  if( out_path != FALSE ){
    
    if( ! dir.exists( paste( out_path, "logs", sep = "" ) ) && out_log == TRUE ){
      dir.create( paste( out_path, "logs/", sep = "" ) )
    }
    if( ! dir.exists( paste( out_path, "logs/", name_out, sep = "" ) ) && out_log == TRUE  ){
      dir.create( paste( out_path, "logs/", name_out, sep = "" ) )
    }
    if( ! dir.exists( paste( out_path, "plots", sep = "" ) ) && out_plots == TRUE ){
      dir.create( paste( out_path, "plots/", sep = "" ) )
      
    }
    if( ! dir.exists( paste( out_path, "plots/", name_out, sep = "" ) ) && out_log == TRUE  ){
      dir.create( paste( out_path, "plots/", name_out, sep = "" ) )
    }
  }
  
  # 2. Fit a ST distribution to each node. Then save in the lists previously
  # created both the object output by sn::st.mple and only the "dp" pars
  for ( i in 1:length( colnames( divtimes ) ) ){
    
    # 2.1. Start working with node `colnames( divtimes )[i]`
    cat( "Working with node", colnames( divtimes )[i], "...\n\n" )
    if( out_log == TRUE ){
      write( paste( "Working with node ", colnames( divtimes )[i], "...\n\n", sep = "" ),
             file = paste( out_path, "logs/", name_out,
                           "/log_file_convergence_BFGS.txt", sep = "" ),
             sep = "\n", append = TRUE )
    }
    tmp_node <- sn::st.mple( y = divtimes[,i], opt.method = "BFGS" )
    
    # 2.2. Check for convergence, otherwise keep trying
    count_tries_conv <- 1
    while( tmp_node$opt.method$convergence != 0 ){
      
      count_tries_conv <- count_tries_conv + 1
      cat( "Convergence has not been reached with node", colnames( divtimes )[i],
           "...\nSEARCH NUMBER", count_tries_conv, "...\n",
           "The parameters found in the previous search:\n",
           tmp_node$dp[1], "|", tmp_node$dp[2], "|",tmp_node$dp[3], "|",tmp_node$dp[4],
           "\nare now used\as starting values now\n\n" )
      if( out_log == TRUE ){
        write( paste( "Convergence has not been reached with node ", colnames( divtimes )[i],
                      "...\nSEARCH NUMBER ", count_tries_conv, "...\n",
                      "The parameters found in the previous search are used\n",
                      "as starting values now:\n",
                      tmp_node$dp[1], "|", tmp_node$dp[2], "|",tmp_node$dp[3], "|",tmp_node$dp[4],
                      "\n\n", sep = "" ),
               file = paste( out_path, "logs/", name_out,
                             "/log_file_convergence_BFGS.txt", sep = "" ), 
               sep = "\n", append = TRUE )
      }
      
      tmp_node <- sn::st.mple( y = divtimes[,i], opt.method = "BFGS",
                               dp = tmp_node$dp.complete )
      
      if( tmp_node$opt.method$convergence == 0 ){
        cat( "Convergenced reached now!\n" )
        cat( "Final parameters for node", colnames( divtimes )[i], "are:\n",
             tmp_node$dp[1], "|", tmp_node$dp[2], "|",tmp_node$dp[3], "|",tmp_node$dp[4],"\n\n" )
        if( out_log == TRUE ){
          write( paste( "Convergenced reached now!\n\n", 
                        "Final parameters for node ", colnames( divtimes )[i], "are:\n",
                        tmp_node$dp[1], "|", tmp_node$dp[2], "|",tmp_node$dp[3], "|",
                        tmp_node$dp[4],"\n\n", sep = "" ),
                 file = paste( out_path, "logs/", name_out,
                               "/log_file_convergence_BFGS.txt", sep = "" ),
                 sep = "\n", append = TRUE )
        }
        
        
      }
      
      if( count_tries_conv == 50 ){
        
        cat( "You have tried 50 times\n",
             "We will try the optimizing approach within this function...\n" )
        if( out_log == TRUE ){
          write( paste( "You have tried 50 times\n",
                        "We will try the optimizing approach within this function...\n", sep = "" ),
                 file = paste( out_path, "logs/", name_out,
                               "/log_file_convergence_BFGS.txt", sep = "" ),
                 sep = "\n", append = TRUE )
        }
        count_tries_conv_FUN <- 0
        
        while( tmp_node$opt.method$convergence != 0 ){
          
          count_tries_conv_FUN <- count_tries_conv_FUN + 1
          tmp_node <- sn::st.mple( y = divtimes[,i], dp = tmp_node$dp.complete )
          if( tmp_node$opt.method$convergence == 0 ){
            cat( "Convergenced reached with their method now!\n" )
            cat( "Final parameters for node", colnames( divtimes )[i], "are :\n",
                 tmp_node$dp[1], "|", tmp_node$dp[2], "|",tmp_node$dp[3], "|",
                 tmp_node$dp[4],"\n\n" )
            if( out_log == TRUE ){
              write( paste( "Convergenced reached with their method now!\n\n", 
                            "Final parameters for node ", colnames( divtimes )[i], "are:\n",
                            tmp_node$dp[1], "|", tmp_node$dp[2], "|",tmp_node$dp[3], "|",
                            tmp_node$dp[4],"\n\n", sep = "" ),
                     file = paste( out_path, "logs/", name_out,
                                   "/log_file_convergence_BFGS.txt", sep = "" ),
                     sep = "\n", append = TRUE )
            }
            
          }
          if( count_tries_conv_FUN == 50 ){
            cat( "You have tried 50 times with their approach,",
                 "this is going to be killed\n" )
            if( out_log == TRUE ){
              write( paste( "You have tried 50 times with their approach,",
                            "this is going to be killed\n", sep = "" ),
                     file = paste( out_path, "logs/", name_out,
                                   "/log_file_convergence_BFGS.txt", sep = "" ),
                     sep = "\n", append = TRUE )
            }
            break
          }
          
        }
        
      }
      
    }
    
    # 2.3. Get data 
    ST_empty_obj[[ i ]] <- tmp_node
    ST_empty_dist[[ i ]]   <- tmp_node$dp
    
    # 2.4. (OPTIONAL, DEPENDS ON USER'S DECISION) 
    #    Plot fitted ST for each node previously 
    #    computed using the values sampled during the MCMC
    if( out_plots == TRUE ){
      png( filename = paste( out_path, "plots/", name_out,
                             "/Fit_ST_", colnames( divtimes )[i], ".png", sep = "" ),
           width = 1024, height = 768 )
      # 2.4.1. Find limit axis
      max_x_st    <- round( max( density( rst( n = 1000, dp = tmp_node$dp.complete ) )$x ) + 0.5 )
      max_x_chain <- round( max( density( divtimes[,i] )$x ) + 0.5 )
      x_lim       <- max( max_x_chain, max_x_st )
      max_y_st    <- round( max( density( rst( n = 1000, dp = tmp_node$dp.complete ) )$y ) )
      max_y_chain <- round( max( density( divtimes[,i] )$y ) + 0.5 )
      y_lim       <- max( max_y_chain, max_y_st )
      if( out_log == TRUE ){
        write( paste( colnames( divtimes )[i], max_x_chain, max_y_chain, sep = "\t" ),
               file = paste( out_path, "logs/", name_out,
                             "/log_limaxis.txt", sep = "" ),
               sep = "\n", append = TRUE )
      }
      # 2.4.2. Plot
      plot( density( divtimes[,i], adj = 1 ),
            xlim = c( 0, x_lim ), ylim = c( 0, y_lim ), 
            main = paste( colnames( divtimes )[i], " = ",
                          "ST(", paste0( round(tmp_node$dp, 2),
                                         collapse = "," ),
                          ")", sep = "" ) )
      curve( dst( x, xi = tmp_node$dp[1], omega = tmp_node$dp[2],
                  alpha = tmp_node$dp[3], nu = tmp_node$dp[4] ),
             from = 0, to = x_lim,
             n = 1e4, add = TRUE, col = "red" )
      dev.off()
    }
  }
  
  # 3. Return R objects
  return( list( ST_obj = ST_empty_obj,
                ST_dists = ST_empty_dist ) )
}


# Function that takes as input the ST distribution objects generated with
# function `fitSt()` and outputs a matrix that can be used to specify ST 
# distribtuions as priors to date the subtrees.
# 
# Arguments
#
# ST_dists  Object output by `fitST()` in-house function.
# out_path         Character, path to the main directory where output files will be
#                  written out if enabled. It needs to have a "/" at the end of the
#                  path.
# out_file         Character, name for the output file in tsv format
write_STmat <- function( ST_dists, out_path, out_file )
{
  
  # Create output dir if non-existent
  if( ! dir.exists( paste( out_path, "Rout", sep = "" ) ) ){
    dir.create( paste( out_path, "Rout", sep = "" ) )
  }
  # Create empty matrix to later be filled out
  mat_ST <- matrix( 0, nrow = length( ST_dists ), ncol = 6 )
  #names( ST_dists[[1]] )
  colnames( mat_ST ) <- c( "xi-location", "omega-scale", "alpha-shape", "nu-df",
                           "MCMCtree-calib", "MCMCtree-calib-rounded" )
  rownames( mat_ST ) <- names( ST_dists )
  # Loop over ST distributions in object and extract the four parameter values
  # for the fitted ST distribution to each node, then fill matrix with those
  for ( i in 1:length( ST_dists ) ){
    
    mat_ST[i,1] <- ST_dists[[ i ]][1]
    mat_ST[i,2] <- ST_dists[[ i ]][2]
    mat_ST[i,3] <- ST_dists[[ i ]][3]
    mat_ST[i,4] <- ST_dists[[ i ]][4]
    
    tmp_ST_dists_rounded <- round( ST_dists[[ i ]], 3 )
    mat_ST[i,5] <- paste( "ST(", ST_dists[[ i ]][1], ",", ST_dists[[ i ]][2],
                          ",", ST_dists[[ i ]][3], ",", ST_dists[[ i ]][4], ")",
                          sep = "" )
    mat_ST[i,6] <- paste( "ST(", tmp_ST_dists_rounded[1], ",",
                          tmp_ST_dists_rounded[2], ",", tmp_ST_dists_rounded[3],
                          ",", tmp_ST_dists_rounded[4], ")", sep = "" )
    
  }
  # Output matrix in tsv format
  write.table( mat_ST, file = paste( out_path, "Rout/", out_file, ".tsv",
                                     sep = "" ),
               sep = "\t", quote = F )
}

# Function that integrates all the previous functions to run a QC check
#
# Parameters
# num_dirs       Numeric. Default is 36 as the number of chains run in a first
#                instance.
# delcol         Numeric, number of columns that are to be deleted as they do
#                not contain divtimes. You can read the header of the `mcmc.txt`
#                files and count the number of `mu*` and `sigma2*` elements. Do
#                not count the `lnL` value as this is automatically ignored in
#                the function `load_dat`. Assuming an MCMC run under a
#                relaxed-clock model with no partitions, we would see `mu` and
#                `sigma2` columns. Therefore, `delcol = 2`. The default value is
#                2, but you will need to change it if you run `MCMCtree` with
#                different settings.
# path           Character. Path to the directory where the analyses for this
#                dataset ran and the `mcmc.txt` file can be found.
# num_divt       Numeric. Number of columns of the `mcmc.txt` that correspond
#                to the samples collected for the times.
# node_calib     Character. CSV file with two columns, where the first column
#                has the names of the calibrations and in the second the
#                corresponding node.
# dataset        Character, name of the dataset that will be used to  
#                generate the name of the output directory.
# perc           Numeric. Percentile to calculate the quantiles. Default: 0.975.
# def_samples    Numeric. Number of samples that the user defined through the
#                `MCMCtree` option `nsample. 
# prior          Boolean. TRUE if `MCMCtree` has been sampling from the prior,
#                FALSE otherwise.
# out_dat        Character, path to the directory to save output data. Do not
#                add the last "/".
# time_unit      Integer, value that you need to multiply the estimated divtimes
#                with
# out_file_pdf   Character, name to be given to the convergence plot that will
#                be generated. If filtering is required, "_FILT" will be 
#                appended to the name given by the user.
# out_title_pdf  Character, title to be used for the convergence plot that will
#                be generated. If filtering is required, "(FILT)" will be 
#                appended to the title given by the user.
# th             Numeric, threshold over which chains will be evaluated as
#                "problematic".
# outchecks_dir  Character, path to the `plots` directory created at the
#                beginning of the script. Same object name has been used to 
#                avoid problems
QC_conv <- function( num_dirs, delcol, path, num_divt, node_calib, dataset, perc, 
                     def_samples, prior = FALSE, out_dat, time_unit, 
                     out_file_pdf, out_title_pdf, 
                     th, outchecks_dir ){
  # 1. Obtain summarise object and generate first convergence plot
  cat( "\n[[ DATASET TO ANALYSE: ", dataset, " ]]\n[[ PATH:", path, " ]]\n" )
  cat( "\n    * * *    \n\n")
  obj_sum <- sum_MCMC( num_dirs = num_dirs, delcol = delcol,
                       data_dir = path,
                       num_divt = num_divt, node_calib = node_calib,
                       dataset = dataset,
                       perc = perc, def_samples = def_samples,
                       prior = prior,
                       # In this case, I want the output dir to be
                       # created within the same data dir, but it
                       # could be changed!
                       out_dat = path,
                       time_unit = time_unit # time unit in Ma
  )
  cat( "\n[[ RUNNING QC CHECKS ]]\n" )
  if( length( num_dirs ) == 1 ){
    num_chains   <- num_dirs
  }else{
    num_chains   <- length( num_dirs )
  }
  cat( ">> There are a total of", num_chains, "chains\n" )
  half_chains <- trunc( num_chains/2 )
  if( half_chains == 1 ){
    post_half1 <- obj_sum$mean[1,]
  }else{
    post_half1  <- apply( X = obj_sum$mean[1:half_chains,],
                          MARGIN = 2, FUN = mean )
  }
  if( c(half_chains+1) == num_chains ){
    post_half2  <- obj_sum$mean[num_chain,]
  }else{
    post_half2  <- apply( X = obj_sum$mean[(half_chains+1):num_chains,],
                          MARGIN = 2, FUN = mean )
  }
  # Create directory in case users have forgotten to run the script
  # from the beginning or deleted the `plots` directory accidentally
  if( ! dir.exists( outchecks_dir ) ){
    dir.create( outchecks_dir )
  }
  if( ! dir.exists( paste( outchecks_dir, "ESS_and_chains_convergence/",
                           sep = "" ) ) ){
    dir.create( paste( outchecks_dir, "ESS_and_chains_convergence/",
                       sep = "" ) )
  }
  pdf( paste( outchecks_dir,
              "ESS_and_chains_convergence/", out_file_pdf, ".pdf",
              sep = "" ),
       paper = "a4" )
  plot_convergence( name_dir_dat = out_title_pdf,
                    mean_divt1 = post_half1,
                    mean_divt2 = post_half2, num_runs = num_chains )
  dev.off()
  # Plot traces -- uncomment below and run only if you have enough space!
  # Otherwise, you can always quickly run Tracer or another graphical software
  # to inspect the final `mcmc.txt` files that we will generate with the chains
  # that pass the filters (see README.md file).
  #
  # plot_traces( name_dir_dat = "flex_conc_GBM",
  #              sum_dat = obj_sum,
  #              divt = c(start_divt:stop_divt), n_chains = num_chains,
  #              out_dir = path )
  
  #> CHECK: Looking for q97.5% and q2.5% issues across the chains. Basically, this 
  #> function compares the divergence times estimated for each node across all 
  #> chains. The difference between time estimates is then computed. Differences
  #> larger than the threshold set by the user will be flagged, and hence the 
  #> corresponding chain will be flagged as problematic.
  #> 
  #> If you are working with deep phylogenies, you may want to be more generous 
  #> with the threshold (e.g., >0.4). Otherwise, a threshold lower than 0.4 can 
  #> be too stringent and many chains will be flagged as "problematic". Bear in
  #> mind that this threshold should not bee too vague either (e.g., you would
  #> accept a different of +-threshold in time estimates at the 97.5% and 2.5%
  #> quantile) as the larger the threshold, the larger the differences between
  #> the estimates collected across chains for the same nodes.
  #> For shallow datasets, you may try a threshold `th = 0.2` so you can better
  #> refine which chains you keep. Once you run the ESS checks, you can make sure
  #> that Rhat is not larger than 1.05 too.
  #> If your analysis has returned flagged chains, please check the output file
  #> `check_chains.txt` that will be generated before deciding whether a chain is
  #> to be kept or deleted from your analysis.
  #> 
  #> 1. Set threshold and run preliminary comparison anlaysis across MCMC runs
  cat( ">> Threshold value to run QC on quantiles: th =", th, "\n\n" )
  th <- th
  sum_quantiles <- check_quantiles( dat = obj_sum, threshold = th,
                                    num_chains = num_chains, compare_all = TRUE )
  #> 2. Find out which is the one that should be used as "main chain" against
  #> which the rest should be compared (i.e., the one with less differences in
  #> q97.5% and q2.5%)
  chain_ind <- which( sum_quantiles$sum_chains[,1] %in% min( sum_quantiles$sum_chains[,1] ) )
  main_ch   <- sum_quantiles$sum_chains[chain_ind,2]
  if( length( main_ch ) > 1 ){
    # Use by default one of the chains that ran first
    # i.e., starting from chain labelled "1"
    main_ch <- min( main_ch )
  }
  #> 3. Find out which chains should be removed by using this main chain
  check_quantiles( dat = obj_sum, threshold = th, main_chain = main_ch,
                   num_chains = num_chains, compare_all = FALSE,
                   outdir = path )
  #>
  #> 2. Chech if there are chains to be deleted
  if( file.exists( paste( path, "check_chains.txt", sep = "" ) ) ){
    cat( "\n>> WARNING: some chains have not passed QC checks\n")
    # Get filtered summary object and generate convergence plots with filtered
    # dataset
    chains_f  <- readLines( paste( path, "check_chains.txt", sep = "" ) )
    ind_lines <- grep( pattern = "in chain", x = chains_f )
    rm_chains <- as.numeric( unique( gsub( pattern = "..*in chain ",
                                           replacement = "" ,
                                           x = chains_f[ind_lines] ) ) )
    filt_chains <- setdiff( c( 1:num_chains ), rm_chains )
    cat( ">> Chains kept: ", filt_chains,
         "\n>> Chains removed: ", rm_chains, "\n" )
    cat( "\n>> Summarising only filtered chains now!\n\n" )
    write.table( x = t( filt_chains ),
                 file = paste( path, "chains_kept.txt", sep = "" ),
                 quote = FALSE, row.names = FALSE, col.names = FALSE )
    # Abort if only one chain has passed the filters
    if( length( filt_chains ) == 1 ){
      stop( "You need to increase the threshold value for this\n",
            "filtering step to continue. You may have very divergent\n",
            "sequences, too few age constraints, or other technical\n",
            "problems may have occurred during the MCMC\n\n",
            "Current threshold value: ", th, "\n" )
    }
    post_FILT_sum <- sum_MCMC( num_dirs = filt_chains, delcol = delcol,
                               data_dir = path,
                               num_divt = num_divt, node_calib = node_calib,
                               dataset = paste( dataset, "_FILT", sep = "" ),
                               perc = perc, def_samples = def_samples,
                               prior = prior,
                               # In this case, I want the output dir to be
                               # created within the same data dir, but it
                               # could be changed!
                               out_dat = path,
                               time_unit = time_unit # time unit
    )
    # Convergence plot
    half_chains <- trunc( length(filt_chains)/2 )
    if( half_chains == 0 ){
      # Double check
      stop( "You need to increase the threshold value for this\n",
            "filtering step to continue. You may have very divergent\n",
            "sequences, too few age constraints, or other technical\n",
            "problems may have occurred during the MCMC\n\n",
            "Current threshold value: ", th, "\n" )
    }else if( half_chains == 1 ){
      post_half1  <- post_FILT_sum$mean[1,]
    }else{
      post_half1  <- apply( X = post_FILT_sum$mean[1:half_chains,],
                            MARGIN = 2, FUN = mean )
    }
    if( c(half_chains+1) == length(filt_chains) ){
      post_half2  <- post_FILT_sum$mean[length(filt_chains),]
    }else{
      post_half2  <- apply( X = post_FILT_sum$mean[(half_chains+1):length(filt_chains),],
                            MARGIN = 2, FUN = mean )
    }
    pdf( paste( outchecks_dir,
                "ESS_and_chains_convergence/", out_file_pdf, "_filt.pdf",
                sep = "" ),
         paper = "a4" )
    plot_convergence( name_dir_dat = paste( out_title_pdf, " (FILT)", sep = "" ),
                      mean_divt1 = post_half1,
                      mean_divt2 = post_half2,
                      num_runs = length(filt_chains) )
    dev.off()
    #> Delete old prior object with unfiltered chains and keep the new one as main
    rm( obj_sum )
    obj_sum <- post_FILT_sum
    rm( post_FILT_sum )
    # Update `num_chain` object so that now only includes num chains filtered
    num_chains <- length(filt_chains)
  }
  
  # 3. Compute ESS with RStan
  ##> Each column is assumed to be an MCMC. Rows are iterations for parameter X
  ##> Source explaining why it is preferable than the function in coda:
  ##> https://nature.berkeley.edu/~pdevalpine/MCMC_comparisons/nimble_MCMC_comparisons.html
  ##> 
  ##> We will compute the ESS taking into account all the final filtered chains
  ##> NOTE: The same number of rows are required to compute the tail-ESS and 
  ##> bulk-ESS. The second argument of the in-house function `sum_MCMC_ESS` used
  ##> below is a vector with the number of samples collected for each independent
  ##> chain. Once the chain with less samples collected has been identified, then
  ##> we can crop the number of samples in each element of the array to fit 
  ##> that number. This is essentially what the in-house function `sum_MCMC_ESS`
  ##> does below. In that way, only the minimum number of samples collected
  ##> across all independent chains are used for all the chains, even though 
  ##> more samples have been collected -- more conservative, but only way the 
  ##> function can be used properly to my knowledge.
  cat( "\n[[ CALCULATING ESS STATS ]]\n")
  ESS <- sum_MCMC_ESS( x = obj_sum$arr4stan,
                       samp_per_chain = obj_sum$samp_per_chain )
  # Create a list
  ESS_results <- vector( "list", 8 )
  names( ESS_results ) <- c( "median", "min", "max", "num_samples_for_ESS",
                             "ESS", "minRhat", "maxRhat", "total_samples" )
  ESS_results[[1]] <- median( obj_sum$samp_per_chain  ) # Median samples per chain
  ESS_results[[2]] <- min( obj_sum$samp_per_chain  )  # Minimum samples per chain
  ESS_results[[3]] <- max( obj_sum$samp_per_chain  ) # Maximum samples per chain
  ESS_results[[4]] <- min( obj_sum$samp_per_chain  ) * num_chains # Number of samples used to compute the ESS
  ESS_results[[5]] <- ESS$tab # Show tail-ESS and bulk-ESS
  ESS_results[[6]] <- min( ESS$stats$Rhat ) # Calculate Rhat, min and max. Good if max Rhat <= 1.05
  ESS_results[[7]] <- max( ESS$stats$Rhat ) # Calculate Rhat, min and max. Good if max Rhat <= 1.05
  ESS_results[[8]] <- dim( obj_sum$all_mcmc ) # Number of samples collected throughout
  # Check which nodes have Rhat >= 1.05
  if( min( ESS$stats$Rhat ) > 1.05 | max( ESS$stats$Rhat ) > 1.05 ){
    cat( "\n[[ CHECKING WHICH NODES HAVE NOT CONVERGED ]]\n")
    tmp_ESS  <- ESS$stats
    tmp_ind  <- which( ESS$stats$Rhat > 1.05 )
    tmp_divt <- as.numeric( gsub( pattern = "t_n", replacement = "", 
                                  x = names( ESS$stats[tmp_ind,10] ) ) )
    count <- 0
    for( i in tmp_divt ){
      count <- count + 1
      tmp_df <- which( node_calib[,2] %in% i )
      if( length( tmp_df ) > 0 ){
        rownames(tmp_ESS)[tmp_ind[count]] <- paste( i, 
                                                    "-",
                                                    node_calib[tmp_df,1],
                                                    sep = "" )
      }else{
        cat( "Node", i, "was not calibrated\n" )
      }
    }
    tmp_ESS[tmp_ind,10]
    tmp_dataset <- gsub( pattern = "_", replacement = "", x = dataset )
    write.table( x = tmp_ESS[tmp_ind,10],
                 file = paste( path,
                               "problem_nodes_conv_", tmp_dataset, ".txt",
                               sep = "" ),
                 quote = FALSE, row.names = TRUE, col.names = c( "Rhat" ) )
  }
  
  cat( "\n[[ RETURNING OBJECTS, END OF TASKS ]]\n")
  cat( "\n    * * *    \n\n")
  # Return objects!
  if( min( ESS$stats$Rhat ) > 1.05 | max( ESS$stats$Rhat ) > 1.05 ){
    return( list( obj_sum = obj_sum, ESS_obj = ESS, 
                  ESS_results = ESS_results,
                  not_conv_nodes = tmp_ESS[tmp_ind,10],
                  num_chains = num_chains ) )
  }else{
    return( list( obj_sum = obj_sum, ESS_obj = ESS, 
                  ESS_results = ESS_results, num_chains = num_chains ) )
  }
  
}

# Function to incorporate node age constraints following PAML notation
#
# Arguments:
#
# tt            Phylo, object with the tree, previously generated
# calibrations  Matrix, element of the list generated with the calibration
#               info. E.g., `all_calibs[[1]]`
# out_name      Character, name of the dataset analysed
# out_dir_raw   Character, abs/rel path to the output directory for raw
#               trees
# out_dir_inp   Character, abs/rel path to th eoutput directory for inp data
#
add_node_const <- function( tt, calibrations, out_name, out_dir_raw, out_dir_inp )
{
  # Print message regarding the data
  cat( "\n\n* * *\nCalibrated tree to be generated: ", out_name, "\n* * *\n" )
  
  #-------------------------------#
  # PART 1: GENERAT TEMPLATE FILE #
  #-------------------------------#
  # Generate empty vector
  keep_indexes <- matrix( 0, nrow = length(rownames(calibrations)), ncol = 3 )
  # Generate empty vector with as many entries as nodes in the tree
  tt$node.label <- rep( NA, tt$Nnode )
  for( i in 1:length(rownames(calibrations)) ){
    ## Build MCMCtree calib
    ind_nas  <- which( calibrations[i,] == "" )
    if( c(8,9,10) %in% ind_nas && length( ind_nas ) == 1 ){ # if B calib
      node_lab <- paste( "'B(", calibrations[i,4], ",",
                         calibrations[i,6], ",", calibrations[i,5],
                         ",", calibrations[i,7], ")'", sep = "" )
    }else if( length( ind_nas ) == 3 ){ # col8 + two cols for U or L
      # if lower bounds missing, then it is an upper bound calib
      if( 4 %in% ind_nas || 5 %in% ind_nas ){
        node_lab <- paste( "'U(", calibrations[i,6], ",",
                           calibrations[i,7], ")'", sep = "" )
      } # if upper bounds missing, then it is a lower bound calib
      else if( 6 %in% ind_nas || 7 %in% ind_nas ){
        node_lab <- paste( "'L(", calibrations[i,4], ",",
                           calibrations[i,5], ")'", sep = "" )
      }
      else{
        node_lab <- ""
      }
    }else{
      # For inequality calibrations that start with '#[0-9]' or for other MCMCtree cals
      # Only column 8 is allowed
      node_lab <- calibrations[i,8]
    }
    
    # Get MRCA for these two tips
    mrca <- ape::getMRCA( phy = tt, tip = c( calibrations[i,2],
                                             calibrations[i,3]) )
    keep_indexes[i,1] <- mrca-ape::Ntip(tt)
    keep_indexes[i,2] <- calibrations[i,1]
    keep_indexes[i,3] <- paste( calibrations[i,2], "-", calibrations[i,3], "-",
                                node_lab, sep = "" )
    
    cat( "Node label with calibration: ", mrca-ape::Ntip(tt), "\n" )
    # Replace node label accordingly
    tt$node.label[mrca-ape::Ntip(tt)] <- paste0( "[", calibrations[i,1],
                                                 "]", collapse = "" )
    ## While here, I will not use `tt_cals` as the format is wrong
    ## once the tree is written out in the next lines outside the `loop`
    #tt_cals$node.label[mrca-ape::Ntip(tt)] <- node_lab
  }
  ## Find duplicates
  ind_dup   <- which( duplicated(keep_indexes[,1]) == TRUE )
  nodes_dup <- which( keep_indexes[,1] %in% as.numeric( keep_indexes[ind_dup,1] ) )
  #keep_indexes[nodes_dup,]
  if( length( nodes_dup ) > 0 ){
    stop( "There are duplicated nodes, you need to check!" )
  }
  # Remove "NA" from the labs
  ind_na_bools <- is.na( x = tt$node.label )
  ind_na       <- which( ind_na_bools == TRUE )
  tt$node.label[ind_na] <- ""
  # Write PHYLIP header, then the calibrated tree
  writeLines( text = paste( length(tt$tip.label), " 1", sep = "" ), 
              con = paste( out_dir_raw, "cals_only_", out_name, ".tree",
                           sep = "" ) )
  ape::write.tree( phy = tt,
                   file = paste( out_dir_raw, "cals_only_", out_name, ".tree",
                                 sep = "" ),
                   append = TRUE )
  
  #-----------------------------------------#
  # PART 2: READ TREE AND CALIBRATIONS FILE #
  #-----------------------------------------#
  # Read tree and get phylip header
  # NOTE: Make always sure that there is at least one blank line at the 
  # end of the tree file! Otherwise, you will get an error telling you that 
  # there is an incomplete final line in these files.
  tt_name       <- paste( out_dir_raw, "cals_only_", out_name, ".tree",
                          sep = "" )
  tt            <- readLines( tt_name )
  phylip.header <- tt[1]
  tt            <- tt2 <- tt3 <- tt[2]
  
  #----------------------------------------#
  # PART 3: REPLACE TAGS WITH CALIBRATIONS #
  #----------------------------------------#
  # Replace calibration names with corresponding calibration
  for( j in 1:length(rownames(calibrations)) ){
    # Build MCMCtree calib that will be replaced
    ind_nas  <- which( calibrations[j,4:7] == "" )
    if( length( ind_nas ) == 0 ){ # For soft bounds when cols 4-7 are given...
      node_lab <- paste( "B(", calibrations[j,4], ",",
                         calibrations[j,6], ",", calibrations[j,5], 
                         ",", calibrations[j,7], ")", sep = "" )
    }else if ( length( ind_nas ) == 4 ){ # If only MCMCtree is given in col8...
      node_lab <- calibrations[j,8]
    }else{ # For lower/upper bounds when cols 4-5 or 6-7 are given...
      # NOTES: 1 = minage | 2 = mintail | 3 = maxage | 4 = maxtal
      # If upper bound because NAs are in 1 and 2...
      if( 1 %in% ind_nas || 2 %in% ind_nas ){
        node_lab <- paste( "U(", calibrations[j,6], ",",
                           calibrations[j,7], ")", sep = "" )
      } # If lower bound because NAs are in 3 and 4...
      else if( 3 %in% ind_nas || 4 %in% ind_nas ){
        node_lab <- paste( "L(", calibrations[j,4], ",", 
                           calibrations[j,5], ")", sep = "" )
      }
      else{
        node_lab <- ""
      }
    }
    # Now, `nodel_lab` will have the calibration in `MCMCtree` format already
    # ready to replace the corresponding flag generated in the previous step
    # for each calibrated node
    #
    # Conditional is used so that the single quotation marks are only kept 
    # in the upper-bound calibration for the root. Inequality calibrations
    # do not require single quotation marks
    tmp_calib <- gsub( x = node_lab, pattern = "\\(..*",
                       replacement = "" )
    tmp_calib <- gsub( x = tmp_calib, pattern = "[0-9]..*",
                       replacement = "" )
    if( tmp_calib == 'B' || tmp_calib == 'U' || tmp_calib == 'L' || tmp_calib == 'ST' ){
      tt <- gsub( pattern = paste0("\\[",calibrations[j,1],"\\]"),
                  x = tt,
                  replacement = paste( "'", node_lab, "'", sep = "" ) )
    }else{ # For cross-braced nodes
      tt <- gsub( pattern = paste0("\\[",calibrations[j,1],"\\]"),
                  x = tt,
                  replacement = paste( node_lab, sep = "" ) )
    }
    # Copy to visualise in FigTree
    reps <- gsub( x = gsub( x = gsub( x = gsub( x = gsub( x = node_lab,
                                                          pattern = "\\{",
                                                          replacement = "(" ),
                                                pattern = "\\}",
                                                replacement = ")" ), 
                                      pattern = "\\[|\\]", replacement = "" ),
                            pattern = "\\#", replacement = "flag" ),
                  pattern = " ", replacement = "-" )
    # For cross-braced calibrations without fossil
    if( tmp_calib == '#' ){
      reps <- gsub( x = gsub( x = reps, pattern = "\\#", replacement = "flag" ),
                    pattern = "\\]", replacement = "" )
      tt2 <- gsub( pattern = paste0("\\[",calibrations[j,1],"\\]"),
                   x = tt2,
                   replacement = paste0( "'", reps, "-", calibrations[j,1], "'", 
                                         collapse = "" ) )
    }else{ # For the rest of calibrations
      tt2 <- gsub( pattern = paste0("\\[",calibrations[j,1],"\\]"),
                   x = tt2,
                   replacement = paste0( "'", reps, "-", calibrations[j,1], "'",
                                         collapse = "" ) )
    }
    # Generate an uncalibrated tree for CODEML!
    tt3 <- gsub( pattern = paste0("\\[",calibrations[j,1],"\\]"),
                 x = tt3,
                 replacement = "" )
  }
  
  #---------------------------------------#
  # PART 4: WRITE CALIBRATED TREE IN FILE #
  #---------------------------------------#
  out_dir <- out_dir_inp
  if( ! dir.exists( out_dir_inp ) ){
    dir.create( out_dir_inp )
  }
  write( x = phylip.header, file = paste( out_dir, out_name, "_calib_MCMCtree.tree", sep = "" ) )
  write( x = tt, file = paste( out_dir, out_name, "_calib_MCMCtree.tree", sep = "" ),
         append = TRUE )
  write( x = phylip.header, file = paste( out_dir_raw, out_name,
                                          "_fordisplay_calib_MCMCtree.tree",
                                          sep = "" ) )
  write( x = tt2, file = paste( out_dir_raw, out_name,
                                "_fordisplay_calib_MCMCtree.tree", sep = "" ),
         append = TRUE )
  # write( x = phylip.header, file = paste( out_dir, out_name, "_uncalib.tree", sep = "" ) )
  # write( x = tt3, file = paste( out_dir, out_name, "_uncalib.tree", sep = "" ),
  #        append = TRUE )
  
  # Return objects
  return( list( tt_calib=tt, keep_indexes=keep_indexes ) )
}


