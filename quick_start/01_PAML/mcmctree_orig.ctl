          seed = -1                * Used timestamp to define seed number
       seqfile = mtCDNApri123.txt  * Path to input sequence file
      treefile = mtCDNApri.trees   * Path to input tree file
      mcmcfile = mcmc.txt          * Path to output file with MCMC samples
       outfile = out.txt           * Path to log output file

         ndata = 3        * Number of partitions
       seqtype = 0        * 0: nucleotides; 1:codons; 2:AAs
       usedata = 2        * 0: no data (prior); 1:exact likelihood; 2:approximate likelihood; 3:out.BV (in.BV)
         clock = 2        * 1: global clock; 2: independent rates; 3: correlated rates

         model = 4        * 0:JC69, 1:K80, 2:F81, 3:F84, 4:HKY85
         alpha = 0.5      * alpha for gamma rates at sites
         ncatG = 5        * No. categories in discrete gamma

     cleandata = 0        * remove sites with ambiguity data (1:yes, 0:no)?

       BDparas = 1 1 0.1  * birth, death, sampling
   rgene_gamma = 2 20     * gammaDir prior for rate for genes
  sigma2_gamma = 1 10     * gammaDir prior for sigma^2     (for clock=2 or 3)

         print = 1        * 0: no mcmc sample; 1: everything except branch rates 2: everything
        burnin = 2000     * Number of iterations that will be discarded as part of burn-in
      sampfreq = 100      * Sampling frequency
       nsample = 20000    * Number of samples to be collected

*** Note: Make your window wider (100 columns) before running the program.