          seed = -1
       seqfile = ALN
      treefile = TREE
      mcmcfile = MCMC
       outfile = out.txt

         ndata = 1
       seqtype = 2    * 0: nucleotides; 1:codons; 2:AAs
       usedata = 0
                      * 2:approximate likelihood; 3:out.BV (in.BV)
         clock = 1

         model = 0
                      *     0:poisson, 1:proportional,2:Empirical,3:Empirical+F
                      *     6:FromCodon, 8:REVaa_0, 9:REVaa(nr=189)
         alpha = 0.5  * alpha for gamma rates at sites
         ncatG = 5    * No. categories in discrete gamma

     cleandata = 0    * remove sites with ambiguity data (1:yes, 0:no)?

       BDparas = 1 1 0.1    * birth, death, sampling

   rgene_gamma = 2 4   * gammaDir prior for rate for genes
  sigma2_gamma = 1 10  * gammaDir prior for sigma^2     (for clock = 1

         print = -1       * 0: no mcmc sample; 1: everything except branch rates 2: everything
        burnin = 100000
      sampfreq = 1000 
       nsample = 20000
