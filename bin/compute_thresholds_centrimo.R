# This script will automatically calculate threshold for the motif score and 
# the enrichment zone defined as the distance from the peak summit based on 
# probability density functions and the maximum entropy of the PDF

compute_thresholds = function(data.file, score.type, seq.length.column, distance.column, score.column, window.size, hadb.resid.fold=2, hadb.score.diff=0.20, output.file.name = NULL, print.hits = T){ 

  ### DEBUG
  # setwd("/home/mariughe/Desktop/Projects/ENCODE/thresholding/")
  #data.file = "JAMM_REST_peaks_500bp.fa.score.ext"
  #score.type = "relative"
  #seq.length.column = 7
  #distance.column = 8
  #score.column = 10
  #window.size = 500
  #hadb.resid.fold=2
  #hadb.score.diff=0.20
  #output.file.name = "TEST"
  #print.hits = T
  
  #read data
  file.name = data.file
  data.file = read.delim(file=data.file, as.is=T, header=F )

  #make sure the column indexes are numeric as it may come from command line arguments
  score.column = as.numeric(score.column)
  distance.column = as.numeric(distance.column)
 
 
  #----------------------------------------------------------
  # Calculate centrality p-value
  #----------------------------------------------------------
  adj.pvalue = "NA"
  pvalue = "NA"
  
  pvalues = c() 
  nb.hits = c()
  windows = c()

  seq.length = data.file[1, as.numeric(seq.length.column)]
  window.size = as.numeric(window.size)
  end_loop = window.size - floor(seq.length / 2. - 0.5)
    
#    cat(sprintf("%f, %f, %d, %d\n", as.double(seq.length), as.double(seq.length) / 2. -
#            0.5, floor(seq.length / 2. - 0.5), end_loop))
    for (i in 0:(end_loop-1)){ # One has to stop before end_loop to be able to compute p-value
       # cat(sprintf("%d\n", i))
       enrichment.zone.thresh = i
    
#compute p-value for the enrichment zone (ratio based and not discrete)
       nb.top.hits = length(which((abs(as.numeric(data.file[, distance.column])) <= enrichment.zone.thresh)))
   
       if (nb.top.hits != 0){
            #calculate probability of the entire enrichment zone + peak summit divided by the double window size (both sides of peak summit) + peak summit
            # half of the scored sequence should be removed from upstream and from downstream as it is ignored while sliding the window
            # the "- 0.5" is to compensate for the probability window when dealing with sequences of even length 
            p = ((enrichment.zone.thresh * 2) + 1) / (((window.size * 2) + 1) - floor(seq.length / 2 - 0.5) * 2)
 #           ratio = nb.top.hits / nrow(data.file) * 100 
#            pvalues[i] = log_ibeta(ratio, 100, p)
            #using discrete values to be consistent with JASPAR
            pvalues[i] = log_ibeta(nb.top.hits, nrow(data.file), p)
        }
    }

    #get the minimum p-value corresponding to the enrichment zone threshold
    enrichment.zone.thresh = which(pvalues == min(pvalues))
    pvalue = pvalues[enrichment.zone.thresh]
    #here we do not need window.size * 2 because it is the number of tries, and NOT hits
    adj.pvalue = logev(log(window.size + 1 - (floor(seq.length / 2 - 0.5) * 2)), pvalue)
    write.table(adj.pvalue, paste(file.name, "centrimo", "pval", sep = "."), row.names = F, col.names = F, quote = F)
    
    # calculated p-value is 0 then it is not significant
    #print(sprintf("nb hits in enrichment zone %d is %d", enrichment.threshold, nb.top.hits))
    
    cat('OK')
  # signal too poor to calculate thresholds. Notify the user
   
   #save the results to file for later use
   # plot the hits and the thresholds
   svg(file=paste(file.name, "centrimo","svg", sep = "."), width = 5, height = 5)
#   distance.freq =density(data.file[, distance.column])
#   plot(distance.freq$x,nrow(data.file)*distance.freq$y,type="l",xlab="Distance to peak centre", ylab="Number of motif occurrences", lwd = 2, cex.lab = 1.5)
   
    hist(data.file[,distance.column], breaks = length(table(data.file[,distance.column])), col = "#010116", 
          ylim = c(0,round(max(table(data.file[,distance.column]))*1.1)), xlab = "Distance to peak centre", ylab = "Number of motif occurrences", main = "")
   box(which = "plot")

   bla = dev.off()

 }
  
  #----------------- UTILITIES FOR CENTRIMO P-VALUE CALCULATION -------------------#

  # Compute the log p-value of the extreme value of n independent
  # observations
  # given the single-trial probability p of the observed extreme;
  #  if (n*p < 1e-6) pv = n*p
  #  else if (n*p > A) pv = 1.0
  #  else pv = 1 - (1-A)**(n*p/A)
  #
  logev = function(logn, logp){
    
    A = 1e-6
    LOGMINPROD = -13.8155105579643 # log(1e-6) 
    LOGMAXPROD = 4.60517018598809  # log(100) 
    
    if ((logn)+(logp) < LOGMINPROD){
      return ((logn)+(logp))
    }else if ((logn)+(logp) > LOGMAXPROD){
      return (0)
    }else{
      return (log(1 - (1-A)^exp((logn)+(logp)-LOGMINPROD)))
    }
  }
  
  # Computes the natural logarithm of centrimo p-values based on a binomial distribution.
  
  # Regularized incomplete beta function.
  ricbeta = function(a, b, x) {
    
    # Do not iterate more than 100 times. We should be precise enough.
    MAXIT = 1000;
    EPS = 3.0e-7;
    # Minimum resolution.
    FPMIN = 1.0e-300;
    
    #Initialize constants.
    # double aa, c, d, del, h, qab, qam, qap;
    qab = a + b;
    qap = a + 1.0;
    qam = a - 1.0;
    c = 1.0;
    d = 1.0 - qab*x/qap;
    if (abs(d) < FPMIN) d = FPMIN;
    d = 1.0/d;
    h = d;
    # Everyone likes recursion.
    for (m in 1:MAXIT) {
      m2 = 2*m
      # Even recurrence.
      aa = m*(b - m)*x / ((qam + m2) * (a + m2))
      d = 1.0 + aa * d
      if (abs(d) < FPMIN){
        d = FPMIN
      }
      c = 1.0 + aa/c
      if (abs(c) < FPMIN){
        c = FPMIN
      }
      d = 1.0 / d
      h = h * (d*c)
      # Odd recurrence.
      aa = -(a+m)*(qab + m)*x/((a + m2) * (qap + m2))
      d = 1.0 + aa*d
      if (abs(d) < FPMIN){
        d = FPMIN
      } 
      c = 1.0 + aa/c
      if (abs(c) < FPMIN){
        c = FPMIN
      }
      d = 1.0 / d
      del = d*c
      h = h * del
      if (abs(del - 1.0) < EPS){
        break # Are we done?
      }  
    }
    # This really shouldn't happen at all (please don't happen). 
    if (m > MAXIT){
      stop("a or b too big, or MAXIT too small in ricbeta")
    }
    
    return (h)
    
  }
  
  # Regularized incomplete beta function, log-ed.
  log_ibeta = function(a, b, x) {
    
    #Deal immediately with the number of trials; convert to failures + 1. 
    b = b + 1.0 - a
    
    #Verify probability is valid.
    if ((x < 0.0) || (x > 1.0))  {
      stop("Check probability; must be between 0 and 1.")
    }
    
    #Deal with edge cases. lgamma(x) = log(gamma(x)) = log((x-1)!). Exploit lgamma for precision.
    log_bt = c()
    if ((x == 0.0) || (x == 1.0)) {
      log_bt = 1.0
    } else {
      log_bt = lgamma(a + b) - lgamma(a) - lgamma(b)+ a*log(x)+ b*log(1.0 - x)
    }
    
    # Converge.
    if (x < ((a + 1.0) / (a + b + 2.0))) {
      return (log_bt + log(ricbeta(a, b, x)/a))
    } else {
      return (log(1.0 - exp(log_bt) * ricbeta(b, a, 1.0 - x) / b))
    }
    
  }

  #wrapper for the java function 
  get_threshold = function(int.arr, method){
  rJava::.jinit("/storage/scratch/marius/ENCODE/bin")
  rim <- range(int.arr, na.rm = TRUE)
  im.hist <- as.vector(table(factor(int.arr, levels = rim[1]:rim[2])))
  entropy.class <- rJava::.jnew("MaxEntropy")
  rJava::.jcall(entropy.class, "I", method, im.hist) + rim[1]
 }

  
