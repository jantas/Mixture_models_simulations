# Fits a mixture model of K Markov chains using EM algorithm
# The formulas for the maximization step were derived
# analytically in "Model-based biclustering of clickstream data" by Volodymyr Melnykov

rm(list=ls())
# library(beepr)

##### Functions ########################
source("EM_MC.R")
source("simulate_MC_data.R")

#####  Parameters for Data Simulation #####
P1 <- matrix(rep(1/4,4*4),ncol = 4)  # uniform
P2 <- matrix(c(0.3,0.5,0.1,0.1,0.5,0.3,0.1,0.1,0.4,0.4,0.1,0.1,0.4,0.4,0.1,0.1),byrow = T, ncol=4)  # hangs out at 1,2
P3 <- matrix(c(0.1,0.1,0.5,0.3, 0.1,0.1,0.3,0.5,0.1,0.1,0.4,0.4,0.1,0.1,0.4,0.4),byrow = T, ncol=4)  # hangs out at 3,4
a1 <- a2 <- a3 <- matrix(c(1/4,1/4,1/4,1/4),ncol=1)
Gamma.init <- list(P1,P2,P3)
alpha.init <- list(a1,a2,a3)
J <- 4


##### Change this ###########
S <- 100  # the number of simulations
NI <- 50  # number of random initializations

X.T <- 20  # the length of MCs (10, 20, 50, 100)
N <- 600  # the nimber of MCs (600, 1500)
K_start <- 1

EM.CEM <- 1  # 0... standard EM, 1... classification EM (CEM)
EM.CEM.stoch <- 0  # 0... classification EM (CEM), 1... stochastic EM (SEM)

# stopping criterions
tol <- 1e-5  # 1e-7 for EM is fine but for CEM and SEM 1e-5 is better
iter.max <- 100  #  the actual criterion is K*iter.max so it considers the incresing complexity with increasing K 

##### Global Variables
y.levels <<- CombSet(1:J, m=2,repl=TRUE, ord=TRUE)  # y.levels are different combinations of the MC states for finding trans. prob.
y.levels <<- paste(y.levels[,1],y.levels[,2],sep="")
y.levels <<- as.array(as.numeric(y.levels))
for (K in K_start:5){
  
  results.model <- vector(mode = "list", length=S)
  results.fit_indic <- vector(mode = "numeric", length=S)
  
  for(s in 1:S){  # sth simulation 
    set.seed(s)
    
    print(sprintf('K = %d/5, s = %d/%d ', K,s,S))
    print('end, nr of iter,    BIC    , fit_indicator')
    
    # sth dataset
    dataset <- simulate_MC_data(Gamma.init,alpha.init,X.T,N)
    X.T <- length(dataset[1,])-1  # the length of the MC
    ID.true <- dataset[,X.T+1]
    X <- dataset[,1:(X.T)]
    
    results.RI.model <- vector(mode = "list", length=NI)  # results for a given dataset but for different initial values
    results.RI.fit_indic <- vector(mode = "numeric", length=NI)
    
    ##### Start the clock!
    ptm <- proc.time()
    
    if(K==1){NI.K <- 1}else{NI.K <- NI}  # when K=1 it is just one iteration, no need for EM
    
    for(i in 1:NI.K){  # ith random initialization
      
      # Random Initial Conditions (Works for ALL orders)
      random_init.output <- random_init(K,J)
      Gamma <- random_init.output[[1]]
      alpha <- random_init.output[[2]]
      tau <- random_init.output[[3]]
      
      param.init <- list(tau, Gamma, alpha)
      
      EM.output <- EM_MC(X, param.init, tol,K*iter.max,EM.CEM, EM.CEM.stoch,ID.plot=ID.true, ID.true)
      
      results.RI.model[[i]] <- EM.output[[1]]
      results.RI.fit_indic[i] <- EM.output[[2]]  # for simulated datasets
    }
    #### Stop the clock
    tm <- proc.time() - ptm
    print(sprintf('%2.00f sec', tm[3]))
    
    # pick the best model (w.r.t. different initialization) (with min BIC)
    models.s.BICs <- unlist(lapply(results.RI.model, "[[",4))
    model.best.ind <- which(models.s.BICs==min(models.s.BICs,na.rm = TRUE))
    print(sprintf('coefficient of variation of BICs for different initializations = %f',sd(models.s.BICs,na.rm = TRUE)/mean(models.s.BICs,na.rm = TRUE)))
    
    # save the best model (w.r.t. different initialization)
    results.model[[s]] <- results.RI.model[[model.best.ind[1]]]
    results.fit_indic[s] <- results.RI.fit_indic[[model.best.ind[1]]]
  }  # end of S cycle
  
  results <- list(results.model,results.fit_indic)
  
  # name a file based on thr parameters used
  if(EM.CEM==1){
    if(EM.CEM.stoch==1){
      filename <- sprintf('MC_simul_SEM_K-%d_N-%d_T-%d.RData',K, N, X.T)
    }else{
      filename <- sprintf('MC_simul_CEM_K-%d_N-%d_T-%d.RData',K, N, X.T)
    }
  }else{
    filename <-    sprintf('MC_simul_EM_K-%d_N-%d_T-%d.RData',K, N, X.T)
  }
  
  save(results, file=filename)
}
