# Fits a mixture model of K Markov chains using EM algorithm
# The formulas for the maximization step were derived 
# analytically in "Model-based biclustering of clickstream data" by Volodymyr Melnykov

rm(list=ls())
library(beepr)

##### Functions ##############################################
source("EM_MC.R")
source("simulate_MC_data.R")
##### SIMULATE DATA ##########################################
################## CHANGE THIS ###############################
P1 <- matrix(rep(1/4,4*4),ncol = 4)  # uniform
P2 <- matrix(c(0.3,0.5,0.1,0.1,0.5,0.3,0.1,0.1,0.4,0.4,0.1,0.1,0.4,0.4,0.1,0.1),byrow = T, ncol=4)  # hangs out at 1,2
P3 <- matrix(c(0,0.1,0.4,0.5, 0.1,0.1,0.4,0.4,0.1,0,0.5,0.4,0.1,0.2,0.6,0.1),byrow = T, ncol=4)  # hangs out at 3,4
a1 <- a2 <- a3 <- matrix(c(1/4,1/4,1/4,1/4),ncol=1)
##############################################################

P <- list(P1,P2,P3)
alpha <- list(a1,a2,a3)

dataset <- simulate_MC_data(P,alpha,X.T=20,N=600)

###### PARAMETERS ############################################
################## CHANGE THIS ###############################
X.T <- length(dataset[1,])-1  # the length of the MC
ID.true <- dataset[,X.T+1]
X <- dataset[,1:(X.T)]
n <- dim(X)[1];
J <- 4  # number of states
K <- 3  # number of components
S <- 3  # number of random initializations
method <- "EM"  # EM, CEM, or SEM
##############################################################

results <- vector(mode = "list", length=S)
for(s in 1:S){
  
  ############# PARAMETERS ##################
  # Random Initial Conditions (Works for ALL orders)
  random_init.output <- random_init(K,J)
  Gamma <- random_init.output[[1]]
  alpha <- random_init.output[[2]]
  tau <- random_init.output[[3]]
  
  param.init <- list(tau, Gamma, alpha)
  
  ##### Global Variables
  y.levels <<- CombSet(1:J, m=2,repl=TRUE, ord=TRUE)  # y.levels are different combinations of the MC states for finding trans. prob.
  y.levels <<- paste(y.levels[,1],y.levels[,2],sep="")
  y.levels <<- as.array(as.numeric(y.levels))
  
  ##### Start the clock!
  ptm <- proc.time()
  if (method == "EM") {EM.CEM=0; EM.CEM.stoch=0} else if (method == "CEM") {EM.CEM=1; EM.CEM.stoch=0} else if (method == "SEM") {EM.CEM=1; EM.CEM.stoch=1};
  EM.output <- EM_MC(X, param.init, tol=1e-6,iter.max=400,EM.CEM, EM.CEM.stoch,ID.plot=NaN, ID.true=NaN)
  model.final <- EM.output[[1]]
  fit_indicator <- EM.output[[2]]  # for simulated datasets
  
  #### Stop the clock
  proc.time() - ptm
  
  results[[s]] <- model.final
  
}

# pick the best model (w.r.t. different initialization) (with min BIC)
models.s.BICs <- unlist(lapply(results, "[[",4))
model.best.ind <- which(models.s.BICs==min(models.s.BICs,na.rm = TRUE))

model.best <- results[[model.best.ind[1]]]
model.best
