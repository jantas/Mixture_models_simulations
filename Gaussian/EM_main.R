# Fits a mixture model of K mutlivariate (p-dimensional) normal distributions N_p(mu, Sigma) using the EM algorithm
# The formulas for the maximization step were derived analytically
# Gaussian mixture pdf: f(x, theta) = sum_{k=1}^K(tau_k * f_k(x; mu_k, Sigma_k)), 
# where K is the number of components (clusters), f_k is the pdf of N_p(mu_k, Sigma_k), theta = (tau, mu, Sigma) is the vector of parameters
# tau = (tau_1, ..., tau_K) are weights (or inclusion probabilities), s.t. 0 < tau_k <= 1, sum(tau) = 1

rm(list=ls())
library(randtoolbox)
library(EMCluster)
library(MixSim)

##### Functions ########################
source("EM_alg_mvnormal_fun.R")
source("Load_dataset.R")

### PARAMETERS THAT CAN BE CHANGED #####
set.seed(100)
K = 3  # nr of components of the model (for synthetic data)
Mu.n <- 2  # nr of random initializations 
gen.params <- list(n = 800, K, p = 2, AvgOmega = 0.01)
dataset_type <- "gen"  # dataset: "iris", "crabs", "gen"

K.fit <- 3;  # nr of clusters the EM algortithm will try to identify
method = "EM";  # EM (standard EM), CEM (classification EM), SEM (stochastic EM) 
K.max <- K.fit  # for initialization of tables for results
##### Load Data ########################
dataset <- load_dataset(dataset = dataset_type, gen.params)  
X <- dataset[[1]]; ID.true <- dataset[[2]];
n <- dim(X)[1]; p <- dim(X)[2]
# plot(X[,1], X[,2], col=ID.true)
pairs(X,col=ID.true)
ID.plot <- ID.true;

############# PARAMETERS ##################
# initial parameter guesses

# initialization of the covariance matrices Sigma
S.init <-array(0,dim=c(p,p,K.max)) 
for(k in 1:K.max){S.init[,,k] <- 2*diag(p)} # 0.01 for random data

# random initialization of Means by sampling the dataset
# random observations are chosen to serve as initial guesses for means
Mu.arr <- array(rep(NaN, Mu.n*K.max*p),dim=c(Mu.n,K.max,p))
X.shaked <- X[sample(1:n,replace=F),]
for(i in 1:Mu.n){
  Mu.arr[i,,] <- X.shaked[sample(1:n,K.max, replace=F),]  
}

##### Start the clock!
ptm <- proc.time()

tol=1e-5; S.det.ratio=6; freq.ratio=5; make_plots=1; iter.max = 400;  # can be changed 
if (method == "EM") {EM.CEM=0; EM.CEM.stoch=0} else if (method == "CEM") {EM.CEM=1; EM.CEM.stoch=0} else if (method == "SEM") {EM.CEM=1; EM.CEM.stoch=1};
EM.output <- EM(X, K=K.fit, K.min, K.max, Mu.arr,tol,iter.max,iter.overfit,EM.CEM, EM.CEM.stoch,S.det.ratio,freq.ratio, make_plots, ID.plot, ID.true)
model.final <- EM.output[[1]]
fit_indicator <- EM.output[[2]]

#### Stop the clock
proc.time() - ptm
