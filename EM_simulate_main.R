# Fits a mixture model of K univariate normal distributions using EM algorithm
# The formulas for maximization step was derived analytically

rm(list=ls())
library(randtoolbox)
library(EMCluster)
library(MixSim)
library(ellipse)
library(e1071)
library(mvtnorm)


##### Functions ########################
source("EM_alg_mvnormal_fun.R")
source("Load_dataset.R")

### PARAMETERS THAT CAN BE CHANGED ###
sim.start <- 1
sim.n <- 50
# s = sim.start, sim.start+1, ..., sim.n
K <- 10
p <- 2
n <- 5000
AvgOmega = 0.01
Mu.n <- 30

EM.CEM <- 0  # 0: classical EM, 1: CEM
EM.CEM.stoch <- 0  # only if EM.CEM <- 1. Then if 0: CEM, if 1: SEM

tol = 1e-6
make_plots <- 0
iter.max <- 250
iter.overfit <- Inf
K.max <- K  # for initialization of tables for results

############################################################

data.simul <- list()
data.simul[[1]] <- rep(NaN,sim.n)
data.simul[[2]] <- rep(NaN,sim.n)
data.simul[[3]] <- vector("list", length = sim.n)
# data.simul[[4]] <- rep(NaN,sim.n)
# data.simul[[5]] <- rep(NaN,sim.n)
# data.simul[[6]] <- vector("list", length = sim.n)
# data.simul[[7]] <- rep(NaN,sim.n)
# data.simul[[8]] <- rep(NaN,sim.n)
# data.simul[[9]] <- vector("list", length = sim.n)

data.simul[[10]] <- rep(NaN,sim.n)

for(s in sim.start:sim.n){
  
  print(sprintf("__________________ %1.0f _____________________",s))
  
  ##### Load Data ########################
  
  # set.seed(100+s)
  gen.params <- list(n, K, p, AvgOmega)
  dataset <- load_dataset(dataset = "gen", gen.params)  # dataset: "iris", "crabs", "gen"
  X <- dataset[[1]]; ID.true <- dataset[[2]];
  n <- dim(X)[1]; p <- dim(X)[2]
  plot(X[,1], X[,2], col=ID.true, xlim=c(-0.75,1.5), ylim=c(-0.75,1.5), xlab = expression('x'[1]), ylab =  expression('x'[2]))
  
  # pairs(X,col=ID.true)
  ID.plot <- ID.true;
  
  ############# PARAMETERS ##################
  # initial conditions
  p <- dim(X)[2]
  
  S.init <-array(0,dim=c(p,p,K.max))
  for(k in 1:K.max){S.init[,,k] <- 0.02*diag(p)} # 0.01 for random data
  
  set.seed(100+s)
  Mu.arr <- array(rep(NaN, Mu.n*K.max*p),dim=c(Mu.n,K.max,p))
  X.shaked <- X[sample(1:n,replace=F),]
  for(i in 1:Mu.n){
    Mu.arr[i,,] <- X.shaked[sample(1:n,K.max, replace=F),]  
  }
  
  # #for(kk in 1:K.max){
  # #  Mu.arr[,kk,] <- halton(Mu.n, dim = p, init=T, normal=F, usetime=F, mixed=T, method="C", mexp=19937) # for each K we have D Mu's
  # #}
  # Mu.arr <- halton(K.max*p*Mu.n, dim = p, init=T, normal=F, usetime=F, mixed=T, method="C", mexp=19937) # for each K we have D Mu's
  # # Mu.arr <- runif(K.max*p*Mu.n)
  # Mu.arr <- array(Mu.arr,dim=c(Mu.n,K.max,p))
  # 
  # for(pp in 1:p){
  #   R <- max(X[,pp])-min(X[,pp])
  #   Mu.arr[,,pp] <- R*Mu.arr[,,pp]+min(X[,pp])      # scale x1 coordinate
  # }

  # pairs(Mu.arr[,1,])
  
  # par(mfrow=c(1,2))
  # plot(Mu.arr[,1,1],Mu.arr[,2,1],pch=3)
  # plot(X[,1],X[,3],pch=3)
  # 
  
  ##### Start the clock!
  ptm <- proc.time()
  
  # DEBUG
  # K.min=NaN; tol=1e-6; EM.CEM=0; EM.CEM.stoch=0; S.det.ratio=Inf; freq.ratio=10; make_plots=0;
  EM.output <- EM(X,K, K.min=NaN, K.max ,Mu.arr,tol,iter.max,iter.overfit,EM.CEM, EM.CEM.stoch,S.det.ratio=Inf,freq.ratio=Inf, make_plots,ID.plot, ID.true)

  #### Stop the clock
  stoptime <- proc.time() - ptm

  print(sprintf('!!! EM Done in %1.0f mins !!!', stoptime[3]/60))
  
  model.final <- EM.output[[1]]
  fit_indicator <- EM.output[[2]]

  data.simul[[1]][s] <- fit_indicator
  data.simul[[2]][s] <- stoptime[3]
  data.simul[[3]][[s]] <- model.final
  
  
  
# 
#   ##### Start the clock!
#   ptm <- proc.time()
# 
#   print("!!! CEM !!!")
#   # DEBUG
#   K.min=NaN; tol=1e-4; EM.CEM=1; EM.CEM.stoch=0; S.det.ratio=Inf; freq.ratio=Inf;
#   EM.output <- EM(X,K, K.min=NaN, K.max ,Mu.arr,tol=1e-4,iter.max,iter.overfit,EM.CEM=1, EM.CEM.stoch=1,S.det.ratio=Inf,freq.ratio=Inf, make_plots,ID.plot, ID.true)
# 
#   #### Stop the clock
#   stoptime <- proc.time() - ptm
# 
#   model.final <- EM.output[[1]]
#   fit_indicator <- EM.output[[2]]
# 
#   data.simul[[4]][s] <- fit_indicator
#   data.simul[[5]][s] <- stoptime[3]
#   data.simul[[6]][[s]] <- model.final
  
  data.simul[[10]][s] <- s
  
  if(EM.CEM==0){
    filename <- sprintf('EM_K_%1.0f_p_%1.0f_n_%1.0f_Omega_%0.2f_Mu.n_%1.0f_s_%1.0f-to-%1.0f.RData',K,p,n,AvgOmega,Mu.n,sim.start,s)
  }else if(EM.CEM==1 && EM.CEM.stoch == 0){
    filename <- sprintf('CEM_K_%1.0f_p_%1.0f_n_%1.0f_Omega_%0.2f_Mu.n_%1.0f_s_%1.0f-to-%1.0f.RData',K,p,n,AvgOmega,Mu.n,sim.start,s)
  }else if(EM.CEM==1 && EM.CEM.stoch == 1){
    filename <- sprintf('SEM_K_%1.0f_p_%1.0f_n_%1.0f_Omega_%0.2f_Mu.n_%1.0f_s_%1.0f-to-%1.0f.RData',K,p,n,AvgOmega,Mu.n,sim.start,s)
  }
  save(data.simul, file=filename)
  
}


