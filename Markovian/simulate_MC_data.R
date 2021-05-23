simulate_MC_data <- function(P,alpha,X.T,N){
  # onlu for K = 3 right now
  # P, alpha lists of transition matrices and initial distributions
  # N <- 600  # number of MC sequences (simulations)
  # T <- 20  # the length of each MC sequence
  
  P.list <- P
  alpha.list <- alpha
  
  MC_t <- function(X0,P){
    # makes a transition from current to next step
    # X0 ... current state
    # P   ... transition matrix
    # return X1 ... state in the next step
    p <- P[X0,]  # prob distribution of X1|X0 
    X1 <- rdiscrete_rv(p)
    return(X1)
  }
  
  rdiscrete_rv <- function(p){
    # returns a number from a discrete distribution
    # p ... vector of probabilites (prob distribution of a rv)
    CDF <- cumsum(p)
    U <- runif(1)
    i <- 1
    while(U>CDF[i]){
      i <- i+1
    }
    return(i)
  }
  
  ID.true <- c(rep(1,N/3), rep(2,N/3), rep(3,N/3))
  X <- matrix(rep(0,N*(X.T)),ncol=X.T)  # MC
  for(n in 1:N){
    P <- P.list[[ID.true[n]]]
    alpha <- alpha.list[[ID.true[n]]]
    
    # simulate the MC
    X[n,1] <- rdiscrete_rv(alpha)  # X0 initialization
    for(t in 2:X.T){
      X[n,t] <- MC_t(X[n,t-1],P) 
    }
  }
  
  data <- cbind(X,ID.true)
  
  return(data)
}
