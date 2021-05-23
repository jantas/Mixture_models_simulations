## Function for EM_main

library(e1071)  # for ClassAgreement()
library(numbers)
#library(plyr)
library(dplyr)
library(DescTools)

random_init <- function(K,J){
  # inputs: K (the number of clusters), J (the number of MC states)
  # output: a list of: Gamma (list of K random matrices of size JxJ), Alpha (list of initial distributions), tau (weights for MCs)
  
  Gamma <- vector(mode="list",length=K)
  Alpha <- vector(mode="list",length=K)
  tau <- random_tuple(K)
  
  for (k in 1:K){
    if (J == 2){  # probabilities must add up to 1
      rn <- runif(1)
      Alpha[[k]] <- c(rn,1-rn)
      
      Gamma[[k]] <- matrix(rep(0,J*J), byrow = T, ncol = J)
      Gamma[[k]][,1] <- runif(2)
      Gamma[[k]][,2] <- c(1,1) - Gamma[[k]][,1]
    } else if(J>2){  # for 3 and more MC states, use random_tuple to generate a vector that adds up to 1
      Alpha[[k]] <- random_tuple(J)
      
      gmat <- matrix(rep(0,J*J),ncol=J)
      for(j in 1:J){  # each row adds up to 1
        gmat[j,] <- random_tuple(J)
      }
      Gamma[[k]] <- gmat
    }
  }
  return(list(Gamma,Alpha,tau))
}

random_tuple <- function(J){
  # generate a Dirichlet vector with alpha_i = 1 thus uniform (and the vector adds up to 1)
  tuple <- vector(mode="numeric",length = J)
  rvec <- rexp(J,rate=1)
  for(i in 1:J){
    tuple[i] <- rvec[i]/sum(rvec)
  }
  return(tuple)
}  


pdf_mixed <- function(y, XX, tau, Gamma, alpha){
  # inputs: y (a list OR a matrix of MC data), XX (the list of matrces of numbers of transitions), tau (weights), 
  #         Gamma (the list of transition matrices), Alpha (the list of initial distributions)
  # output: values of the mixed model pdf (for each MC, i.e. for each row of y)
  # y can be one MC or n MCs juxtaposed in a list
  J <- length(alpha[[1]])  # the number of states
  
  # y can be a list or a matrix
  if(class(y)=="list"){
    n <- length(y)
    K <- length(tau)
    f <- matrix(rep(0,n),ncol=1)
    
    # find the matrix of numbers of transitions
    for (i in 1:n){
      for(k in 1:K){
        f[i] <- f[i] + tau[k]*alpha[[k]][y[[i]][1]]*prod(Gamma[[k]]^XX[[i]])
      }
    }
  }else{  # if y is a matrix
    n <- dim(y)[1]  # the number of MCs in y
    K <- length(tau)
    f <- matrix(rep(0,n),ncol=1)
    
    # find the matrix of numbers of transitions
    for (i in 1:n){
      for(k in 1:K){
        f[i] <- f[i] + tau[k]*alpha[[k]][y[i,1]]*prod(Gamma[[k]]^XX[[i]])
      }
    }
  }
  
  return(f)
}


get_log_like <- function(X, XX, tau, Gamma, alpha){
  # computes log-likelihood using the mixed pdf
  # inputs: X (data, ncol=K), tau (weights), Gamma (the list of transition matrices), alpha (the list of initial distributions)
  # X is the data vector
  l <- sum(log(pdf_mixed(X, XX, tau, Gamma, alpha)))
  return(l)
}

# max.index <- function(a){which(a==max(a))}

get_XX <- function(X,J) {
  # input: X data, J= number of states
  # output: XX is a list of matrices of numbers of transitions from j to j'
  
  
  if (class(X)=="list") {  # X can be a list or a matrix
    n <- length(X)
    XX <- vector(mode = "list", length = n)
    for (i in 1:n) {
      y <- X[[i]]  # ith MC
      y.T <- length(y)  # the length of ith MC
      XX_i <- matrix(rep(0,J*J),ncol=J)  # matrix JxJ of numbers of observed transitions x_ijj'
      # count the number of factors formed by taking the sequence and lagged sequence (sliding window of length two)
      y.factor <- factor(paste(c('0',as.character(y)),c(as.character(y),'0'), sep="")[2:y.T],levels=y.levels)
      y.freq <- table(y.factor)  # makes a table s.t. columns = levels, values = numbers of observation for each level

      # factors y.freq$x are for example: 11,12,21,... and y.freq$freq are the numbers of observation
      # we want to store the observations into a matrix into y.freq$x entries
      # i_index is then 1,1,2 and j_index 1,2,1
      i_index <- as.integer(floor(as.numeric(as.character(row.names(y.freq)))/10))
      j_index <- as.double(as.character(row.names(y.freq))) - 10*i_index
      for(k in 1:length(i_index)){  # store the number of observations in a matrix XX_i
        XX_i[i_index[k],j_index[k]] <- y.freq[k]  # (i,j)-entry of XX_i is the corresponding nr of obs to factor ij
      }
      XX[[i]] <- XX_i  # put the XX_i matrix into a list of XX matrices

      # y.weights <- exp(0.5*seq(-1,0,1/(y.T+1)))
      # 
      # i_index <- y[1:(y.T-1)]
      # j_index <- y[2:y.T]
      # for (j in 1:(y.T-1)) {
      #   XX_i[i_index[j], j_index[j]] <- XX_i[i_index[j], j_index[j]] + y.weights[j]
      # }
      # XX[[i]] <- XX_i  # put the XX_i matrix into a list of XX matrices
    }
    
  } else {  # if X is a matrix (the idea is the same)
    n <- dim(X)[1]  # the number of MCs in y
    XX <- vector(mode = "list", length = n)
    for(i in 1:n){
      y <- X[i,]  # ith MC
      y.T <- length(y)  # the length of ith MC
      XX_i <- matrix(rep(0,J*J),ncol=J)  # matrix JxJ of numbers of observed transitions x_ijj'
      # y.freq <- count(factor(paste(c('0',as.character(y)),c(as.character(y),'0'), sep="")[2:y.T],levels=y.levels))  # count the number of factors formed by taking the sequence and lagged sequence
      y.factor <- factor(paste(c('0',as.character(y)),c(as.character(y),'0'), sep="")[2:y.T],levels=y.levels)
      y.freq <- table(y.factor)  # makes a table s.t. columns = levels, values = numbers of observation for each level
      
      i_index <- as.integer(floor(as.numeric(as.character(row.names(y.freq)))/10))
      j_index <- as.double(as.character(row.names(y.freq))) - 10*i_index
      for(k in 1:length(i_index)){
        XX_i[i_index[k],j_index[k]] <- y.freq[k]
      }
      XX[[i]] <- XX_i
    }
  }
  return(XX)
}

get_BIC <- function(l,n,K,J){
  # computes BIC from log-likelihood l
  # n length of the data vector X
  # K number of components in the mixed model
  M <- K*(J*J-J)+K*(J-1)+(K-1)  # number of estimated parameters k transition matrices + k initial state probs + k mixing proportions
  BIC <- -2*l + M*log(n)  # -2log(L)+Mlog(n)
  return(BIC)
}

get_pi_ik <- function(X,XX, tau,Gamma, alpha){
  # computes posterior probabilities for each MC
  # TODO: using a trick for numerical stability 
  # X data matrix (each row = MC)
  # all the formulas are from Model-based biclustering of clickstream data by Volodymyr Melnykov
  # all the equation references in this method are from this article
  ###############################################################################################
  # if (class(X)=="list") {  # X can be a list or a matrix
  #   n <- length(X)  # the number of MCs in y
  #   K <- length(tau)  # the number of components
  #   J <- length(alpha[[1]])  # the number of states
  # 
  #   pi_ik <- matrix(rep(0,n*K),ncol=K)
  # 
  #   # computation of pi_ik (alt eqn (3), p.6)
  #   for (i in 1:n) {
  #     denom <- 0
  # 
  #     for (kk in 1:K) {
  #       denom <- denom + tau[kk]*alpha[[kk]][X[[i]][1]]*prod(Gamma[[kk]]^XX[[i]])
  #     }
  # 
  #     for (k in 1:K) {
  #       numer <- tau[k]*alpha[[k]][X[[i]][1]]*prod(Gamma[[k]]^XX[[i]])
  #       pi_ik[i,k] <- numer/denom
  #     }
  #   }
  # } else {  # when X is a matrix
  #   n <- dim(X)[1]
  #   K <- length(tau)
  #   J <- length(alpha[[1]])
  # 
  #   pi_ik <- matrix(rep(0,n*K),ncol=K)
  # 
  #   for(i in 1:n){
  #     denom <- 0
  # 
  #     for(kk in 1:K){
  #       denom <- denom + tau[kk]*alpha[[kk]][X[i,1]]*prod(Gamma[[kk]]^XX[[i]])
  #     }
  # 
  #     for(k in 1:K){
  #       numer <- tau[k]*alpha[[k]][X[i,1]]*prod(Gamma[[k]]^XX[[i]])
  #       pi_ik[i,k] <- numer/denom
  #     }
  #   }
  # }
  ########################################################################
  # log formula for pi_ik, p.15 for numerical stability
  # T1, T2, T3 are the three terms inside exp
  # k' is kk, j' is jj
  # I[y_{ij} = j] implemented as as.numeric(X[[i]][1] == j)
  
  if (class(X)=="list") {  # X can be a list or a matrix
    n <- length(X)  # the number of MCs in y
    K <- length(tau)  # the number of components
    J <- length(alpha[[1]])  # the number of states
    pi_ik <- matrix(rep(0,n*K),ncol=K)
    
    for (i in 1:n) {
      for (k in 1:K) {
        
        exp_formula_ik <- 0  # this is the main sum over k'
        
        for (kk in 1:K) {
          T1 <- log(tau[kk]/tau[k])
          T2 <- 0
          T3 <- 0
          for (j in 1:J) {
            # sometimes some parameters can be NaN, then break the cycle (throw it out and start over with different initial values)
            if (is.nan(alpha[[k]][j]) || is.nan(sum(Gamma[[k]]))) {break}  
            if (alpha[[k]][j] != 0){  # some values can be 0, so we cannot divide by them, if 0, skip the addition
              T2 <- T2 + as.numeric(X[[i]][1] == j) * log(alpha[[kk]][j] / alpha[[k]][j])
            }
            for (jj in 1:J) {
              if ((Gamma[[k]][j,jj] != 0)  || is.na(Gamma[[k]][j,jj])  ) {
                T3 <- T3 + XX[[i]][j, jj]*log(Gamma[[kk]][j,jj] / Gamma[[k]][j,jj]) 
              }
            }
          }
          
          exp_formula_ik <- exp_formula_ik + exp(T1 + T2 + T3)
          pi_ik[i,k] <- 1/exp_formula_ik
          
        }
      }
    }
  } else {  # when X is a matrix, for comments, look for the case when class(X) == 'list'
    n <- dim(X)[1]
    K <- length(tau)
    J <- length(alpha[[1]])
    pi_ik <- matrix(rep(0,n*K),ncol=K)
    
    for (i in 1:n) {
      for (k in 1:K) {
        
        exp_formula_ik <- 0
        
        for (kk in 1:K) {
          T1 <- log(tau[kk]/tau[k])
          T2 <- 0
          T3 <- 0
          for (j in 1:J) {
            if (is.nan(alpha[[k]][j])) {break}
            if (alpha[[k]][j] != 0){
              T2 <- T2 + as.numeric(X[i,1] == j) * log(alpha[[kk]][j] / alpha[[k]][j])
            }
            for (jj in 1:J) {
              if ((Gamma[[k]][j,jj] != 0)  || is.na(Gamma[[k]][j,jj])  ) {
                T3 <- T3 + XX[[i]][j, jj]*log(Gamma[[kk]][j,jj] / Gamma[[k]][j,jj])   
              }
            }
          }
          exp_formula_ik <- exp_formula_ik + exp(T1 + T2 + T3)
        }
        pi_ik[i,k] <- 1/exp_formula_ik
      }
    }
  }
  
  return(pi_ik)
}

estimate_params_MC <- function(X,XX,pi_ik,ID.est, EM.CEM, EM.CEM.stoch){
  # computes the mixed model parameters: tau, Gamma, alpha and also ID.est
  # ID.est... estimated membership (model membership) of observations
  # all the formulas are from Model-based biclustering of clickstream data by Volodymyr Melnykov
  # all the equation references in this method are from this article
  
  n <- dim(pi_ik)[1]
  K <- dim(pi_ik)[2]
  
  if (EM.CEM == 1) {
    # in CEM
    if (EM.CEM.stoch == 1) {
      # in SEM, the memberships ID.est are assigned randomly (following the distribution with probs pi_ik)
      
      rand <- runif(n)
      ID.est <- rep(NaN,n)
      pi_ik.cs <- t(apply(pi_ik,1,cumsum))
      
      for (k in 2:K) {
        ID.est[pi_ik.cs[,k-1] < rand & rand < pi_ik.cs[,k]] <- k
      }
      
      ID.est[rand < pi_ik.cs[,1]] <- 1
      
    } else if (EM.CEM.stoch == 0) {
      ID.est <- apply(pi_ik, 1, which.max)  # estimated ID of each observation, which.max returns the index of the max element in a vector
    }
    
    tau <- 1/n * colSums(pi_ik)
    
    for(k in 1:K){
      X.k <- X[ID.est==k,]
      XX.k <- XX[ID.est==k]
      if (is.null(dim(X.k))) {  # if there is no XX[ID.est==k], the dim is null or sometimes 0
        Gamma[[k]] <- matrix(rep(1/J,J*J),ncol=J)
        alpha[[k]] <- matrix(rep(1/J,J  ),ncol=1)
      }else{ 
        if (dim(X.k)[1]==0) {  # when it is not null but zero
          Gamma[[k]] <- matrix(rep(1/J,J*J),ncol=J)
          alpha[[k]] <- matrix(rep(1/J,J  ),ncol=1)
        } else { # in a regular case when we do not have to deal with nulls and empty or NAN
          XX.k.byrow <- matrix(unlist(XX.k),byrow = T,ncol=J*J)  # for better manipulation
          XX.k.summed <- matrix(apply(XX.k.byrow,2,sum),byrow = T, ncol=J)  # sum of matrices of transitions corresp to kth group
          XX.k.avg <- XX.k.summed/apply(XX.k.summed,1, sum)  # normalization to achieve rows to sum up to 1
          Gamma[[k]] <- XX.k.avg
          
          a<- hist(X.k[,1],breaks = 0:J, plot = F)
          alpha[[k]] <- matrix(a$density,ncol=1)  # frequencies of the first states of MCs corresp to kth group
        }
      }
      
    }
    
    return(list(tau, Gamma, alpha, ID.est))
    
  }else if(EM.CEM==0){
    
    tau <- 1/n * colSums(pi_ik)
    
    for(k in 1:K){
      alpha_denom <- sum(pi_ik[,k])
      
      if(class(X)=="list"){
        for(j in 1:J){
          alpha[[k]][j] <- sum(  pi_ik[,k]*as.numeric(as.numeric(lapply(X, `[[`, 1))==j)  )  /  alpha_denom
        }
      }else{
        for(j in 1:J){
          alpha[[k]][j] <- sum(pi_ik[,k]*as.numeric(X[,1]==j))  /  alpha_denom
        }
      }
      
      
      
      Gamma[[k]] <- matrix(rep(0,J^2),ncol=J)
      for(j in 1:J){
        for(jj in 1:J){  # jj for j'
          numer <- 0
          denom <- 0
          for(i in 1:n){
            numer <- numer + pi_ik[i,k]*XX[[i]][j,jj]
            denom <- denom + pi_ik[i,k]*sum(XX[[i]][j,])
          }
          Gamma[[k]][j,jj] <- numer/denom
        }
      }
    }
    
    return(list(tau, Gamma, alpha, NaN))
  }
}

plot_MC <- function(X,ID){
  plot(apply(X,1,mean),col=ID)
}

EM_MC <- function(X, param.init, tol,iter.max,EM.CEM, EM.CEM.stoch,ID.plot, ID.true){
  # EM_MC <- function(X, K, K.min, K.max, Mu.arr,tol,iter.max,iter.overfit,EM.CEM, EM.CEM.stoch,S.det.ratio,freq.ratio, make_plots, ID.plot, ID.true){
  # Inputs
  # X               data
  # param.init      Initial parameters. List of tau, Gamma, alpha
  # ID.plot
  # K               if it's a number run the lagorith only for that particular K (we assume we know K)
  # K.max           finding K until K.max
  # tol             tolerance for maximization of likelihood in the Maximization Step
  # EM.CEM          boolean value for inclusion of classification step
  # EM.CEM.stoch    boolean value for randomization of classification step
  # make_plots      if 1 draw plots if 0 no plots
  ###########################################
  
  # for debugging
  # K.max <- 3;
  # K <- 3
  # tol=1e-5
  # EM.CEM=0;
  # EM.CEM.stoch=0;
  # freq.ratio=5;
  # make_plots=0
  # iter.max <- 100
  make_plots <- 0
  
  tau <- param.init[[1]]
  Gamma <- param.init[[2]]
  alpha <- param.init[[3]]
  
  n <- dim(X)[1]
  if (!is.integer(n)) {n <- length(X)}
  
  K <- length(alpha)
  
  XX <- get_XX(X,J)
  
  l <- get_log_like(X, XX, tau, Gamma, alpha)
  tmp <- c(tau, unlist(Gamma), unlist(alpha), l, NaN)
  tab <- matrix(tmp, ncol=length(tmp))  # stores results of each iteration
  
  if (make_plots == 1) {
    plot_MC(X, ID.plot)
  }
  
  l.change <- 2*tol
  iter <- 0;
  ############## EM algorithm #######################
  repeat {  # EM algorithm is always increasing
    
    iter.starttime <- proc.time()
    
    if (is.nan(l.change)) {print("l is NaN"); break;}
    
    if (l.change<=tol) {
      cat(" tol break")
      break;
    }
    
    # compute posterior probabilities
    pi_ik <- get_pi_ik(X, XX, tau, Gamma, alpha)
    if (EM.CEM == 1) {
      # if(sum(is.nan(pi_ik) | pi_ik<1e-6)){print("pi_ik adjustment.")}
      pi_ik[is.nan(pi_ik) | pi_ik < 1e-6 ] <- 1e-6  # normalize to 1?
      pi_ik <- pi_ik / apply(pi_ik,1,sum)
    }
    
    # estimate tau, Gamma, alpha (maximizers of Q function)
    output <- estimate_params_MC(X,XX,pi_ik,ID.est = NaN, EM.CEM, EM.CEM.stoch)  # stoch = 1 or stoch = 0
    tau <- output[[1]]; Gamma <- output[[2]]; alpha <- output[[3]]; ID.est <- output[[4]]
    
    # only one class problem
    # if(EM.CEM==1 & mean(ID.est)%%1 == 0){print("only one class break"); break;}  # if ID.est includes only one class
    
    # compute log-likelihood
    l.old <- l
    l <- get_log_like(X, XX, tau, Gamma, alpha)
    l.change <- abs(l-l.old)/abs(l.old)
    
    tab <- rbind(tab, c(tau, unlist(Gamma), unlist(alpha),l,l.change))
    
    iter.stoptime <- proc.time() - iter.starttime
    # cat(sprintf(' %1.0f: %1.0f sec|', iter+1, round(iter.stoptime[3],0)  ))
    # cat(sprintf(' %1.0f|', iter+1))
    
    # if(sum(is.nan(Mu))>0){print("NaNs in Mu"); break;}
    
    if (iter>iter.max) {
      cat("Iter break ")
      break;
    }
    
    ID.model <- apply(pi_ik, FUN = which.max, MARGIN = 1)
    points_min_number = (J*J-2*J)+(J-1)+(K-1)
    
    for (k in 1:K) {
      if (sum(ID.model==k) < points_min_number || class(ID.model) == "list") {  # class(ID.model) is a list when there are NaNs in pi_iks
        # print('not enough points for estimation in a component')
        break
      }
    }
    
    if (make_plots == 1) {
      plot_MC(X,ID.model)
    }
    # if(to_break){break} 
    
    
    # overfitting
    # if(iter>=round(iter.overfit,0)){
    #   if(S.det.ratio<Inf){
    #     # against overfitting (if there is a component that is compared to others two small (narrow), stop)
    #     S.det <- rep(0,K)
    #     for(k in 1:K){S.det[k] <- det(S[,,k])}
    #     # print(max(S.det)/min(S.det))
    #     if(max(S.det)/min(S.det)>S.det.ratio){
    #       print("det(S) large break")
    #       break
    #     }
    #   }
    #   
    #   ID.model <- apply(pi_ik, FUN = which.max, MARGIN = 1)
    #   group.freq <- rep(0,K)
    #   for(k in 1:K){group.freq[k] <- sum(ID.model==k)}
    #   # print(max(group.freq)/min(group.freq))
    #   if(max(group.freq)/min(group.freq) > freq.ratio){  # if the biggest group contains more than 3 times the points that the smallest group contains, stop
    #     print("very small component break")
    #     break
    #   }
    # }
    
    ###################################################
    #plot_mixed_pdf(X,Mu,S,ID.plot)
    
    # compute BIC
    BIC <- get_BIC(l,n,K,J)
    
    # print(c(iter,l,BIC))
    
    iter <- iter+1
    
  }
  
  log_l <- get_log_like(X, XX, tau, Gamma, alpha)
  
  if (make_plots == 1) {
    plot_MC(X,ID.plot)
  }
  
  pi_ik <- get_pi_ik(X,XX, tau,Gamma, alpha)
  pi_ik[is.nan(pi_ik)] <- 1e-10
  
  if (!is.nan(sum(ID.true))) {  # if we know what the true class membertships are
    ID.model <- apply(pi_ik, FUN = which.max, MARGIN = 1)
    fit_indicator <- classAgreement(table(ID.true, ID.model))$crand  # calculate fit_indicator which is the adjusted Rand Index  
  } else {
    fit_indicator <- NaN  
  }
  
  model.final <- list(tau, Gamma, alpha , BIC, log_l, ID.model)
  cat(sprintf(',   %d,    %1.2f,     %1.4f,     \n',iter,BIC, fit_indicator))
  # print(c(iter,BIC, fit_indicator))
  
  return(list(model.final, fit_indicator))
}
