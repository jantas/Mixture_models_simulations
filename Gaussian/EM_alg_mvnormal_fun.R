## FunctionS for EM_main.R and EM_simulate_main.R
library(randtoolbox)
library(EMCluster)
library(MixSim)
library(ellipse)
library(e1071)
library(mclust)
library(numbers)

pdf_mixed <- function(x, tau, Mu, S){
  # inputs: x (data, ncol=K), tau, mu (list of mean vectors), Sigma (list of covar matrices)
  # output: values of the mixed model pdf
  # if x is a vector, also the output y is a vector of values corresponding to values of x
  n <- dim(x)[1]
  p <- dim(x)[2]
  if(is.null(n)){n=1; p=length(x)}
  K <- length(tau)
  y <- matrix(rep(0,n),ncol=1);
  x <- array(x,dim=c(n,p))
  for (i in 1:n){
    for(k in 1:K){
      y[i] <- y[i] + tau[k]*dmvnorm(t(x[i,]), mean=t(Mu[k,]), sigma = S[,,k])
    }
  } 
  return(y)
}


get_log_like <- function(X,tau,Mu,S){
  # computes log-likelihood of the mixed pdf
  # inputs: x (data, ncol=K), tau, mu (list of mean vectors), Sigma (list of covar matrices)
  # X is the data vector
  l <- sum(log(pdf_mixed(X, tau, Mu, S)))
  return(l)
}

# max.index <- function(a){which(a==max(a))}


get_BIC <- function(l,n,K){
  # computes BIC from log-likelihood l
  # n length of the data vector X
  # K number of components in the mixed model
  M <- 3*K-1  # number of estimated parameters (K-times mu, K-times sigma2, (K-1)-times tau since sum(tau)=1)
  BIC <- -2*l + M*log(n)  # -2log(L)+Mlog(n)
  return(BIC)
}



get_pi_ik_old <- function(X,tau,mu,S){
  # compute posterior probabilities
  # X data vector
  n <- dim(X)[1]
  K <- length(tau)
  pi_ik <- matrix(rep(0,n*K),ncol=K)
  for (i in 1:n){
    denom <- 0
    for (k in 1:K){
      denom <- denom + tau[k] * dmvnorm(t(X[i,]), t(Mu[k,]), S[,,k])
    }
    for (k in 1:K){
      pi_ik[i,k] <- (tau[k] * dmvnorm(t(X[i,]), t(Mu[k,]), S[,,k]))/denom
    }
  }
  return(pi_ik)
}



get_pi_ik <- function(X,tau,Mu,S){
  # computes posterior probabilities
  # using trick for numerical stability
  # X data vector
  n <- dim(X)[1]
  K <- length(tau)
  pi_ik <- matrix(rep(0,n*K),ncol=K)
  for (i in 1:n){
    for (k in 1:K){
      f_k <- dmvnorm(t(X[i,]), t(Mu[k,]), S[ , ,k])
      # if(f_k < 1e-6){f_k <- 1e-6}  # prevents overflowing
      for (kk in 1:K){
        f_kk <- dmvnorm(t(X[i,]), t(Mu[kk,]), S[ , ,kk])
        # if(f_kk < 1e-6){f_kk <- 1e-6}  # prevents overflowing
        pi_ik[i,k] <- pi_ik[i,k] + exp( log(tau[kk]) - log(tau[k]) + log(f_kk) - log(f_k) )
      }
      pi_ik[i,k] <- 1/pi_ik[i,k]  # might underflow
    }
  }
  return(pi_ik)
}

# estimate_params_mvnormal <- function(X,pi_ik){
#   p <- dim(X)[2]
#   n <- dim(pi_ik)[1]
#   K <- dim(pi_ik)[2]
#   tau <- 1/n * colSums(pi_ik)  # correct
#   for(k in 1:K){
#     Mu[k,] <- colSums(pi_ik[,k]*X)  /  sum(pi_ik[,k])
#     S[,,k] <- matrix(rep(0,p^2),ncol=p)
#     for(i in 1:n){
#       S[,,k] <- S[,,k] + pi_ik[i,k]*(X[i,] - Mu[k,])%*%t(X[i,] - Mu[k,])
#     }
#     S[,,k] <- S[,,k]/sum(pi_ik[,k])
#   }
#   return(list(tau, Mu, S))
# }

# # for Mu
# nom <- 0
# denom <- 0
# k <- 3
# for(i in 1:n){
#   nom <- nom + pi_ik[i,k] * X[i,] 
#   denom <- denom + pi_ik[i,k]
# }
# nom/denom


plot_mixed_pdf <- function(X,Mu,Mu.initial,S,ID){
  #TODO: better
  K <- dim(Mu)[1]
  p <- dim(Mu)[2]
  # if(length(arg.ID)>1){
  #   ID <- arg.ID  # estimated (for CEM)
  # }else{
  #   ID <- plot.ID
  # }
  # par(mfrow=c(p,p))
  # for(j in 1:p){
  #   for(i in 1:p){
  #     if(i == j){
  #       plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
  #       text(x = 0.5, y = 0.5, paste(c("X    ",i)), 
  #            cex = 1.6, col = "black")
  #     }else{
  #       plot(X[,i], X[,j], col=trueID, xlab = paste("X", i, sep = ""), ylab = paste("X", j, sep = ""))
  #       for(k in 1:K){
  #         points(ellipse(x = S[c(i,j),c(i,j),k], centre = Mu[k,c(i,j)], level=0.95), type="l", col=k)  # contains 95% points
  #       }
  #     }
  #   }
  # }
  
  if(p==2){par(mfrow=c(1,1))
  }else{
    p.div <- divisors(choose(p,2))
    plot.n <- p.div[length(p.div)-1]
    plot.m <- choose(p,2)/plot.n
    par(mfrow=c(plot.m, plot.n))
  }
  
  for(i in 1:p){
    for(j in (i+1):p){
      if(j>p){break}
      plot(X[,i], X[,j], col=ID, xlab = expression('x'[1]), ylab = expression('x'[2]),xlim=c(min(X[,i]),max(X[,i])),ylim=c(min(X[,j]),max(X[,j])))
      for(k in 1:K){
        points(ellipse(x = S[c(i,j),c(i,j),k], centre = Mu[k,c(i,j)], level=0.95), type="l", col=1)  # contains 95% points
        if(!is.nan(Mu.initial)){points(Mu.initial[k,i],Mu.initial[k,j], pch=8, col=5, cex=3)}
      }
    }
  }
}



estimate_params_mvnormal <- function(X,pi_ik,ID.est, EM.CEM, EM.CEM.stoch){
  
  p <- dim(X)[2]
  n <- dim(pi_ik)[1]
  K <- dim(pi_ik)[2]
  
  Mu <- matrix(rep(NaN,K*p),ncol = p)
  S <-  array(0,dim=c(p,p,K))
  
  if(EM.CEM == 1){
    if(EM.CEM.stoch == 1){
      
      rand <- runif(n)
      ID.est <- rep(NaN,n)
      pi_ik.cs <- t(apply(pi_ik,1,cumsum))
      
      for(k in 2:K){
        ID.est[pi_ik.cs[,k-1]<rand & rand<pi_ik.cs[,k]] <- k
      }
      
      ID.est[rand<pi_ik.cs[,1]] <- 1
      
    }else if(EM.CEM.stoch == 0){
      ID.est <- apply(pi_ik,1,which.max)  # estimated ID of each observation, which.max returns the index of the max element in a vector
    }
    
    tau <- 1/n * colSums(pi_ik)
    
    for(k in 1:K){
      X.k <- X[ID.est==k,]
      if(is.null(dim(X.k))){
        Mu[k,] <- X.k
        S[,,k] <- diag(p)
      }else{
        Mu[k,] <- apply(X.k,2,mean)
        S[,,k] <- var(X.k)
      }
      
    }
    
    return(list(tau, Mu, S, ID.est))
    
  }else if(EM.CEM==0){
    
    tau <- 1/n * colSums(pi_ik)
    
    for(k in 1:K){
      Mu[k,] <- colSums(pi_ik[,k]*X)  /  sum(pi_ik[,k])
      S[,,k] <- matrix(rep(0,p^2),ncol=p)
      for(i in 1:n){
        S[,,k] <- S[,,k] + pi_ik[i,k]*(X[i,] - Mu[k,])%*%t(X[i,] - Mu[k,])
      }
      S[,,k] <- S[,,k]/sum(pi_ik[,k])
    }
    
    return(list(tau, Mu, S, NaN))
  }
}

EM <- function(X, K, K.min, K.max, Mu.arr,tol,iter.max,iter.overfit,EM.CEM, EM.CEM.stoch,S.det.ratio,freq.ratio, make_plots, ID.plot, ID.true){
  # Inputs
  # X               data
  # ID.plot
  # K               if it's a number run the algorithm only for that particular K (we assume we know K)
  # K.max           finding K until K.max
  # Mu.arr          initial points for means of ellipsoids
  # tol             tolerance for maximization of likelihood in the Maximization Step
  # EM.CEM          boolean value for inclusion of classification step
  # EM.CEM.stoch    boolean value for randomization of classification step
  # S.det.ratio     stop when the ellipsoid is too narrow, iris 6
  # freq.ratio      stop when the biggest group contains more than 'freq.ratio' times the points that the smallest group contains
  # make_plots      if 1 draw plots if 0 no plots
  ###########################################
  
  if(!is.nan(K)){
    
    n <- dim(X)[1]; p <- dim(X)[2]; Mu.n <- length(Mu.arr)
    
    tab.main <- matrix(rep(NaN,((K*(1+p+p^2))+3)*Mu.n), ncol=K*(1+p+p^2)+3)  # stores results for each rand mu
    col.dim <- dim(tab.main)[2]
    Mu.n <- dim(Mu.arr)[1]
    
    for(ii in 1:Mu.n){  # loop for different mu's
      
      Mu <- Mu.arr[ii,1:K,]
      S <- S.init[,,1:K]
      tau <- rep(1/K,K)
      l <- get_log_like(X,tau,Mu,S)
      tmp <- c(tau, as.numeric(Mu), as.numeric(S),l,NaN)
      tab <- matrix(tmp,ncol=length(tmp))  # stores results of each iteration
      
      if(make_plots == 1){
        plot_mixed_pdf(X,Mu,Mu.arr[ii,1:K,],S,ID.plot)
      }
      
      l.change <- 2*tol
      
      #print(sprintf('SIM #%1.0f - MEAN #%1.0f: ',s,ii))
      
      iter <- 0;
      ############## EM algorithm #######################
      repeat{  # EM algorithm is always increasing
        
        iter.starttime <- proc.time()
        
        if(l.change<=tol){
          print("tol break")
          break;
        }
        
        # compute posterior probabilities
        pi_ik <- get_pi_ik(X,tau,Mu,S)
        # if(sum(is.nan(pi_ik) | pi_ik<1e-6)){print("pi_ik adjustment.")}
        pi_ik[is.nan(pi_ik) | pi_ik<1e-6 ] <- 1e-6  # normalize to 1?
        pi_ik <- pi_ik/apply(pi_ik,1,sum)
        
        # estimate tau, mu, sigma2 by ML method
        output <- estimate_params_mvnormal(X, pi_ik, ID.est, EM.CEM, EM.CEM.stoch)  # stoch = 1 or stoch = 0
        tau <- output[[1]]; Mu <- output[[2]]; S <- output[[3]]; ID.est <- output[[4]]
        
        # only one class problem
        # if(EM.CEM==1 & mean(ID.est)%%1 == 0){print("only one class break"); break;}  # if ID.est includes only one class
        
        # compute log-likelihood
        l.old <- l
        l <- get_log_like(X,tau,Mu,S)
        l.change <- abs(l-l.old)/abs(l.old)
        
        tab <- rbind(tab, c(tau, as.numeric(Mu), as.numeric(S),l,l.change))
        
        iter.stoptime <- proc.time() - iter.starttime
        cat(sprintf(' %1.0f: %1.0f sec|', iter+1, round(iter.stoptime[3],0)  ))
        # cat(sprintf(' %1.0f|', iter+1))
        
        if(make_plots == 1){
          if(EM.CEM == 1){plot_mixed_pdf(X,Mu,Mu.arr[ii,1:K,],S,ID.est)}else{plot_mixed_pdf(X,Mu,Mu.arr[ii,1:K,],S,ID.plot)}
        }
        
        if(sum(is.nan(Mu))>0){print("NaNs in Mu"); break;}
        
        if(iter>iter.max){
          print("Iter break")
          break;
        }
        
        
        ID.model <- apply(pi_ik, FUN = which.max, MARGIN = 1)
        points_min_number = p + (p+1)*p/2 + K
        to_break <- 0
        for(k in 1:K){
          if(sum(ID.model==k) < points_min_number){
            print('not enough points for estimation in a component break')
            to_break <- 1
            break
          }
        }
        if(to_break){break} 
        
        
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
        BIC <- get_BIC(l,n,K)
        
        # print(c(K,l,BIC))
        
        # update tab.main
        tab.main[ii,1:K] <- tau
        tab.main[ii,(K+1):(p*K+K)] <- Mu
        tab.main[ii,((p+1)*K+1):((p+1)*K+p^2*K)] <- S
        tab.main[ii,((1+p+p^2)*K+1):((1+p+p^2)*K+3)] <- c(K, l,BIC)
        
        iter <- iter+1
        
      }
    }
    
  }else{
    
    n <- dim(X)[1]; p <- dim(X)[2]
    
    tab.main <- matrix(rep(NaN,((K.max*(1+p+p^2))+3)*Mu.n), ncol=K.max*(1+p+p^2)+3)  # stores results for each rand mu
    col.dim <- dim(tab.main)[2]
    Mu.n <- dim(Mu.arr)[1]
    
    for(ii in 1:Mu.n){  # loop for different mu's
      # K <- 1;
      # # When K is 1, only MLE est #TODO: always the same ==> move outside the loop
      # tau <- 1
      # Mu <- colMeans(X)
      # Mu <- matrix(Mu,ncol=p)
      # S <- var(X)*(n-1)/n  # MLE var matrix
      # S <-array(S,dim=c(p,p,1))
      # 
      # l <- get_log_like(X,tau,Mu,S)
      # BIC <- get_BIC(l,n,K)
      # 
      # if(make_plots == 1){
      #   plot_mixed_pdf(X,Mu,S,ID.plot)
      #   print(c(ii,'-------------'))
      # }
      
      K <- K.min
      BIC <- Inf
      
      # when K > 1, do EM alg
      repeat{  # loop to find the best K w.r.t. BIC (smallest), while BIC gets smaller
        
        Mu <- Mu.arr[ii,1:K,]
        S <- S.init[,,1:K]
        tau <- rep(1/K,K)
        l <- get_log_like(X,tau,Mu,S)
        tmp <- c(tau, as.numeric(Mu), as.numeric(S),l,NaN)
        tab <- matrix(tmp,ncol=length(tmp))  # stores results of each iteration
        
        if(make_plots == 1){
          plot_mixed_pdf(X,Mu,Mu.arr[ii,1:K,],S,ID.plot)
        }
        
        l.change <- 2*tol
        iter <- 0;
        ############## EM algorithm #######################
        repeat{  # EM algorithm is always increasing
          
          if(l.change<=tol){
            print("tol break")
            break;
          }
          
          # compute posterior probabilities
          pi_ik <- get_pi_ik(X,tau,Mu,S)
          # if(sum(is.nan(pi_ik) | pi_ik<1e-6)){print("pi_ik adjustment.")}
          pi_ik[is.nan(pi_ik) | pi_ik<1e-6 ] <- 1e-6  # normalize to 1?
          pi_ik <- pi_ik/apply(pi_ik,1,sum)
          
          # estimate tau, mu, sigma2 by ML method
          output <- estimate_params_mvnormal(X, pi_ik, ID.est, CEM=EM.CEM, stoch=EM.CEM.stoch)  # stoch = 1 or stoch = 0
          tau <- output[[1]]; Mu <- output[[2]]; S <- output[[3]]; ID.est <- output[[4]]
          
          # to avoid problem that all points are assigned to only one class 
          # if(mean(ID.est)%%1 == 0){break;}  # if ID.est includes only one class
          
          # compute log-likelihood
          l.old <- l
          l <- get_log_like(X,tau,Mu,S)
          l.change <- abs(l-l.old)/abs(l.old)
          
          tab <- rbind(tab, c(tau, as.numeric(Mu), as.numeric(S),l,l.change))
          
          if(sum(is.nan(Mu))>0){print("NaNs in Mu"); break;}
          
          if(make_plots == 1){
            if(EM.CEM == 1){plot_mixed_pdf(X,Mu,Mu.arr[ii,1:K,],S,ID.est)}else{plot_mixed_pdf(X,Mu,Mu.arr[ii,1:K,],S,ID.plot)}
          }
          
          if(iter>30){
            print("Iter break")
            break;
          }
          iter <- iter+1
          
        }
        ###################################################
        #plot_mixed_pdf(X,Mu,S,ID.plot)
        
        # compute BIC
        BIC.old <- BIC
        BIC <- get_BIC(l,n,K)
        
        # print(c(K,l,BIC))
        
        # Stopping criterion
        if(BIC >= BIC.old){  # we want minimal BIC, if new BIC < BIC.old, continue. If otherwise, break
          print("BIC break")
          break
        }
        
        
        # against overfitting (if there is a component that is compared to others two small (narrow), stop)
        S.det <- rep(0,K)
        for(k in 1:K){S.det[k] <- det(S[,,k])}
        # print(max(S.det)/min(S.det))
        if(K>3 && max(S.det)/min(S.det)>S.det.ratio){
          print("det(S) large break")
          break
        }
        
        ID.model <- apply(pi_ik, FUN = which.max, MARGIN = 1)
        group.freq <- rep(0,K)
        for(k in 1:K){group.freq[k] <- sum(ID.model==k)}
        # print(max(group.freq)/min(group.freq))
        if(max(group.freq)/min(group.freq) > freq.ratio){  # if the biggest group contains more than 3 times the points that the smallest group contains, stop
          print("very small component break")
          break
        }
        
        # update tab.main (K can be different)
        tab.main[ii,1:K] <- tau
        tab.main[ii,(K.max+1):(p*K+K.max)] <- Mu
        tab.main[ii,((p+1)*K.max+1):((p+1)*K.max+p^2*K)] <- S
        tab.main[ii,((1+p+p^2)*K.max+1):((1+p+p^2)*K.max+3)] <- c(K, l,BIC)
        
        K <- K + 1
        
      }
      
      
    }
  }
  ########################## END of EM #############################
  
  # tab.main
  tab.main <- tab.main[!is.nan(tab.main[,dim(tab.main)[2]]),]
  
  # my.max <- function(x) ifelse( !all(is.na(x)), max(x, na.rm=T), NA)
  # my.min <- function(x) ifelse( !all(is.na(x)), min(x, na.rm=T), NA)
  # 
  # S <- tab.main[,(2*K.max+1):(3*K.max)]  # all sigmas (including NaNs)
  # sigmas <- matrix(sigmas,byrow = F,ncol=K.max)
  # sigmas.max <- matrix(apply(sigmas,1,my.max),ncol=1)
  # sigmas.min <- matrix(apply(sigmas,1,my.min),ncol=1)
  # # cbind(sigmas, sigmas.max, sigmas.min)  # check
  # sigmas.ratio <- sigmas.max/sigmas.min
  
  # tab.main.restr <- tab.main[sigmas.ratio<sigmas.ratio.limit, ]  # restricts solution to sigmas than are not too different
  # tab.main.minBIC <- tab.main.restr[tab.main.restr[,3*K.max+3]==min(tab.main.restr[,3*K.max+3]),]
  tab.main.minBIC <- tab.main[tab.main[,col.dim]==min(tab.main[,col.dim]),]
  
  # if(dim(tab.main.minBIC)[1]>1){
  #   print("More than one optimum!")
  #   tab.main.minBIC <- tab.main.minBIC[1,]
  # }
  
  K <- tab.main.minBIC[col.dim-2]
  tau <- tab.main.minBIC[1:K]
  Mu <- matrix(tab.main.minBIC[(K.max+1):(p*K+K.max)],byrow =F, ncol=p)
  S <- array(tab.main.minBIC[((p+1)*K.max+1):((p+1)*K.max+p^2*K)], dim=c(p,p,K))
  BIC <- tab.main.minBIC[col.dim]
  log_l <- get_log_like(X,tau,Mu,S)
  
  if(make_plots == 1){
    plot_mixed_pdf(X,Mu,Mu.arr[ii,1:K,],S, ID.plot)
  }
  
  pi_ik <- get_pi_ik(X,tau,Mu,S)
  pi_ik[is.nan(pi_ik)] <- 1e-10
  
  ID.model <- apply(pi_ik, FUN = which.max, MARGIN = 1)
  fit_indicator <- classAgreement(table(ID.true,ID.model))$crand
  
  model.final <- list(K, tau, Mu, S, BIC, log_l, ID.model)
  
  return(list(model.final, fit_indicator))
}
