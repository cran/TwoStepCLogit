
#############
## Function that computes the
## two-step estimator
## proposed in Craiu, Duchesne, Fortin, Baillargeon
##
## Calls coxph() from the library survival via
## the function Fit.models defined further below.
########
Ts.estim <- function(Donnees,Proposed.model=c(1,2),random=1:length(Proposed.model),All.m.1=FALSE,
  D="UN(1)",itermax=2000,tole=0.000001){
 # Donnees: Dataset with cluster number in column 1, stratum number
 #   in column 2, Y values in column 3, covariate values in columns
 #   4, 5, ..., p+3 
 # Proposed.model: lists the covariates that should be in the model,
 #   for instance c(1,3,4) if we want x1, x3 and x4
 # random: lists the covariates with random coefficients among the covariates
 #   listed in Proposed.model. For instance, if Proposed.model=c(1,3,4) and
 #   random=c(1,2), then x1 and x3 have random coefficients. The default is 
 #   to consider all coefficients random.   
 # All.m.1: true if sum of Y's in all strata is 1, false otherwise
 #   ... when in doubt use FALSE (always works, but slower than necessary if all
 #   stratum sums are 1
 # D: "UN" for unstructured matrix D, "UN(1)" for diagonal matrix D
 # itermax: maximal number of EM iterations
 # tole: maximal distance between successive EM iterations tolerated
 #   before declaring convergence

 Clusters <- unique(Donnees[,1]) # cluster identification
 k <- length(Clusters) # number of clusters
 p <- length(Proposed.model)  # number of beta coefficients
 q <- length(random)  # number of random coefficients
 Y <- rep(9999,k*p) # vector of beta_{ij}, i=1,...,k, j=1,...,p
 R <- matrix(0,ncol=k*p,nrow=k*p) # variance matrix of Y
 R.sim <- R
 # ML estimation separately for each cluster
 for(i in Clusters){
     Cluster.i <- Donnees[(Donnees[,1]==i),]
     Estimation <- Fit.models(Cluster.i,Proposed.model,all.m.1=All.m.1) #call to coxph()
     Y[((i-1)*p+1):(i*p)] <- Estimation[[1]]
     R[(((i-1)*p+1):(i*p)),(((i-1)*p+1):(i*p))] <- Estimation[[2]] 
     R.sim[(((i-1)*p+1):(i*p)),(((i-1)*p+1):(i*p))] <- diag(diag(Estimation[[2]]))
 }

 # code for REML estimation
 W1 <- kronecker(rep(1,k),diag(rep(1,q))) # fixed effects design matrix
 M <- diag(rep(1,k*q))-W1%*%t(W1)/k # matrix to "kill" the fixed effects
 M <- M[,(1:(q*(k-1)))]
 gammA <- t(M)%*%Y[rep(p*(0:(k-1)),each=q)+rep(random,k)] # response with fixed effects removed
 if (q==p) { RR <- R.sim } else { 
     RR <- diag(diag(R.sim)[rep(p*(0:(k-1)),each=q)+rep(random,k)]) # subset of R.sim corresponding to the random regression coefficients
 }
  
 # EM algo for either diagonal (UN(1)) or 
 # unstructured (UN) form for D
 if(D=="UN(1)"){
   # EM-algorithm 
   # Initial values
   eta.0 <- rep(1,q)
   eta.1 <- eta.0
   iterations <- 0
   abs.error <- 99999
   # E & M steps ... stops when itermax reach or tole reached
   while((iterations < itermax)*(abs.error>tole)){
      D.0 <- diag(rep(eta.0,k))
      MMRM.inv <- M%*%solve(t(M)%*%RR%*%M)
      Sigma <- solve(MMRM.inv%*%t(M)+solve(D.0))
      Mu <- Sigma%*%MMRM.inv%*%gammA
      for(j in 1:q){
          A <- diag(Sigma+Mu%*%t(Mu))
          eta.1[j] <- mean(A[((1:k)-1)*q+j])
      }
      abs.error <- max(abs(eta.0-eta.1))
      iterations <- iterations+1
      eta.0 <- eta.1
   }
   DD <- diag(rep(eta.0,k))
   DD.block <- diag(eta.0)
 }
 if(D=="UN"){
   # EM-algorithm 
   # Initial values
   D0.block <- diag(rep(1,q))
   D1.block <- D0.block
   D.0 <- kronecker(diag(rep(1,k)),D0.block)
   iterations <- 0
   abs.error <- 99999
   # E & M steps ... stops when itermax reach or tole reached
   while((iterations < itermax)*(abs.error>tole)){
      MMRM.inv <- M%*%solve(t(M)%*%RR%*%M)
      Sigma <- solve(MMRM.inv%*%t(M)+solve(D.0))
      Mu <- Sigma%*%MMRM.inv%*%gammA
      A <- Sigma + Mu%*%t(Mu)
      SUM <- 0
      for(cc in 1:k){
         Vcc <- A[(((cc-1)*q+1):(cc*q)),(((cc-1)*q+1):(cc*q))]
         SUM <- SUM + Vcc
      }
      D1.block <- SUM/k
      abs.error <- max(abs(D0.block-D1.block))
      iterations <- iterations+1
      D.0 <- kronecker(diag(rep(1,k)),D1.block)    
      #print(cat(paste(c(iterations,abs.error))))    
   }
   DD <- D.0
   DD.block <- D1.block
 }
 # Now estimation of the regression coefficients
 if (q==p) { D.sim <- DD; D.sim.block <- DD.block } else {
     D.sim.block <- matrix(0,p,p)
     D.sim.block[random,random] <- DD.block
     D.sim <- matrix(0,k*p,k*p)
     for(cc in 1:k){
          D.sim[random+(cc-1)*p,random+(cc-1)*p] <- DD[(q*(cc-1)+1):(cc*q),(q*(cc-1)+1):(cc*q)]
     }
 }
 V.sim <- R.sim + D.sim # Normally R + ZDZ', but here Z=identity ...
 Q <- kronecker(rep(1,k),diag(rep(1,p)))
 Var.Beta.sim <- solve(t(Q)%*%solve(V.sim)%*%Q)
 BetaWLS.simREML <- Var.Beta.sim%*%t(Q)%*%solve(V.sim)%*%Y
 se <- sqrt(diag(Var.Beta.sim))
 outp <- list(beta=c(BetaWLS.simREML),se=se,D=D.sim.block)
 return(outp)
}



#####
## Function that calls coxph() in the two-step method
####
Fit.models<-function(simul.object,proposed.model,all.m.1=FALSE)
{
 # simul.object: Dataset with cluster number in column 1, stratum number
 #   in column 2, Y values in column 3, covariate values in columns
 #   4, 5, ..., p+3
 # proposed.model: lists the covariates that should be in the model,
 #   for instance c(1,3,4) if we want x1, x3 and x4
 # all.m.1 is true if the sum of Y's in all strata is 1 and is false otherwise
 #   ... when in doubt use FALSE (always works, but slower than necessary if all
 #   stratum sums are 1 ...

 # First we must create dummy failure times. The 
 # censoring indicator is 1 for cases and 0 for controls, so the same as 
 # the responses. 
 cens <- simul.object[,3]
 x <- as.matrix(simul.object[,c(-1,-2,-3)])
 NN <- length(cens)
 fail.times <- 2-cens    
 SS <- max(simul.object[,2])
 Clust <- simul.object[,1]
 Strat <- (Clust-1)*SS + simul.object[,2]

 # We now fit the proposed model with coxph(), from the library survival.
 if(all.m.1==1){
     model.fit <-
     coxph(Surv(fail.times,cens)~as.matrix(x[,proposed.model])+strata(Strat)) #faster
 }
 else{
     model.fit <- 
     coxph(Surv(fail.times,cens)~as.matrix(x[,proposed.model])+strata(Strat),method="exact")  #slower
 }

 # Preparing the output necessary for the two-step method ...
 betas <- model.fit$coefficients # the cluster-level coefficient estimates
 V.indep <- model.fit$var # the R_i matrix
 output <- list(betas,V.indep)
 return(output)
}


