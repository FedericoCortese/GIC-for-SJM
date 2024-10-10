# Function to check if a package is installed and install it if necessary
check_and_install <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

# Check and install R packages
check_and_install("RcppHMM")
check_and_install("reticulate")
check_and_install("pdfCluster")
check_and_install("boot")
check_and_install("xtable")

library(RcppHMM)
library(reticulate)
library(pdfCluster)
library(boot)
library(xtable)

# Run the following only one time to install the python package
# py_install("scipy")

# Import the python module
import("scipy")

# Import python functions for SJM estimation
source_python('SJ.py')

order_states=function(states){
  
  # This function organizes states by assigning 1 to the first observed state and sequentially numbering each new state as 2, 3, etc., incrementing by 1 for each newly observed state.
  # states is a vector of observed states
  
  N=length(states)
  states_temp=rep(0,N)
  new=1
  states_temp[1]=new
  for(i in 2:N){
    if(sum(states[i]==states[1:(i-1)])==0){
      # we enter this is if-stat. whenever a new state appeares
      states_temp[i]=new+1
      new=new+1
    }
    else{
      states_temp[i]=states_temp[which(states[1:(i-1)]==states[i])[1]]
    }
  }
  return(states_temp)
}

SJM_est=function(Y,K,kappa,lambda,Ksat=6,alpha0=NULL,K0=NULL,pers0=NULL){
  
  # This function estimates the parameters of the SJM model
  
  # Arguments:
  # Y: data matrix with N times (rows) and P variables (columns)
  # K: number of states
  # kappa: maximum number of features
  # lambda: jump penalty
  # Ksat: number of states for the saturated model
  # alpha0: prior number of features
  # K0: prior number of states
  # pers0: prior for the persistence
  
  # Value:
  # A list with the following elements:
  # res_df: data frame with the estimated states and  input features
  # est_weights: estimated weights
  # est_states: estimated states
  # FTIC: FTIC criterion
  # BIC: BIC criterion
  # AIC: AIC criterion
  
  library(dplyr)
  library(tidyverse)
  library(tidyr)
  
  if(is.null(alpha0)){
    alpha0=dim(Y)[2]
  }
  
  if(is.null(K0)){
    K0=K
  }
  
  if(is.null(pers0)){
    pers0=.95
  }
  
  YY=apply(Y,2,scale)
  PP=dim(Y)[2]
  N=dim(Y)[1]
  
  res=sparse_jump(Y=as.matrix(YY), 
                  n_states=as.integer(K), 
                  max_features=kappa, 
                  jump_penalty=lambda,
                  max_iter=as.integer(10), 
                  tol=1e-4, 
                  n_init=as.integer(10), 
                  verbose=F)
  
  
  est_weights=res[[2]]
  norm_est_weights=est_weights/sum(est_weights)
  
  est_states=order_states(res[[1]])
  
  res_df=data.frame(state=as.factor(est_states), Y)
  
  indx=which( est_weights!=0)
  XX=YY[,indx]
  ### First version: compute BCSS as L1 norm
  Ln=sum(get_BCSS(as.matrix(XX),est_states))
  
  pen=sum(est_states[1:(N-1)]!=est_states[2:N])
  alphak=length(which(est_weights!=0))
  
  CKp=-log(K)-log(2*pi)*kappa/2
  CKp_sat=-log(Ksat)-log(2*pi)*sqrt(PP)/2
  
  loglik=N*CKp-Ln/2
  
  anFTIC=log(log(N))*log(PP)
  anAIC=rep(2,length(N))
  anBIC=log(N)
  
  pen0=(1-pers0)*N*(K0-1)
  
  TotalPenalty=(alpha0+pen0)*K+K0*(alphak-alpha0+pen-pen0) 
  
  res_sat=sparse_jump(Y=as.matrix(YY), n_states=as.integer(Ksat), 
                      max_features=sqrt(dim(YY)[2]), 
                      jump_penalty=0,
                      max_iter=as.integer(10), 
                      tol=1e-4, 
                      n_init=as.integer(10), 
                      verbose=F)
  
  est_weights_sat=res_sat[[2]]
  norm_est_weights_sat=est_weights_sat/sum(est_weights_sat)
  
  est_states_sat=order_states(res_sat[[1]])
  
  Lnsat=sum(get_BCSS(as.matrix(YY),est_states_sat))
  
  Ln_diff=Lnsat-Ln
  CKp_diff=CKp_sat-CKp
  
  FTIC=2*CKp_diff+(Ln_diff+anFTIC*TotalPenalty)/N
  BIC=2*CKp_diff+(Ln_diff+anBIC*TotalPenalty)/N
  AIC=2*CKp_diff+(Ln_diff+anAIC*TotalPenalty)/N
  
  return(list(res_df=res_df,
              est_weights=est_weights, 
              est_states=est_states, 
              FTIC=FTIC, 
              BIC=BIC, 
              AIC=AIC, 
              pen0=pen0, 
              pers0=pers0,
              K0=K0, 
              alpha0=alpha0, 
              K=K, 
              kappa=kappa,
              lambda=lambda,
              Ksat=Ksat))
}






sim_data=function(seed,Ktrue,N,P,cors,pers,m){
  
  # Function to simulate data from a multivariate Gaussian HMM

  # Arguments:
  # seed is the seed for the random number generator
  # Ktrue is the number of states
  # N is the number of observations
  # P is the number of features
  # cors is a vector with correlations, one for each state
  # pers is the self-transition probability 
  # m is an integer number such that m*P is the number of false features
  
  # Value:   
  # A list with the true states and the simulated data
  
  states_names = as.character(seq(1,Ktrue))
  
  Sigma <- array(0, dim =c(P,P,Ktrue))
  
  for(i in 1:Ktrue){
    Sigma[,,i] <- matrix(cors[i], ncol = P,  
                         byrow = TRUE)
    diag(Sigma[,,i])=1
  }
  
  set.seed(seed)
  Mu <- matrix(runif(P*Ktrue,min=-2,max=2), 
               nrow = P, 
               byrow = TRUE)

    A <- matrix(rep((1-pers)/(Ktrue-1),Ktrue*Ktrue), 
              ncol = Ktrue,
              byrow = TRUE)
  diag(A)=rep(pers,Ktrue)
  
  # initial probabilities
  Pi <- rep(1/Ktrue,Ktrue)
  
  HMM.cont.multi <- verifyModel(list( "Model" = "GHMM",
                                      "StateNames" = states_names,
                                      "A" = A, 
                                      "Mu" = Mu, 
                                      "Sigma" = Sigma, 
                                      "Pi" = Pi))
  set.seed(seed)
  observationSequence <- generateObservations(HMM.cont.multi, N)
  Y=t(observationSequence$Y)
  Y=apply(Y,2,scale)
  true_states=order_states(as.factor(observationSequence$X))
  
  if(m==0){
    return(list(true_states=true_states,
                YY=Y))
  }
  else{
    set.seed(seed)
    Ytil=Y[sample(N),]
    
    if(m>1){
      for(i in 2:m){
        Ytil=cbind(Ytil,Y[sample(N),])
      }
    }
    
    YY=cbind(Y,Ytil)
  }
  
  return(list(true_states=true_states,
              YY=YY))
  
}
