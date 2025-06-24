# Function to check if a package is installed and install it if necessary
check_and_install <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

# Check and install R packages
check_and_install("RcppHMM")
check_and_install("Rcpp")
check_and_install("pdfCluster")
check_and_install("boot")
check_and_install("xtable")

library(RcppHMM)
library(Rcpp)
library(pdfCluster)
library(boot)
library(xtable)

# Run the following only one time to install the python package
# py_install("scipy")

# # Import the python module
# import("scipy")
# 
# # Import python functions for SJM estimation
# source_python('SJ.py')

# Source cpp files
Rcpp::sourceCpp("sjm.cpp")

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

SJM_est=function(Y,K,kappa,lambda,
                 verbose=F,
                 max_iter=10,
                 tol=1e-8,
                 n_init=10,
                 Ksat=6,alpha0=NULL,K0=NULL,pers0=NULL){
  
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
  
  ## Saturated model ##
  
  res_sat=sparse_jump(Y_in=as.matrix(YY), n_states=as.integer(Ksat), 
                      max_features=sqrt(dim(YY)[2]), 
                      jump_penalty=0,
                      max_iter=max_iter, 
                      tol=tol, 
                      n_init=max_iter, 
                      verbose=F)
  
  est_weights_sat=res_sat$feat_w
  norm_est_weights_sat=est_weights_sat/sum(est_weights_sat)
  est_states_sat=order_states(res_sat$states)
  Lnsat=sum(get_BCSS(as.matrix(YY),est_states_sat,Ksat))
  
  ##
  
  
  ## Non saturated model ##
  
  res=sparse_jump(
    Y_in=as.matrix(YY),
    n_states=as.integer(K),
    max_features=kappa, 
    jump_penalty=lambda,
    max_iter=max_iter, 
    tol=tol, 
    n_init=n_init, 
    verbose=verbose
  )
  
  
  est_weights=res$feat_w
  norm_est_weights=est_weights/sum(est_weights)
  
  est_states=order_states(res$states)
  
  res_df=data.frame(state=as.factor(est_states), Y)
  
  indx=which( est_weights!=0)
  XX=YY[,indx]
  
  Ln=sum(get_BCSS(as.matrix(XX),est_states,K))
  
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

grid_search_SJM <- function(
    Y,
    K_grid=NULL,
    kappa_grid=NULL,
    lambda_grid=NULL,
    parallel   = F,
    n_cores=NULL,
    max_iter   = 10,
    n_init     = 10,
    tol        = 1e-8,
    Ksat       = 6,
    K0=NULL,
    pers0=NULL,
    alpha0=NULL,
) {
  
  library(parallel)
  
  # grid_search_SJM: Perform grid search over (K, kappa, lambda) for Sparse Jump Models
  #
  # Arguments:
  #   Y           - data matrix (N × P)
  #   K_grid      - vector of candidate numbers of states
  #   kappa_grid  - vector of candidate maximum features
  #   lambda_grid - vector of candidate jump penalties
  #   parallel    - logical; TRUE for parallel execution (default: TRUE)
  #   max_iter    - maximum iterations per model fit (default: 10)
  #   n_init      - number of random initializations per fit (default: 10)
  #   tol         - convergence tolerance (default: 1e-8)
  #   Ksat        - number of states for the saturated (zero‐penalty) model (default: 6)
  #   K0          - prior number of states 
  #   pers0       - prior for the persistence 
  #   alpha0      - prior number of relevant features 
  
  # Value:
  #   A data.frame with one row per (K, kappa, lambda) combination, containing:
  #     K           - number of states tested
  #     kappa       - maximum features used
  #     lambda      - jump penalty value
  #     FTIC        - FTIC criterion for that fit
  #     BIC         - BIC criterion for that fit
  #     AIC         - AIC criterion for that fit
  #     n_features  - number of features selected by the model
  #     n_jumps     - number of state transitions in the estimated sequence
  
  # Scale each column of Y
  YY <- apply(Y, 2, scale)
  N  <- nrow(YY)
  P  <- ncol(YY)
  
  if(is.null(K_grid)) {
    K_grid <- seq(2, 4, by = 1)  # Default range for K
  }
  
  if(is.null(kappa_grid)) {
    kappa_grid <- seq(1, sqrt(P), length.out=5)  # Default range for kappa
  }
  if(is.null(lambda_grid)) {
    lambda_grid <- c(.1,.5,1,5,10,50,100,1000)  # Default range for lambda
  }
  
  # Compute the “saturated” fit ONCE (jump_penalty = 0)
  res_sat <- sparse_jump(
    Y_in         = as.matrix(YY),
    n_states     = as.integer(Ksat),
    max_features = sqrt(P),
    jump_penalty = 0,
    max_iter     = max_iter,
    tol          = tol,
    n_init       = n_init,
    verbose      = FALSE
  )
  
  est_states_sat <- order_states(res_sat$states)
  CKp_sat        <- -log(Ksat) - ((log(2 * pi) * sqrt(P)) / 2)
  
  # Build the parameter list (one element per (K, κ, λ))
  param_grid <- expand.grid(
    K      = K_grid,
    kappa  = kappa_grid,
    lambda = lambda_grid,
    KEEP.OUT.ATTRS  = FALSE,
    stringsAsFactors = FALSE
  )
  
  # Convert each row to a named list: list(K=..., kappa=..., lambda=...)
  param_list <- apply(param_grid, 1, function(row) {
    list(
      K      = as.integer(row["K"]),
      kappa  = as.numeric(row["kappa"]),
      lambda = as.numeric(row["lambda"])
    )
  })
  
  # Define helper function to compute (FTIC, BIC, AIC) for one triple
  compute_metrics <- function(params) {
    Kval     <- params$K
    kappaval <- params$kappa
    lamval   <- params$lambda
    
    # Fit the non‐saturated model with (Kval, kappaval, lamval)
    res_non <- sparse_jump(
      Y_in         = as.matrix(YY),
      n_states     = as.integer(Kval),
      max_features = kappaval,
      jump_penalty = lamval,
      max_iter     = max_iter,
      tol          = tol,
      n_init       = n_init,
      verbose      = FALSE
    )
    
    est_w_non  <- res_non$feat_w
    est_states <- order_states(res_non$states)
    
    # Build matrix XX from selected features (where weight != 0)
    nonzero_idx <- which(est_w_non != 0)
    if (length(nonzero_idx) == 0) {
      XX <- matrix(0, nrow = N, ncol = 0)
    } else {
      XX <- YY[, nonzero_idx, drop = FALSE]
    }
    
    # Compute within‐cluster BCSS for the non‐saturated model
    Ln_non <- sum(get_BCSS(as.matrix(XX), est_states, Kval))
    
    # Count jumps in the estimated state sequence, and # features
    pen_switch <- sum(est_states[-N] != est_states[-1])
    alpha_k    <- length(nonzero_idx)
    
    # Compute CKp for non‐saturated model
    CKp_non <- -log(Kval) - ((log(2 * pi) * kappaval) / 2)
    
    # Recompute Lnsat for this Kval using the SATURATED state‐assignment
    Lnsat_for_K <- sum(get_BCSS(as.matrix(YY), est_states_sat, Ksat))
    
    # Differences Δ_Ln and Δ_CKp
    delta_Ln  <- Lnsat_for_K - Ln_non
    delta_CKp <- CKp_sat - CKp_non
    
    # Penalty components (defaults: alpha0 = P, K0 = Kval, pers0 = 0.95)
    
    if(is.null(K0)) K0 <- Kval
    if(is.null(alpha0)) alpha0 <- P
    if(is.null(pers0)) pers0 <- 0.95
    
    pen0         <- (1 - pers0) * N * (K0 - 1)
    TotalPenalty <- (alpha0 + pen0) * Kval +
      K0 * (alpha_k - alpha0 + pen_switch - pen0)
    
    # Compute FTIC, BIC, AIC
    anFTIC <- log(log(N)) * log(P)
    anAIC  <- 2
    anBIC  <- log(N)
    
    FTIC_val <- 2 * delta_CKp + (delta_Ln + anFTIC * TotalPenalty) / N
    BIC_val  <- 2 * delta_CKp + (delta_Ln + anBIC  * TotalPenalty) / N
    AIC_val  <- 2 * delta_CKp + (delta_Ln + anAIC  * TotalPenalty) / N
    
    data.frame(
      K          = Kval,
      kappa      = kappaval,
      lambda     = lamval,
      FTIC       = FTIC_val,
      BIC        = BIC_val,
      AIC        = AIC_val,
      n_features = alpha_k,
      n_jumps    = pen_switch,
      stringsAsFactors = FALSE
    )
  }
  
  # Run grid search: parallel or sequential
  if (parallel) {
    if(is.null(n_cores)){
      n_cores <- detectCores() - 1
    }
    results_raw <- mclapply(param_list, FUN = compute_metrics, mc.cores = n_cores)
  } else {
    results_raw <- lapply(param_list, FUN = compute_metrics)
  }
  
  # Combine list of data.frames into one data.frame
  results_df <- do.call(rbind, results_raw)
  
  # Return results as a data.frame (or convert to matrix if strictly needed)
  return(results_df)
}


cv_sparse_jump <- function(
    Y,
    true_states,
    K_grid=NULL,
    kappa_grid=NULL,
    lambda_grid=NULL,
    n_folds = 5,
    parallel=F,
    n_cores=NULL,
    cv_method="forward-chaining"
) {
  
  # cv_sparse_jump: Cross-validate Sparse Jump Model parameters (K, kappa, lambda)
  
  # Arguments:
  #   Y           - data matrix (N × P)
  #   true_states - vector of true states for ARI computation
  #   K_grid      - vector of candidate numbers of states
  #   kappa_grid  - vector of candidate kappas
  #   lambda_grid - vector of candidate lambdas
  #   n_folds     - number of folds for cross-validation (default: 5)
  #   parallel    - logical; TRUE for parallel execution (default: FALSE)
  #   n_cores     - number of cores to use for parallel execution (default:  NULL)
  #   cv_method   - method for cross-validation: "blocked-cv" or "forward-chain"
  
  
  # Value:
  #   A data.frame with one row per (K, kappa, lambda) combination, containing:
  #     K      - number of states tested
  #     kappa  - sparsity hyperparameter value
  #     lambda - jump penalty value
  #     ARI    - mean Adjusted Rand Index across folds
  
  
  if(is.null(K_grid)) {
    K_grid <- seq(2, 4, by = 1)  # Default range for K
  }
  
  if(is.null(kappa_grid)) {
    kappa_grid <- seq(1, sqrt(P), length.out=5)  # Default range for kappa
  }
  if(is.null(lambda_grid)) {
    lambda_grid <- c(.1,.5,1,5,10,50,100,1000)  # Default range for lambda
  }
  
  # Libreria per ARI
  library(mclust)
  
  N <- nrow(Y)
  P <- ncol(Y)
  
  # Suddivido gli N campioni in n_folds blocchi contigui
  if(method=="blocked-cv"){
  fold_indices <- split(
    1:N,
    rep(1:n_folds, each = ceiling(N / n_folds), length.out = N)
  )
  
  }
  else if(method=="forward-chain"){
    fold_indices <- lapply(seq_len(n_folds), function(k) {
      idx_end <- N - (k - 1)
      1:(idx_end-1)
    })
    names(fold_indices) <- as.character(seq_len(n_folds))
  }
  else{
    stop("cv_method must be either 'blocked-cv' or 'forward-chain'")
  }
  
  # Funzione che, per una tripla (K, kappa, lambda) e un fold (train_idx, val_idx),
  # calcola l’ARI sui punti di validazione
  fold_ari <- function(K, kappa, lambda, train_idx, val_idx) {
    # 1) Fit del modello sparse_jump su soli dati di TRAIN
    res <- sparse_jump(
      Y_in         = as.matrix(Y[train_idx, , drop = FALSE]),
      n_states     = as.integer(K),
      max_features = kappa,
      jump_penalty = lambda,
      max_iter     = 10,
      tol          = 1e-8,
      n_init       = 10,
      verbose      = FALSE
    )
    states_train <- res$states
    feat_idx     <- which(res$feat_w != 0)
    
    # Se non vengono selezionate feature, restituisco ARI = 0
    if (length(feat_idx) == 0) {
      return(0)
    }
    
    # 2) Calcolo dei centroidi (medi anche solo sulle feature selezionate)
    #    Ogni riga di "centroids" è il centro per uno stato k = 1..K
    centroids <- matrix(NA, nrow = K, ncol = length(feat_idx))
    Y_train_feats <- as.matrix(Y[train_idx, feat_idx, drop = FALSE])
    
    for (k in 1:K) {
      idx_k <- which(states_train == k)
      if (length(idx_k) == 0) {
        # se uno stato non ha punti in TRAIN, rimane NA
        centroids[k, ] <- NA
      } else {
        centroids[k, ] <- colMeans(Y_train_feats[idx_k, , drop = FALSE])
      }
    }
    
    # 3) Assegno ciascun punto in VAL a uno stato: 
    #    lo stato k che minimizza la distanza euclidea su feature selezionate
    Y_val_feats     <- as.matrix(Y[val_idx, feat_idx, drop = FALSE])
    pred_val_states <- integer(length(val_idx))
    
    for (i in seq_along(val_idx)) {
      dists <- rep(Inf, K)
      for (k in 1:K) {
        if (!any(is.na(centroids[k, ]))) {
          dists[k] <- sum((Y_val_feats[i, ] - centroids[k, ])^2)
        }
      }
      pred_val_states[i] <- which.min(dists)
    }
    
    # 4) Calcolo ARI tra etichette vere e quelle predette sul blocco VAL
    return(adjustedRandIndex(true_states[val_idx], pred_val_states))
  }
  
  # Costruisco la griglia di tutte le combinazioni di (K, kappa, lambda)
  grid <- expand.grid(
    K      = K_grid,
    kappa  = kappa_grid,
    lambda = lambda_grid,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  
  # Data.frame in cui raccogliere (K, kappa, lambda, media_ARI)
  results <- data.frame(
    K      = integer(0),
    kappa  = integer(0),
    lambda = numeric(0),
    ARI    = numeric(0)
  )
  
  # Loop su ciascuna riga della griglia
  if (!parallel) {
    # applichiamo una funzione su ogni riga di 'grid'
    results_list <- lapply(seq_len(nrow(grid)), function(row) {
      K_val     <- as.integer(grid$K[row])
      kappa_val <- as.integer(grid$kappa[row])
      lambda_val<-        grid$lambda[row]
      
      # calcolo ARI su ciascun fold
      ari_vals <- numeric(n_folds)
      for (f in seq_len(n_folds)) {
        val_idx   <- fold_indices[[f]]
        train_idx <- setdiff(seq_len(N), val_idx)
        ari_vals[f] <- fold_ari(K_val, kappa_val, lambda_val, train_idx, val_idx)
      }
      mean_ari <- mean(ari_vals)
      
      # ritorno un data.frame di una sola riga
      data.frame(
        K      = K_val,
        kappa  = kappa_val,
        lambda = lambda_val,
        ARI    = mean_ari,
        stringsAsFactors = FALSE
      )
    })
    
    # combino tutti i data.frame in un unico data.frame
    results <- do.call(rbind, results_list)
  }
  
  
  # 2) VERSIONE PARALLELA: usare mclapply()
  if (parallel) {
    if(is.null(n_cores)){
      n_cores <- detectCores() - 1
    }

    results_list <- mclapply(
      seq_len(nrow(grid)),
      function(row) {
        K_val     <- as.integer(grid$K[row])
        kappa_val <- as.integer(grid$kappa[row])
        lambda_val<-        grid$lambda[row]
        
        # calcolo ARI su ciascun fold
        ari_vals <- numeric(n_folds)
        for (f in seq_len(n_folds)) {
          val_idx   <- fold_indices[[f]]
          train_idx <- setdiff(seq_len(N), val_idx)
          ari_vals[f] <- fold_ari(K_val, kappa_val, lambda_val, train_idx, val_idx)
        }
        mean_ari <- mean(ari_vals)
        
        # ritorno un data.frame di una sola riga
        data.frame(
          K      = K_val,
          kappa  = kappa_val,
          lambda = lambda_val,
          ARI    = mean_ari,
          stringsAsFactors = FALSE
        )
      },
      mc.cores = n_cores
    )
    
    # combino tutti i data.frame in un unico data.frame
    results <- do.call(rbind, results_list)
  }
  
  return(results)
}
