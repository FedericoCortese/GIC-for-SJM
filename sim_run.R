library(RcppHMM)
library(reticulate)
library(pdfCluster)
library(boot)
library(xtable)
#py_install("scipy")
import("scipy")
source_python('SJ.py')

order_states=function(states){
  
  # This function orders states in the following way: we denote with 1  
  # the first observed state, the subsequent with "new" = 2,3,..., with 
  # "new" increasing by 1 every time we observe a new state
  
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

# Function to simulate data
sim_data=function(seed,Ktrue,N,P,cors,pers,m){
  
  # cors is a vector with correlations, one for each state
  # pers is the self-transition probability 
  
  states_names = as.character(seq(1,Ktrue))
  # Ktrue=length(N)
  #P <- 100
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
  # a=pers
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
  
  # False features generation
  
  # P is the number of true features
  # m*P is the number of false features
  
  set.seed(seed)
  Ytil=Y[sample(N),]
  for(i in 2:m){
    Ytil=cbind(Ytil,Y[sample(N),])
  }
  
  YY=cbind(Y,Ytil)
  
  return(list(true_states=true_states,
              YY=YY))
  
}

# Simulate 100 datasets (one for each T=300,600,1000) OUTSIDE any function
# seed=1:100
# YYs=lapply(seed,sim_data,Ktrue=2,N=1000,P=100,cors=c(.8,.4),pers=.8,m=2)
# true_states=simdat$true_states
# YY=simdat$YY

sim_est_SJM=function(seed,K,lambda,kappa,true_data,Ktrue,
                     N=1000,P=100,m=2,cors,pers,K0,alpha0){
  
  # true_data is obtained using the function "sim_data"
  
  YY=true_data[[seed]]$YY
  true_states=true_data[[seed]]$true_states
  
  res=sparse_jump(Y_in=as.matrix(YY), n_states=as.integer(K), 
                  max_features=kappa, 
                  jump_penalty=lambda,
                  max_iter=as.integer(10), 
                  tol=1e-4, 
                  n_init=as.integer(10), 
                  verbose=F)
  
  est_weights=res$feat_w
  norm_est_weights=est_weights/sum(est_weights)
  
  est_states=order_states(res$states)
  true_states=order_states(true_states)
  
  Pfalse=dim(YY)[2]-P
  true_weights_seq=c(rep(1,P),rep(0,Pfalse))
  est_weights_seq=as.numeric(est_weights!=0)
  
  ARI_states=adj.rand.index(true_states,est_states)
  ARI_weights=adj.rand.index(true_weights_seq,est_weights_seq)
  
  corr_sel_feat=sum(est_weights[1:P]!=0)
  wrong_sel_feat=sum(est_weights[-(1:P)]!=0)
  
  indx=which( est_weights!=0)
  XX=YY[,indx]
  ### First version: compute BCSS as L1 norm
  Ln1=sum(get_BCSS(as.matrix(XX),est_states,K))
  ### Second version: compute BCSS as L2 norm
  Ln2=sqrt(sum(get_BCSS(as.matrix(XX),est_states,K)^2))
  
  pen=sum(est_states[1:(N-1)]!=est_states[2:N])
  alphak=length(which(est_weights!=0))
  
  TotalPenalty=K0*alpha0+alpha0*(K-K0)+K0*(alphak-alpha0)+pen
  
  return(list(#YY=YY,
    indx=indx,
    est_weights=est_weights,
    norm_est_weights=norm_est_weights,
    est_weights_seq=est_weights_seq,
    est_states=est_states,
    true_states=true_states,
    seed=seed,
    N=N,
    P=P,
    m=m,
    cors=cors,
    pers=pers,
    K=K,
    lambda=lambda,
    kappa=kappa,
    Ln1=Ln1,
    Ln2=Ln2,
    TotalPenalty=TotalPenalty,
    pen=pen,
    alphak=alphak,
    ARI_states=ARI_states,
    ARI_weights=ARI_weights,
    corr_sel_feat=corr_sel_feat,
    wrong_sel_feat=wrong_sel_feat))
  
}

satmod_est=function(seed,true_data,Ksat=6){
  
  YY=true_data[[seed]]$YY
  true_states=true_data[[seed]]$true_states
  
  res_sat=sparse_jump(Y_in=as.matrix(YY), n_states=as.integer(Ksat), 
                      max_features=sqrt(dim(YY)[2]), 
                      jump_penalty=0,
                      max_iter=as.integer(10), 
                      tol=1e-4, 
                      n_init=as.integer(10), 
                      verbose=F)
  
  est_weights_sat=res_sat$feat_w
  norm_est_weights_sat=est_weights_sat/sum(est_weights_sat)
  
  est_states_sat=order_states(res_sat$states)
  
  ### First version: compute BCSS as L1 norm
  Lnsat1=sum(get_BCSS(as.matrix(YY),est_states_sat,K))
  ### Second version: compute weighted BCSS (some features will disappear)
  Lnsat2=sqrt(sum(get_BCSS(as.matrix(YY),est_states_sat,K)^2))
  
  return(list(Lnsat1=Lnsat1,
              Lnsat2=Lnsat2,
              est_weights_sat=est_weights_sat,
              norm_est_weights_sat=norm_est_weights_sat,
              est_states_sat=est_states_sat))
  
}

P=100
m=2
Ns=c(300,600,1000)
seed=1:100
K=2:4
#kappa=seq(1,floor(sqrt(PP)))
kappa=c(seq(1,10,by=1),12,14,17) # 14 values
lambda=c(0,5,10,25,50,100) # 6 values
hp=expand.grid(seed=seed,K=K,kappa=kappa,lambda=lambda)
hpsat=expand.grid(seed=seed,notimp=1)
#pers1=.8


# pers 90 -----------------------------------------------------------------
pers1=.9

# Ktrue=2 -----------------------------------------------------------------
corsK2=c(.8,.4)
YYs_K2_300=lapply(seed,sim_data,Ktrue=2,N=Ns[1],P=100,cors=corsK2,pers=pers1,m=2)
YYs_K2_600=lapply(seed,sim_data,Ktrue=2,N=Ns[2],P=100,cors=corsK2,pers=pers1,m=2)
YYs_K2_1000=lapply(seed,sim_data,Ktrue=2,N=Ns[3],P=100,cors=corsK2,pers=pers1,m=2)

# T=300
start=Sys.time()
gicsim_K2_T300 <- parallel::mclapply(1:nrow(hp),
                                     function(x)
                                       sim_est_SJM(seed=hp[x,]$seed,
                                                   K=hp[x,]$K,
                                                   lambda=hp[x,]$lambda,
                                                   kappa=hp[x,]$kappa,
                                                   true_data=YYs_K2_300,
                                                   Ktrue=2,
                                                   N=Ns[1],
                                                   P=100,
                                                   m=2,
                                                   cors=corsK2,
                                                   pers=pers1,
                                                   K0=2,
                                                   alpha0=100),
                                     mc.cores = parallel::detectCores())
gicsim_sat_K2_T300<-parallel::mclapply(1:nrow(hpsat),
                                       function(x) 
                                         satmod_est(seed=hpsat[x,]$seed,
                                                    true_data = YYs_K2_300,
                                                    Ksat = 6),
                                       mc.cores = parallel::detectCores())
end=Sys.time()
end-start
save(gicsim_K2_T300,file="GIC_K2_T300_pers90.Rdata")
save(gicsim_sat_K2_T300,file="GICsat_K2_T300_pers90.Rdata")




# T=600
start=Sys.time()
gicsim_K2_T600 <- parallel::mclapply(1:nrow(hp),
                                     function(x)
                                       sim_est_SJM(seed=hp[x,]$seed,
                                                   K=hp[x,]$K,
                                                   lambda=hp[x,]$lambda,
                                                   kappa=hp[x,]$kappa,
                                                   true_data=YYs_K2_600,
                                                   Ktrue=2,
                                                   N=Ns[2],
                                                   P=100,
                                                   m=2,
                                                   cors=corsK2,
                                                   pers=pers1,
                                                   K0=2,
                                                   alpha0=100),
                                     mc.cores = parallel::detectCores())
gicsim_sat_K2_T600<-parallel::mclapply(1:nrow(hpsat),
                                       function(x) 
                                         satmod_est(seed=hpsat[x,]$seed,
                                                    true_data = YYs_K2_600,
                                                    Ksat = 6),
                                       mc.cores = parallel::detectCores())
end=Sys.time()
end-start
save(gicsim_K2_T600,file="GIC_K2_T600_pers90.Rdata")
save(gicsim_sat_K2_T600,file="GICsat_K2_T600_pers90.Rdata")

# T=1000
start=Sys.time()
gicsim_K2_T1000 <- parallel::mclapply(1:nrow(hp),
                                      function(x) 
                                        sim_est_SJM(seed=hp[x,]$seed,
                                                    K=hp[x,]$K,
                                                    lambda=hp[x,]$lambda,
                                                    kappa=hp[x,]$kappa,
                                                    true_data=YYs_K2_1000,
                                                    Ktrue=2,
                                                    N=Ns[3],
                                                    P=100,
                                                    m=2,
                                                    cors=corsK2,
                                                    pers=pers1,
                                                    K0=2,
                                                    alpha0=100),
                                      mc.cores = parallel::detectCores())
gicsim_sat_K2_T1000<-parallel::mclapply(1:nrow(hpsat),
                                        function(x) 
                                          satmod_est(seed=hpsat[x,]$seed,
                                                     true_data = YYs_K2_1000,
                                                     Ksat = 6),
                                        mc.cores = parallel::detectCores())
end=Sys.time()
end-start
save(gicsim_K2_T1000,file="GIC_K2_T1000_pers90.Rdata")
save(gicsim_sat_K2_T1000,file="GICsat_K2_T1000_pers90.Rdata")


# Ktrue=3 -----------------------------------------------------------------
corsK3=c(.8,.6,.3)
YYs_K3_300=lapply(seed,sim_data,Ktrue=3,N=Ns[1],P=100,cors=corsK3,pers=pers1,m=2)
YYs_K3_600=lapply(seed,sim_data,Ktrue=3,N=Ns[2],P=100,cors=corsK3,pers=pers1,m=2)
YYs_K3_1000=lapply(seed,sim_data,Ktrue=3,N=Ns[3],P=100,cors=corsK3,pers=pers1,m=2)

# T=300
start=Sys.time()
gicsim_K3_T300 <- parallel::mclapply(1:nrow(hp),
                                     function(x) 
                                       sim_est_SJM(seed=hp[x,]$seed,
                                                   K=hp[x,]$K,
                                                   lambda=hp[x,]$lambda,
                                                   kappa=hp[x,]$kappa,
                                                   true_data=YYs_K3_300,
                                                   Ktrue=3,
                                                   N=Ns[1],
                                                   P=100,
                                                   m=2,
                                                   cors=corsK3,
                                                   pers=pers1,
                                                   K0=3,
                                                   alpha0=100),
                                     mc.cores = parallel::detectCores())
gicsim_sat_K3_T300<-parallel::mclapply(1:nrow(hpsat),
                                       function(x) 
                                         satmod_est(seed=hpsat[x,]$seed,
                                                    true_data = YYs_K3_300,
                                                    Ksat = 6),
                                       mc.cores = parallel::detectCores())
end=Sys.time()
end-start
save(gicsim_K3_T300,file="GIC_K3_T300_pers90.Rdata")
save(gicsim_sat_K3_T300,file="GICsat_K3_T300_pers90.Rdata")




# T=600
start=Sys.time()
gicsim_K3_T600 <- parallel::mclapply(1:nrow(hp),
                                     function(x) 
                                       sim_est_SJM(seed=hp[x,]$seed,
                                                   K=hp[x,]$K,
                                                   lambda=hp[x,]$lambda,
                                                   kappa=hp[x,]$kappa,
                                                   true_data=YYs_K3_600,
                                                   Ktrue=3,
                                                   N=Ns[2],
                                                   P=100,
                                                   m=2,
                                                   cors=corsK3,
                                                   pers=pers1,
                                                   K0=3,
                                                   alpha0=100),
                                     mc.cores = parallel::detectCores())
gicsim_sat_K3_T600<-parallel::mclapply(1:nrow(hpsat),
                                       function(x) 
                                         satmod_est(seed=hpsat[x,]$seed,
                                                    true_data = YYs_K3_600,
                                                    Ksat = 6),
                                       mc.cores = parallel::detectCores())
end=Sys.time()
end-start
save(gicsim_K3_T600,file="GIC_K3_T600_pers90.Rdata")
save(gicsim_sat_K3_T600,file="GICsat_K3_T600_pers90.Rdata")

# T=1000
start=Sys.time()
gicsim_K3_T1000 <- parallel::mclapply(1:nrow(hp),
                                      function(x) 
                                        sim_est_SJM(seed=hp[x,]$seed,
                                                    K=hp[x,]$K,
                                                    lambda=hp[x,]$lambda,
                                                    kappa=hp[x,]$kappa,
                                                    true_data=YYs_K3_1000,
                                                    Ktrue=3,
                                                    N=Ns[3],
                                                    P=100,
                                                    m=2,
                                                    cors=corsK3,
                                                    pers=pers1,
                                                    K0=3,
                                                    alpha0=100),
                                      mc.cores = parallel::detectCores())
gicsim_sat_K3_T1000<-parallel::mclapply(1:nrow(hpsat),
                                        function(x) 
                                          satmod_est(seed=hpsat[x,]$seed,
                                                     true_data = YYs_K3_1000,
                                                     Ksat = 6),
                                        mc.cores = parallel::detectCores())
end=Sys.time()
end-start
save(gicsim_K3_T1000,file="GIC_K3_T1000_pers90.Rdata")
save(gicsim_sat_K3_T1000,file="GICsat_K3_T1000_pers90.Rdata")


# Ktrue=4 -----------------------------------------------------------------
corsK4=c(.8,.6,.3,0)
YYs_K4_300=lapply(seed,sim_data,Ktrue=4,N=Ns[1],P=100,cors=corsK4,pers=pers1,m=2)
YYs_K4_600=lapply(seed,sim_data,Ktrue=4,N=Ns[2],P=100,cors=corsK4,pers=pers1,m=2)
YYs_K4_1000=lapply(seed,sim_data,Ktrue=4,N=Ns[3],P=100,cors=corsK4,pers=pers1,m=2)

# T=300
start=Sys.time()
gicsim_K4_T300 <- parallel::mclapply(1:nrow(hp),
                                     function(x) 
                                       sim_est_SJM(seed=hp[x,]$seed,
                                                   K=hp[x,]$K,
                                                   lambda=hp[x,]$lambda,
                                                   kappa=hp[x,]$kappa,
                                                   true_data=YYs_K4_300,
                                                   Ktrue=4,
                                                   N=Ns[1],
                                                   P=100,
                                                   m=2,
                                                   cors=corsK4,
                                                   pers=pers1,
                                                   K0=4,
                                                   alpha0=100),
                                     mc.cores = parallel::detectCores())
gicsim_sat_K4_T300<-parallel::mclapply(1:nrow(hpsat),
                                       function(x) 
                                         satmod_est(seed=hpsat[x,]$seed,
                                                    true_data = YYs_K4_300,
                                                    Ksat = 6),
                                       mc.cores = parallel::detectCores())
end=Sys.time()
end-start
save(gicsim_K4_T300,file="GIC_K4_T300_pers90.Rdata")
save(gicsim_sat_K4_T300,file="GICsat_K4_T300_pers90.Rdata")




# T=600
start=Sys.time()
gicsim_K4_T600 <- parallel::mclapply(1:nrow(hp),
                                     function(x) 
                                       sim_est_SJM(seed=hp[x,]$seed,
                                                   K=hp[x,]$K,
                                                   lambda=hp[x,]$lambda,
                                                   kappa=hp[x,]$kappa,
                                                   true_data=YYs_K4_600,
                                                   Ktrue=4,
                                                   N=Ns[2],
                                                   P=100,
                                                   m=2,
                                                   cors=corsK4,
                                                   pers=pers1,
                                                   K0=4,
                                                   alpha0=100),
                                     mc.cores = parallel::detectCores())
gicsim_sat_K4_T600<-parallel::mclapply(1:nrow(hpsat),
                                       function(x) 
                                         satmod_est(seed=hpsat[x,]$seed,
                                                    true_data = YYs_K4_600,
                                                    Ksat = 6),
                                       mc.cores = parallel::detectCores())
end=Sys.time()
end-start
save(gicsim_K4_T600,file="GIC_K4_T600_pers90.Rdata")
save(gicsim_sat_K4_T600,file="GICsat_K4_T600_pers90.Rdata")

# T=1000
start=Sys.time()
gicsim_K4_T1000 <- parallel::mclapply(1:nrow(hp),
                                      function(x) 
                                        sim_est_SJM(seed=hp[x,]$seed,
                                                    K=hp[x,]$K,
                                                    lambda=hp[x,]$lambda,
                                                    kappa=hp[x,]$kappa,
                                                    true_data=YYs_K4_1000,
                                                    Ktrue=4,
                                                    N=Ns[3],
                                                    P=100,
                                                    m=2,
                                                    cors=corsK4,
                                                    pers=pers1,
                                                    K0=4,
                                                    alpha0=100),
                                      mc.cores = parallel::detectCores())
gicsim_sat_K4_T1000<-parallel::mclapply(1:nrow(hpsat),
                                        function(x) 
                                          satmod_est(seed=hpsat[x,]$seed,
                                                     true_data = YYs_K4_1000,
                                                     Ksat = 6),
                                        mc.cores = parallel::detectCores())
end=Sys.time()
end-start
save(gicsim_K4_T1000,file="GIC_K4_T1000_pers90.Rdata")
save(gicsim_sat_K4_T1000,file="GICsat_K4_T1000_pers90.Rdata")

# pers 80 -----------------------------------------------------------------
pers1=.8

# Ktrue=2 -----------------------------------------------------------------
corsK2=c(.8,.4)
YYs_K2_300=lapply(seed,sim_data,Ktrue=2,N=Ns[1],P=100,cors=corsK2,pers=pers1,m=2)
YYs_K2_600=lapply(seed,sim_data,Ktrue=2,N=Ns[2],P=100,cors=corsK2,pers=pers1,m=2)
YYs_K2_1000=lapply(seed,sim_data,Ktrue=2,N=Ns[3],P=100,cors=corsK2,pers=pers1,m=2)

# T=300
start=Sys.time()
gicsim_K2_T300 <- parallel::mclapply(1:nrow(hp),
                                     function(x)
                                       sim_est_SJM(seed=hp[x,]$seed,
                                                   K=hp[x,]$K,
                                                   lambda=hp[x,]$lambda,
                                                   kappa=hp[x,]$kappa,
                                                   true_data=YYs_K2_300,
                                                   Ktrue=2,
                                                   N=Ns[1],
                                                   P=100,
                                                   m=2,
                                                   cors=corsK2,
                                                   pers=pers1,
                                                   K0=2,
                                                   alpha0=100),
                                     mc.cores = parallel::detectCores())
gicsim_sat_K2_T300<-parallel::mclapply(1:nrow(hpsat),
                                       function(x) 
                                         satmod_est(seed=hpsat[x,]$seed,
                                                    true_data = YYs_K2_300,
                                                    Ksat = 6),
                                       mc.cores = parallel::detectCores())
end=Sys.time()
end-start
save(gicsim_K2_T300,file="GIC_K2_T300_pers80.Rdata")
save(gicsim_sat_K2_T300,file="GICsat_K2_T300_pers80.Rdata")




# T=600
start=Sys.time()
gicsim_K2_T600 <- parallel::mclapply(1:nrow(hp),
                                     function(x)
                                       sim_est_SJM(seed=hp[x,]$seed,
                                                   K=hp[x,]$K,
                                                   lambda=hp[x,]$lambda,
                                                   kappa=hp[x,]$kappa,
                                                   true_data=YYs_K2_600,
                                                   Ktrue=2,
                                                   N=Ns[2],
                                                   P=100,
                                                   m=2,
                                                   cors=corsK2,
                                                   pers=pers1,
                                                   K0=2,
                                                   alpha0=100),
                                     mc.cores = parallel::detectCores())
gicsim_sat_K2_T600<-parallel::mclapply(1:nrow(hpsat),
                                       function(x) 
                                         satmod_est(seed=hpsat[x,]$seed,
                                                    true_data = YYs_K2_600,
                                                    Ksat = 6),
                                       mc.cores = parallel::detectCores())
end=Sys.time()
end-start
save(gicsim_K2_T600,file="GIC_K2_T600_pers80.Rdata")
save(gicsim_sat_K2_T600,file="GICsat_K2_T600_pers80.Rdata")

# T=1000
start=Sys.time()
gicsim_K2_T1000 <- parallel::mclapply(1:nrow(hp),
                                      function(x) 
                                        sim_est_SJM(seed=hp[x,]$seed,
                                                    K=hp[x,]$K,
                                                    lambda=hp[x,]$lambda,
                                                    kappa=hp[x,]$kappa,
                                                    true_data=YYs_K2_1000,
                                                    Ktrue=2,
                                                    N=Ns[3],
                                                    P=100,
                                                    m=2,
                                                    cors=corsK2,
                                                    pers=pers1,
                                                    K0=2,
                                                    alpha0=100),
                                      mc.cores = parallel::detectCores())
gicsim_sat_K2_T1000<-parallel::mclapply(1:nrow(hpsat),
                                        function(x) 
                                          satmod_est(seed=hpsat[x,]$seed,
                                                     true_data = YYs_K2_1000,
                                                     Ksat = 6),
                                        mc.cores = parallel::detectCores())
end=Sys.time()
end-start
save(gicsim_K2_T1000,file="GIC_K2_T1000_pers80.Rdata")
save(gicsim_sat_K2_T1000,file="GICsat_K2_T1000_pers80.Rdata")


# Ktrue=3 -----------------------------------------------------------------
corsK3=c(.8,.6,.3)
YYs_K3_300=lapply(seed,sim_data,Ktrue=3,N=Ns[1],P=100,cors=corsK3,pers=pers1,m=2)
YYs_K3_600=lapply(seed,sim_data,Ktrue=3,N=Ns[2],P=100,cors=corsK3,pers=pers1,m=2)
YYs_K3_1000=lapply(seed,sim_data,Ktrue=3,N=Ns[3],P=100,cors=corsK3,pers=pers1,m=2)

# T=300
start=Sys.time()
gicsim_K3_T300 <- parallel::mclapply(1:nrow(hp),
                                     function(x) 
                                       sim_est_SJM(seed=hp[x,]$seed,
                                                   K=hp[x,]$K,
                                                   lambda=hp[x,]$lambda,
                                                   kappa=hp[x,]$kappa,
                                                   true_data=YYs_K3_300,
                                                   Ktrue=3,
                                                   N=Ns[1],
                                                   P=100,
                                                   m=2,
                                                   cors=corsK3,
                                                   pers=pers1,
                                                   K0=3,
                                                   alpha0=100),
                                     mc.cores = parallel::detectCores())
gicsim_sat_K3_T300<-parallel::mclapply(1:nrow(hpsat),
                                       function(x) 
                                         satmod_est(seed=hpsat[x,]$seed,
                                                    true_data = YYs_K3_300,
                                                    Ksat = 6),
                                       mc.cores = parallel::detectCores())
end=Sys.time()
end-start
save(gicsim_K3_T300,file="GIC_K3_T300_pers80.Rdata")
save(gicsim_sat_K3_T300,file="GICsat_K3_T300_pers80.Rdata")




# T=600
start=Sys.time()
gicsim_K3_T600 <- parallel::mclapply(1:nrow(hp),
                                     function(x) 
                                       sim_est_SJM(seed=hp[x,]$seed,
                                                   K=hp[x,]$K,
                                                   lambda=hp[x,]$lambda,
                                                   kappa=hp[x,]$kappa,
                                                   true_data=YYs_K3_600,
                                                   Ktrue=3,
                                                   N=Ns[2],
                                                   P=100,
                                                   m=2,
                                                   cors=corsK3,
                                                   pers=pers1,
                                                   K0=3,
                                                   alpha0=100),
                                     mc.cores = parallel::detectCores())
gicsim_sat_K3_T600<-parallel::mclapply(1:nrow(hpsat),
                                       function(x) 
                                         satmod_est(seed=hpsat[x,]$seed,
                                                    true_data = YYs_K3_600,
                                                    Ksat = 6),
                                       mc.cores = parallel::detectCores())
end=Sys.time()
end-start
save(gicsim_K3_T600,file="GIC_K3_T600_pers80.Rdata")
save(gicsim_sat_K3_T600,file="GICsat_K3_T600_pers80.Rdata")

# T=1000
start=Sys.time()
gicsim_K3_T1000 <- parallel::mclapply(1:nrow(hp),
                                      function(x) 
                                        sim_est_SJM(seed=hp[x,]$seed,
                                                    K=hp[x,]$K,
                                                    lambda=hp[x,]$lambda,
                                                    kappa=hp[x,]$kappa,
                                                    true_data=YYs_K3_1000,
                                                    Ktrue=3,
                                                    N=Ns[3],
                                                    P=100,
                                                    m=2,
                                                    cors=corsK3,
                                                    pers=pers1,
                                                    K0=3,
                                                    alpha0=100),
                                      mc.cores = parallel::detectCores())
gicsim_sat_K3_T1000<-parallel::mclapply(1:nrow(hpsat),
                                        function(x) 
                                          satmod_est(seed=hpsat[x,]$seed,
                                                     true_data = YYs_K3_1000,
                                                     Ksat = 6),
                                        mc.cores = parallel::detectCores())
end=Sys.time()
end-start
save(gicsim_K3_T1000,file="GIC_K3_T1000_pers80.Rdata")
save(gicsim_sat_K3_T1000,file="GICsat_K3_T1000_pers80.Rdata")


# Ktrue=4 -----------------------------------------------------------------
corsK4=c(.8,.6,.3,0)
YYs_K4_300=lapply(seed,sim_data,Ktrue=4,N=Ns[1],P=100,cors=corsK4,pers=pers1,m=2)
YYs_K4_600=lapply(seed,sim_data,Ktrue=4,N=Ns[2],P=100,cors=corsK4,pers=pers1,m=2)
YYs_K4_1000=lapply(seed,sim_data,Ktrue=4,N=Ns[3],P=100,cors=corsK4,pers=pers1,m=2)

# T=300
start=Sys.time()
gicsim_K4_T300 <- parallel::mclapply(1:nrow(hp),
                                     function(x) 
                                       sim_est_SJM(seed=hp[x,]$seed,
                                                   K=hp[x,]$K,
                                                   lambda=hp[x,]$lambda,
                                                   kappa=hp[x,]$kappa,
                                                   true_data=YYs_K4_300,
                                                   Ktrue=4,
                                                   N=Ns[1],
                                                   P=100,
                                                   m=2,
                                                   cors=corsK4,
                                                   pers=pers1,
                                                   K0=4,
                                                   alpha0=100),
                                     mc.cores = parallel::detectCores())
gicsim_sat_K4_T300<-parallel::mclapply(1:nrow(hpsat),
                                       function(x) 
                                         satmod_est(seed=hpsat[x,]$seed,
                                                    true_data = YYs_K4_300,
                                                    Ksat = 6),
                                       mc.cores = parallel::detectCores())
end=Sys.time()
end-start
save(gicsim_K4_T300,file="GIC_K4_T300_pers80.Rdata")
save(gicsim_sat_K4_T300,file="GICsat_K4_T300_pers80.Rdata")




# T=600
start=Sys.time()
gicsim_K4_T600 <- parallel::mclapply(1:nrow(hp),
                                     function(x) 
                                       sim_est_SJM(seed=hp[x,]$seed,
                                                   K=hp[x,]$K,
                                                   lambda=hp[x,]$lambda,
                                                   kappa=hp[x,]$kappa,
                                                   true_data=YYs_K4_600,
                                                   Ktrue=4,
                                                   N=Ns[2],
                                                   P=100,
                                                   m=2,
                                                   cors=corsK4,
                                                   pers=pers1,
                                                   K0=4,
                                                   alpha0=100),
                                     mc.cores = parallel::detectCores())
gicsim_sat_K4_T600<-parallel::mclapply(1:nrow(hpsat),
                                       function(x) 
                                         satmod_est(seed=hpsat[x,]$seed,
                                                    true_data = YYs_K4_600,
                                                    Ksat = 6),
                                       mc.cores = parallel::detectCores())
end=Sys.time()
end-start
save(gicsim_K4_T600,file="GIC_K4_T600_pers80.Rdata")
save(gicsim_sat_K4_T600,file="GICsat_K4_T600_pers80.Rdata")

# T=1000
start=Sys.time()
gicsim_K4_T1000 <- parallel::mclapply(1:nrow(hp),
                                      function(x) 
                                        sim_est_SJM(seed=hp[x,]$seed,
                                                    K=hp[x,]$K,
                                                    lambda=hp[x,]$lambda,
                                                    kappa=hp[x,]$kappa,
                                                    true_data=YYs_K4_1000,
                                                    Ktrue=4,
                                                    N=Ns[3],
                                                    P=100,
                                                    m=2,
                                                    cors=corsK4,
                                                    pers=pers1,
                                                    K0=4,
                                                    alpha0=100),
                                      mc.cores = parallel::detectCores())
gicsim_sat_K4_T1000<-parallel::mclapply(1:nrow(hpsat),
                                        function(x) 
                                          satmod_est(seed=hpsat[x,]$seed,
                                                     true_data = YYs_K4_1000,
                                                     Ksat = 6),
                                        mc.cores = parallel::detectCores())
end=Sys.time()
end-start
save(gicsim_K4_T1000,file="GIC_K4_T1000_pers80.Rdata")
save(gicsim_sat_K4_T1000,file="GICsat_K4_T1000_pers80.Rdata")


# pers 70 -----------------------------------------------------------------

pers1=.7
# Ktrue=2 -----------------------------------------------------------------
corsK2=c(.8,.4)
YYs_K2_1000=lapply(seed,sim_data,Ktrue=2,N=Ns[3],P=100,cors=corsK2,pers=pers1,m=2)

# T=1000
start=Sys.time()
gicsim_K2_T1000 <- parallel::mclapply(1:nrow(hp),
                                      function(x) 
                                        sim_est_SJM(seed=hp[x,]$seed,
                                                    K=hp[x,]$K,
                                                    lambda=hp[x,]$lambda,
                                                    kappa=hp[x,]$kappa,
                                                    true_data=YYs_K2_1000,
                                                    Ktrue=2,
                                                    N=Ns[3],
                                                    P=100,
                                                    m=2,
                                                    cors=corsK2,
                                                    pers=pers1,
                                                    K0=2,
                                                    alpha0=100),
                                      mc.cores = parallel::detectCores())
gicsim_sat_K2_T1000<-parallel::mclapply(1:nrow(hpsat),
                                        function(x) 
                                          satmod_est(seed=hpsat[x,]$seed,
                                                     true_data = YYs_K2_1000,
                                                     Ksat = 6),
                                        mc.cores = parallel::detectCores())
end=Sys.time()
end-start
save(gicsim_K2_T1000,file="GIC_K2_T1000_pers70.Rdata")
save(gicsim_sat_K2_T1000,file="GICsat_K2_T1000_pers70.Rdata")


# Ktrue=3 -----------------------------------------------------------------
corsK3=c(.8,.6,.3)
YYs_K3_1000=lapply(seed,sim_data,Ktrue=3,N=Ns[3],P=100,cors=corsK3,pers=pers1,m=2)

# T=1000
start=Sys.time()
gicsim_K3_T1000 <- parallel::mclapply(1:nrow(hp),
                                      function(x) 
                                        sim_est_SJM(seed=hp[x,]$seed,
                                                    K=hp[x,]$K,
                                                    lambda=hp[x,]$lambda,
                                                    kappa=hp[x,]$kappa,
                                                    true_data=YYs_K3_1000,
                                                    Ktrue=3,
                                                    N=Ns[3],
                                                    P=100,
                                                    m=2,
                                                    cors=corsK3,
                                                    pers=pers1,
                                                    K0=3,
                                                    alpha0=100),
                                      mc.cores = parallel::detectCores())
gicsim_sat_K3_T1000<-parallel::mclapply(1:nrow(hpsat),
                                        function(x) 
                                          satmod_est(seed=hpsat[x,]$seed,
                                                     true_data = YYs_K3_1000,
                                                     Ksat = 6),
                                        mc.cores = parallel::detectCores())
end=Sys.time()
end-start
save(gicsim_K3_T1000,file="GIC_K3_T1000_pers70.Rdata")
save(gicsim_sat_K3_T1000,file="GICsat_K3_T1000_pers70.Rdata")


# Ktrue=4 -----------------------------------------------------------------
corsK4=c(.8,.6,.3,0)
YYs_K4_1000=lapply(seed,sim_data,Ktrue=4,N=Ns[3],P=100,cors=corsK4,pers=pers1,m=2)

# T=1000
start=Sys.time()
gicsim_K4_T1000 <- parallel::mclapply(1:nrow(hp),
                                      function(x) 
                                        sim_est_SJM(seed=hp[x,]$seed,
                                                    K=hp[x,]$K,
                                                    lambda=hp[x,]$lambda,
                                                    kappa=hp[x,]$kappa,
                                                    true_data=YYs_K4_1000,
                                                    Ktrue=4,
                                                    N=Ns[3],
                                                    P=100,
                                                    m=2,
                                                    cors=corsK4,
                                                    pers=pers1,
                                                    K0=4,
                                                    alpha0=100),
                                      mc.cores = parallel::detectCores())
gicsim_sat_K4_T1000<-parallel::mclapply(1:nrow(hpsat),
                                        function(x) 
                                          satmod_est(seed=hpsat[x,]$seed,
                                                     true_data = YYs_K4_1000,
                                                     Ksat = 6),
                                        mc.cores = parallel::detectCores())
end=Sys.time()
end-start
save(gicsim_K4_T1000,file="GIC_K4_T1000_pers70.Rdata")
save(gicsim_sat_K4_T1000,file="GICsat_K4_T1000_pers70.Rdata")


# output ------------------------------------------------------------------

#GIC=function(simest,satmod,Ksat=6){
#GIC=function(simest,satmod,Ksat=6,alpha0,K0,pen0){
GIC=function(simest,satmod,Ksat=6,alpha0,K0){
  
  library(dplyr)
  library(tidyverse)
  library(tidyr)
  
  K=unlist(lapply(simest,function(x)x$K))
  kappa=unlist(lapply(simest,function(x)x$kappa))
  lambda=unlist(lapply(simest,function(x)x$lambda))
  #YY=unlist(lapply(simest,function(x)x$K))
  PP=unlist(lapply(simest,function(x)x$P))
  PP=PP*(1+unlist(lapply(simest,function(x)x$m)))
  N=unlist(lapply(simest,function(x)x$N))
  seed=unlist(lapply(simest,function(x)x$seed))
  ARI.weights=unlist(lapply(simest,function(x)x$ARI_weights))
  ARI.states=unlist(lapply(simest,function(x)x$ARI_states))
  pen=unlist(lapply(simest,function(x)x$pen))
  
  # if(drop_pterm){
  #   CKp=-log(K)-log(2*pi)*1/2
  #   CKp_sat=-log(Ksat)-log(2*pi)*1/2
  # }
  # else{
  CKp=-log(K)-log(2*pi)*kappa/2
  CKp_sat=-log(Ksat)-log(2*pi)*sqrt(PP)/2
  #}
  
  anFTIC=log(log(N))*log(PP)
  anAIC=rep(2,length(N))
  anBIC=log(N)
  
  # TotalPenalty=unlist(lapply(simest,function(x)x$TotalPenalty))
  alphak=unlist(lapply(simest,function(x)x$alphak))
  pers=unlist(lapply(simest,function(x)x$pers))
  #pen0=(1-pers)*N
  pen0=(1-pers)*N*(K0-1)
  
  #TotalPenalty=(alpha0+pen0)*(K-K0)+K0*(alphak-alpha0+pen-pen0)
  # TotalPenalty=(alpha0+pen0)*K+K0*(alphak-alpha0+pen-pen0) #no big changes fortunately
  
  p_hat=pen/(N*(K-1))
  p0_hat=pen0/(N*(K0-1))
  
  #TotalPenalty=alphak*K*(1-(1-p_hat)^N) #first Stas version: NOT WORKING
  TotalPenalty=alphak*(1+(K-1)*(1-(1-p_hat)^N)) #second Stas version: NOT WORKING on real data
  
  #  TotalPenalty=alphak*K+(1+(K-1)*(1-(1-p_hat)^N)) #third Stas version (keep jump penalty sep)
  
  # TotalPenalty=alpha0*(1+(K0-1)*(1-(1-p0_hat)^N))+
  #   (1+(K0-1)*(1-(1-p0_hat)^N))*(alphak-alpha0)+
  #   alpha0*(1-(1-p0_hat)^N)*(K-K0)+
  #   (N*(1-p0_hat)^(N-1)*alpha0*(K0-1))*(p_hat-p0_hat)
  
  Klk=length(unique(K))*length(unique(kappa))*length(unique(lambda))
  setup=rep(1:Klk,
            each=length(unique(seed)))
  
  Ln=unlist(lapply(simest,function(x)x$Ln1))
  Lnsat=unlist(lapply(satmod,function(x)x$Lnsat1))
  Lnsat=(rep(Lnsat,Klk))
  Ln_diff=Lnsat-Ln
  CKp_diff=CKp_sat-CKp
  # if(byN){
  #   FTIC=2*CKp_diff/N+(Ln_diff+anFTIC*TotalPenalty)/N
  #   BIC=2*CKp_diff/N+(Ln_diff+anBIC*TotalPenalty)/N
  #   AIC=2*CKp_diff/N+(Ln_diff+anAIC*TotalPenalty)/N
  # }
  # else{
  FTIC=2*CKp_diff+(Ln_diff+anFTIC*TotalPenalty)/N
  BIC=2*CKp_diff+(Ln_diff+anBIC*TotalPenalty)/N
  AIC=2*CKp_diff+(Ln_diff+anAIC*TotalPenalty)/N
  #}
  
  corr_sel_feat=unlist(lapply(simest,function(x)x$corr_sel_feat))
  wrong_sel_feat=unlist(lapply(simest,function(x)x$wrong_sel_feat))
  est_weights=data.frame(matrix(unlist(lapply(simest,function(x)x$norm_est_weights)),
                                ncol = PP[1],byrow=T))
  
  res=data.frame(seed,FTIC,BIC,AIC,
                 Ln_diff,CKp_diff,pen,
                 lambda,kappa,K,
                 ARI.states,ARI.weights,
                 setup,
                 corr_sel_feat,
                 wrong_sel_feat)
  
  res_est_weights=data.frame(setup,est_weights)
  
  res=res%>%group_by(setup)%>%
    summarise(avFTIC=mean(FTIC),
              avAIC=mean(AIC),
              avBIC=mean(BIC),
              sdFTIC=sd(FTIC),
              sdAIC=sd(AIC),
              sdBIC=sd(BIC),
              avARI_states=mean(ARI.states),
              avARI_weights=mean(ARI.weights),
              avLn_diff=mean(Ln_diff),
              avCKp_diff=mean(CKp_diff),
              avPen=mean(pen),
              lambda=mean(lambda),
              kappa=mean(kappa),
              K=mean(K),
              corr_sel_feat=mean(corr_sel_feat),
              wrong_sel_feat=mean(wrong_sel_feat))
  
  res_est_weights=res_est_weights%>%group_by(setup)%>%
    summarise_all(mean)
  
  res_av=data.frame(res,res_est_weights[,-1])
  
  return(res_av)
}

tab_out=function(res){
  FTIC=round(head(res[order(res$avFTIC),c(2,1,13:15,8:9,16:17)],1),2)
  colnames(FTIC)[1]="IC"
  AIC=round(head(res[order(res$avAIC),c(3,1,13:15,8:9,16:17)],1),2)
  colnames(AIC)[1]="IC"
  BIC=round(head(res[order(res$avBIC),c(4,1,13:15,8:9,16:17)],1),2)
  colnames(BIC)[1]="IC"
  out=rbind(FTIC,AIC,BIC)
  print(xtable(out[,-2], type = "latex"))
  return(out$setup)
}

varPlot=function(datavarofGIC,legend="none"){
  datavarofGIC$K=as.factor(datavarofGIC$K)
  by.kappaK=datavarofGIC%>%group_by(kappa,K)
  dataplotB2=by.kappaK%>%summarise(
    FTIC=min(avFTIC),
    # minFTIC=min(minFTIC),
    # maxFTIC=min(maxFTIC),
    sdFTIC=min(sdFTIC),
    AIC=min(avAIC),
    # minAIC=min(minAIC),
    # maxAIC=min(maxAIC),
    sdAIC=min(sdAIC),
    BIC=min(avBIC),
    # minBIC=min(minBIC),
    # maxBIC=min(maxBIC),
    sdBIC=min(sdBIC)
  ) 
  
  B2=ggplot(dataplotB2, aes(x = kappa, y = FTIC)) + 
    geom_line(aes(color = K),size=.8)+
    #geom_ribbon(aes(ymin=minGIC, ymax=maxGIC,fill = K), alpha=.18, linetype=0)+
    geom_ribbon(aes(ymin=FTIC-2*sdFTIC, ymax=FTIC+2*sdFTIC,fill = K), alpha=.3, linetype=0)+
    geom_point(aes(color = K,shape=K),size=2)+
    scale_x_continuous(name=expression(kappa)) +
    scale_y_continuous(name=" ")+
    #ggtitle(expression(FTIC~", "~K[true]~"= 2"))+
    theme(plot.title = element_text(face = "bold", size = 12),
          axis.ticks = element_line(colour = "grey70", size = 0.2),
          panel.grid.major = element_line(colour = "grey70", size = 0.2),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = 'white', color = 'white'),
          
          legend.position = "none")
  
  B2AIC=ggplot(dataplotB2, aes(x = kappa, y = AIC)) + 
    geom_line(aes(color = K),size=.8)+
    # geom_ribbon(aes(ymin=minAIC, ymax=maxAIC,fill = K), alpha=.18, linetype=0)+
    geom_ribbon(aes(ymin=AIC-2*sdAIC, ymax=AIC+2*sdAIC,fill = K), alpha=.3, linetype=0)+
    geom_point(aes(color = as.factor(K),shape=K),size=2)+
    scale_x_continuous(name=expression(kappa)) +
    scale_y_continuous(name=" ")+
    #ggtitle(expression(AIC~", "~K[true]~"= 2"))+
    theme(plot.title = element_text(face = "bold", size = 12),
          legend.position = "none",
          panel.background = element_rect(fill = 'white', color = 'white'),
          axis.ticks = element_line(colour = "grey70", size = 0.2),
          panel.grid.major = element_line(colour = "grey70", size = 0.2),
          panel.grid.minor = element_blank())
  
  B2BIC=ggplot(dataplotB2, aes(x = kappa, y = BIC)) + 
    geom_line(aes(color = K),size=.8)+
    #geom_ribbon(aes(ymin=minBIC, ymax=maxBIC,fill = K), alpha=.18, linetype=0)+
    geom_ribbon(aes(ymin=BIC-2*sdBIC, ymax=BIC+2*sdBIC,fill = K), alpha=.3, linetype=0)+
    geom_point(aes(color = as.factor(K),shape=K),size=2)+
    scale_x_continuous(name=expression(kappa)) +
    scale_y_continuous(name=" ")+
    #ggtitle(expression(BIC~", "~K[true]~"= 2"))+
    theme(plot.title = element_text(face = "bold", size = 12),
          legend.position = "none",
          panel.background = element_rect(fill = 'white', color = 'white'),
          axis.ticks = element_line(colour = "grey70", size = 0.2),
          panel.grid.major = element_line(colour = "grey70", size = 0.2),
          panel.grid.minor = element_blank())
  
  by.lambdaK=datavarofGIC%>%group_by(lambda,K)
  
  dataplotA2=by.lambdaK%>%summarise(
    FTIC=min(avFTIC),
    # minGIC=min(minGIC),
    # maxGIC=min(maxGIC),
    sdFTIC=min(sdFTIC),
    AIC=min(avAIC),
    # minAIC=min(minAIC),
    # maxAIC=min(maxAIC),
    sdAIC=min(sdAIC),
    BIC=min(avBIC),
    # minBIC=min(minBIC),
    # maxBIC=min(maxBIC),
    sdBIC=min(sdBIC)
  ) 
  
  
  
  A2=ggplot(dataplotA2, aes(x = lambda, y = FTIC)) + 
    geom_line(aes(color = K),size=.8)+
    #geom_ribbon(aes(ymin=minGIC, ymax=maxGIC,fill = K), alpha=.18, linetype=0)+
    geom_ribbon(aes(ymin=FTIC-2*sdFTIC, ymax=FTIC+2*sdFTIC,fill = K), alpha=.3, linetype=0)+
    geom_point(aes(color = as.factor(K),shape=K),size=2)+
    scale_x_continuous(name=expression(lambda)) +
    scale_y_continuous(name=" ")+
    # ggtitle(expression(FTIC~", "~K[true]~"= 2"))+
    ggtitle(" ")+
    theme(plot.title = element_text(face = "bold", size = 12),
          legend.position = "none",
          panel.background = element_rect(fill = 'white', color = 'white'),
          #legend.position = c(0.8, 0.2),
          axis.ticks = element_line(colour = "grey70", size = 0.2),
          panel.grid.major = element_line(colour = "grey70", size = 0.2),
          panel.grid.minor = element_blank()) 
  
  
  A2AIC=ggplot(dataplotA2, aes(x = lambda, y = AIC)) + 
    geom_line(aes(color = K),size=.8)+
    # geom_ribbon(aes(ymin=minAIC, ymax=maxAIC,fill = K), alpha=.18, linetype=0)+
    geom_ribbon(aes(ymin=AIC-2*sdAIC, ymax=AIC+2*sdAIC,fill = K), alpha=.3, linetype=0)+
    geom_point(aes(color = as.factor(K),shape=K),size=2)+
    scale_x_continuous(name=expression(lambda)) +
    scale_y_continuous(name=" ")+
    #ggtitle(expression(AIC~", "~K[true]~"= 2"))+
    ggtitle(" ")+
    theme(plot.title = element_text(face = "bold", size = 12),
          legend.position = "none",
          panel.background = element_rect(fill = 'white', color = 'white'),
          #legend.position = c(0.8, 0.2),
          axis.ticks = element_line(colour = "grey70", size = 0.2),
          panel.grid.major = element_line(colour = "grey70", size = 0.2),
          panel.grid.minor = element_blank())
  
  A2BIC=ggplot(dataplotA2, aes(x = lambda, y = BIC)) + 
    geom_line(aes(color = K),size=.8)+
    #geom_ribbon(aes(ymin=minBIC, ymax=maxBIC,fill = K), alpha=.18, linetype=0)+
    geom_ribbon(aes(ymin=BIC-2*sdBIC, ymax=BIC+2*sdBIC,fill = K), alpha=.3, linetype=0)+
    geom_point(aes(color = as.factor(K),shape=K),size=2)+
    scale_x_continuous(name=expression(lambda)) +
    scale_y_continuous(name=" ")+
    #ggtitle(expression(BIC~", "~K[true]~"= 2"))+
    ggtitle(" ")+
    theme(plot.title = element_text(face = "bold", size = 12),
          legend.position = "none",
          panel.background = element_rect(fill = 'white', color = 'white'),
          #legend.position = c(0.8, 0.2),
          axis.ticks = element_line(colour = "grey70", size = 0.2),
          panel.grid.major = element_line(colour = "grey70", size = 0.2),
          panel.grid.minor = element_blank())
  
  
  
  
  library(ggpubr)
  out=ggarrange(
    A2,A2AIC,A2BIC,
    B2,B2AIC,B2BIC,
    labels=c("FTIC","AIC","BIC",
             " "," "," ",
             " "," "," "),
    label.x = .45,
    ncol=3,nrow=2,
    common.legend = T,
    legend=legend)
  return(out)
}


# pers 90 -----------------------------------------------------------------

### T=300 

library(xtable)
# pen0=(1-pers)*N
#load("allRes_pers90.Rdata)
resK2T300=GIC(gicsim_K2_T300,gicsim_sat_K2_T300,alpha0 = 100,K0=2)
#tab_out(resK2T300)
#resK2T300=GIC(gicsim_K2_T300,gicsim_sat_K2_T300)
resK2T300_constr=resK2T300[-which(resK2T300$avPen>120),]
tab_out(resK2T300_constr)

resK3T300=GIC(gicsim_K3_T300,gicsim_sat_K3_T300,alpha0 = 100,K0=3)
#resK3T300=GIC(gicsim_K3_T300,gicsim_sat_K3_T300)
resK3T300_constr=resK3T300[-which(resK3T300$avPen>120),]
tab_out(resK3T300_constr)

resK4T300=GIC(gicsim_K4_T300,gicsim_sat_K4_T300)
resK4T300_constr=resK4T300[-which(resK4T300$avPen>120),]
tab_out(resK4T300_constr)

### T=600 

resK2T600=GIC(gicsim_K2_T600,gicsim_sat_K2_T600,alpha0 = 100,K0=2)
#resK2T600=GIC(gicsim_K2_T600,gicsim_sat_K2_T600)
resK2T600_constr=resK2T600[-which(resK2T600$avPen>240),]
tab_out(resK2T600_constr)

resK3T600=GIC(gicsim_K3_T600,gicsim_sat_K3_T600,alpha0 = 100,K0=3)
#resK3T600=GIC(gicsim_K3_T600,gicsim_sat_K3_T600)
resK3T600_constr=resK3T600[-which(resK3T600$avPen>240),]
tab_out(resK3T600_constr)

resK4T600=GIC(gicsim_K4_T600,gicsim_sat_K4_T600)
resK4T600_constr=resK4T600[-which(resK4T600$avPen>240),]
tab_out(resK4T600_constr)

### T=1000 

resK2T1000=GIC(gicsim_K2_T1000,gicsim_sat_K2_T1000,alpha0 = 100,K0=2)
tab_out(resK2T1000)
# resK2T1000=GIC(gicsim_K2_T1000,gicsim_sat_K2_T1000)
resK2T1000_constr=resK2T1000[-which(resK2T1000$avPen>400),]
tab_out(resK2T1000_constr)

resK3T1000=GIC(gicsim_K3_T1000,gicsim_sat_K3_T1000,alpha0 = 100,K0=3)
#resK3T1000=GIC(gicsim_K3_T1000,gicsim_sat_K3_T1000)
resK3T1000_constr=resK3T1000[-which(resK3T1000$avPen>400),]
tab_out(resK3T1000_constr)

resK4T1000=GIC(gicsim_K4_T1000,gicsim_sat_K4_T1000,alpha0 = 100,K0=4)
#resK4T1000=GIC(gicsim_K4_T1000,gicsim_sat_K4_T1000)
resK4T1000_constr=resK4T1000[-which(resK4T1000$avPen>400),]
tab_out(resK4T1000_constr)

varPlot(resK2T300_constr)
varPlot(resK3T300_constr)
varPlot(resK4T300_constr,legend="bottom")

varPlot(resK2T600_constr)
varPlot(resK3T600_constr)
varPlot(resK4T600_constr)

varPlot(resK2T1000_constr)
varPlot(resK3T1000_constr)
varPlot(resK4T1000_constr,legend="bottom")

# pers 80 -----------------------------------------------------------------
library(xtable)

plot_weights=function(res){
  stp=tab_out(res)
  plot(x=1:300,
       y=res[which(res$setup==stp[1]),-(1:17)],
       ylab="Weight",xlab="Feature",pch=19)
}

### T=300 

#load("allRes_pers80.Rdata)
resK2T300=GIC(gicsim_K2_T300,gicsim_sat_K2_T300,alpha0=100,K0=2)
resK2T300_constr=resK2T300[-which(resK2T300$avPen>120),]
#stp=tab_out(resK2T300_constr)
tab_out(resK2T300)
tab_out(resK2T300_constr)
plot_weights(resK2T300_constr)

resK3T300=GIC(gicsim_K3_T300,gicsim_sat_K3_T300,alpha0=100,K0=3)
resK3T300_constr=resK3T300[-which(resK3T300$avPen>120),]
#stp=tab_out(resK3T300_constr)
plot_weights(resK3T300_constr)
tab_out(resK3T300)
tab_out(resK3T300_constr)

resK4T300=GIC(gicsim_K4_T300,gicsim_sat_K4_T300,alpha0=100,K0=4)
resK4T300_constr=resK4T300[-which(resK4T300$avPen>120),]
#stp=tab_out(resK4T300_constr)
plot_weights(resK4T300_constr)
tab_out(resK4T300)
tab_out(resK4T300_constr)

### T=600 

resK2T600=GIC(gicsim_K2_T600,gicsim_sat_K2_T600,alpha0=100,K0=2)
resK2T600_constr=resK2T600[-which(resK2T600$avPen>240),]
#stp=tab_out(resK2T600_constr)
plot_weights(resK2T600_constr)
tab_out(resK2T600)
tab_out(resK2T600_constr)

resK3T600=GIC(gicsim_K3_T600,gicsim_sat_K3_T600,alpha0=100,K0=3)
resK3T600_constr=resK3T600[-which(resK3T600$avPen>240),]
tab_out(resK3T600_constr)
plot_weights(resK3T600_constr)
tab_out(resK3T600)
tab_out(resK3T600_constr)

resK4T600=GIC(gicsim_K4_T600,gicsim_sat_K4_T600,alpha0=100,K0=4)
resK4T600_constr=resK4T600[-which(resK4T600$avPen>240),]
tab_out(resK4T600_constr)
plot_weights(resK4T600_constr)
tab_out(resK4T600)
tab_out(resK4T600_constr)

### T=1000 

resK2T1000=GIC(gicsim_K2_T1000,gicsim_sat_K2_T1000,alpha0=100,K0=2)
resK2T1000_constr=resK2T1000[-which(resK2T1000$avPen>400),]
tab_out(resK2T1000_constr)
plot_weights(resK2T1000_constr)
tab_out(resK2T1000)
tab_out(resK2T1000_constr)


resK3T1000=GIC(gicsim_K3_T1000,gicsim_sat_K3_T1000,alpha0=100,K0=3)
resK3T1000_constr=resK3T1000[-which(resK3T1000$avPen>400),]
tab_out(resK3T1000_constr)
plot_weights(resK3T1000_constr)
tab_out(resK3T1000)
tab_out(resK3T1000_constr)

resK4T1000=GIC(gicsim_K4_T1000,gicsim_sat_K4_T1000,alpha0=100,K0=4)
resK4T1000_constr=resK4T1000[-which(resK4T1000$avPen>400),]
tab_out(resK4T1000_constr)
plot_weights(resK4T1000_constr)

png(filename = "k2t300pers80.png",width=800,height=400)
varPlot(resK2T300_constr)
dev.off()

png(filename = "k3t300pers80.png",width=800,height=400)
varPlot(resK3T300_constr)
dev.off()

png(filename = "k4t300pers80.png",width=800,height=400)
varPlot(resK4T300_constr,legend="bottom")
dev.off()


png(filename = "k2t600pers80.png",width=800,height=400)
varPlot(resK2T600_constr)
dev.off()

png(filename = "k3t600pers80.png",width=800,height=400)
varPlot(resK3T600_constr)
dev.off()

png(filename = "k4t600pers80.png",width=800,height=400)
varPlot(resK4T600_constr,legend="bottom")
dev.off()

png(filename = "k2t1000pers80.png",width=800,height=400)
varPlot(resK2T1000_constr)
dev.off()

png(filename = "k3t1000pers80.png",width=800,height=400)
varPlot(resK3T1000_constr)
dev.off()

png(filename = "k4t1000pers80.png",width=800,height=400)
varPlot(resK4T1000_constr,legend="bottom")
dev.off()



# pers 70 -----------------------------------------------------------------

#load("allRes_pers70.Rdata)

### T=1000 

resK2T1000=GIC(gicsim_K2_T1000,gicsim_sat_K2_T1000,alpha0=100,K0=2)
resK2T1000_constr=resK2T1000[-which(resK2T1000$avPen>400),]
tab_out(resK2T1000_constr)

resK3T1000=GIC(gicsim_K3_T1000,gicsim_sat_K3_T1000,alpha0=100,K0=3)
resK3T1000_constr=resK3T1000[-which(resK3T1000$avPen>400),]
tab_out(resK3T1000_constr)

resK4T1000=GIC(gicsim_K4_T1000,gicsim_sat_K4_T1000,alpha0=100,K0=4)
resK4T1000_constr=resK4T1000[-which(resK4T1000$avPen>400),]
tab_out(resK4T1000_constr)


varPlot(resK2T1000_constr)
varPlot(resK3T1000_constr)
varPlot(resK4T1000_constr,legend="bottom")
