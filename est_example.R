source("Utils.R")
source("est.R")

#Simulate data
pers1=.99
Ns=600
seed=123
corsK2=c(.8,.4)
P=100

Y=sim_data(seed,Ktrue=2,Ns,P,corsK2,pers1,m=1)

true_stats=Y[[1]]
Y=Y[[2]]

# Estimate SJM with varying lambda and kappa
Yest_0=SJM_est(Y,K=2,kappa=2,lambda=0,
             Ksat=6,alpha0=100,K0=2,pers0=.95)
Yest_1=SJM_est(Y,K=2,kappa=7,lambda=10,
               Ksat=6,alpha0=100,K0=2,pers0=.95)
Yest_2=SJM_est(Y,K=2,kappa=16,lambda=50,
               Ksat=6,alpha0=100,K0=2,pers0=.95)

# Compare FTIC
Yest_0$FTIC
Yest_1$FTIC
Yest_2$FTIC

# Compare in terms of ARI
pdfCluster::adj.rand.index(Yest_0$est_states,true_stats)
pdfCluster::adj.rand.index(Yest_1$est_states,true_stats)
pdfCluster::adj.rand.index(Yest_2$est_states,true_stats)

# Selected features (the first P are true features)
which(Yest_0$est_weights!=0)
which(Yest_1$est_weights!=0)
which(Yest_2$est_weights!=0)
