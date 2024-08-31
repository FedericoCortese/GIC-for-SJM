source("Utils.R")

#Simulate data
pers=.99
Ns=600
seed=123
corsK2=c(.8,.4)
P=100

Y=sim_data(seed,Ktrue=2,Ns,P,corsK2,pers,m=1)

true_stats=Y[[1]]
Y=Y[[2]]

# Estimate SJM with varying lambda and kappa
Yest_0=SJM_est(Y,K=3,kappa=2,lambda=0)
Yest_1=SJM_est(Y,K=2,kappa=7,lambda=10)
Yest_2=SJM_est(Y,K=2,kappa=16,lambda=50)

# Compare FTIC
Yest_0$FTIC
Yest_1$FTIC
Yest_2$FTIC

# Best model according to FTIC: Yest_1

# Compare in terms of ARI computed between true and estimated states
pdfCluster::adj.rand.index(Yest_0$est_states,true_stats)
pdfCluster::adj.rand.index(Yest_1$est_states,true_stats)
pdfCluster::adj.rand.index(Yest_2$est_states,true_stats)

# Best model according to ARI: Yest_1 and Yest_2

# Selected features (the first P are true features)
which(Yest_0$est_weights!=0)
which(Yest_1$est_weights!=0)
which(Yest_2$est_weights!=0)

# Best model according to correctly selected features: Yest_1
