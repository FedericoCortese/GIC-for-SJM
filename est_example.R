source("Utils.R")

#Simulate data
pers=.99
Ns=600
seed=123
corsK2=c(.8,.4)
P=100

# m*P is the number of false features 
Y=sim_data(seed,Ktrue=2,Ns,P,corsK2,pers,m=1)

true_stats=Y[[1]]
Y=Y[[2]]

dim(Y)
# Final datasets counts P=200 features, 100 of which are false features

# Estimate SJM with varying lambda and kappa
Yest_0=SJM_est(Y,K=3,kappa=2,lambda=0)
Yest_1=SJM_est(Y,K=2,kappa=7,lambda=50)
Yest_2=SJM_est(Y,K=2,kappa=16,lambda=1000)

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


# Automatic selection with GIC

res_grid=grid_search_SJM(
    Y,
    K_grid=NULL,
    kappa_grid=NULL,
    lambda_grid=NULL,
    parallel   = F,
    max_iter   = 10,
    n_init     = 10,
    tol        = 1e-8,
    Ksat       = 6,
    K0=NULL,
    pers0=NULL,
    alpha0=NULL
)

# CV

res_cv=cv_sparse_jump(
    Y=Y,
    true_states=true_stats,
    K_grid=2:4,
    kappa_grid=seq(1,sqrt(P),length.out=5),
    lambda_grid=c(0,.5,1,5,10,50,100,1000,10^4),
    n_folds = 5,
    parallel=F
)


res_cv[which.max(res_cv$ARI),]
best_K=res_cv$K[which.max(res_cv$ARI)]
best_kappa=res_cv$kappa[which.max(res_cv$ARI)]
best_lambda=res_cv$lambda[which.max(res_cv$ARI)]



