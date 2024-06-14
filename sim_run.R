# pers 90 -----------------------------------------------------------------
pers1=.99
Ns=c(300,600,1000)
seed=123

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

