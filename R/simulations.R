

## Simuations to compare the performance against both full data and Ugander Methods

# rm(list = ls())
source("R/ardexp.R") # functions for the main
source("R/simulation_functions.R") # has the additional simulations pieces

## competing methods DM, HT, Oracle Regression, Partial Data Regression, Partial Data Regression (Estimated Clusters)


slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
print(Sys.getenv('SLURM_ARRAY_TASK_ID'))
if(slurm_arrayid == ""){
  id = 1
} else {
  # coerce the value to an integer
  id <- as.numeric(slurm_arrayid)
}

mutual.benefit
cluster.growth
sparse.graph

if(id %% 8 == 1){
  mutual.benefit = T
  cluster.growth = T
  sparse.graph = T

} else if(id %% 8 == 2) {
  mutual.benefit = F
  cluster.growth = T
  sparse.graph = T

} else if(id %% 8 == 3) {
  mutual.benefit = T
  cluster.growth = F
  sparse.graph = T

} else if(id %% 8 == 4) {
  mutual.benefit = F
  cluster.growth = F
  sparse.graph = T

} else if(id %% 8 == 5) {
  mutual.benefit = T
  cluster.growth = T
  sparse.graph = F

} else if(id %% 8 == 6) {
  mutual.benefit = F
  cluster.growth = T
  sparse.graph = F

} else if(id %% 8 == 7) {
  mutual.benefit = T
  cluster.growth = F
  sparse.graph = F

} else if(id %% 8 == 0) {
  mutual.benefit = F
  cluster.growth = F
  sparse.graph = F

}

block = ceiling(id/8)

# True model parameters
a = 1
b = 1
delta = 1

# simulation noise
sigma = 0.5# how much noise are we willing to permit here?


# TODO: make this a function of the arrays
mutual.benefit = T
if(mutual.benefit){
  gamma = 1
} else {
  gamma = -1
}


# Think of a second version of the simulation where there is a competative advantage to treatment
# but that wears off if everyone else around you has it

#Global average treatment effect
gate = gamma + delta

n.sims =  2#

n.seq <- c(100,316,1000,3162,10000) # log scale population growth.
J = length(n.seq) # number of sample sizes to look at

# treatment probability
p.treat = 0.5

#TODO: make this a function of the array id
cluster.growth = F
if(cluster.growth){
  K.seq = rep(10,J)
} else {
  K.seq = round(n.seq/10) # linear growth in number of clusters
}



#TODO: make this a function of the array id
sparse.graph = F
if(sparse.graph) {
  p.in.seq = rep(0.2*100/n.seq, J)
  p.out.seq = rep(0.02*100/n.seq, J)
} else {
  p.in.seq = rep(0.2, J)
  p.out.seq = rep(0.02, J)
}


# other simulation parameters
equal.size.clusters <- T
binomial.randomization <- F

# In order to have some degree heterogeneity, we say that the clusters have an
# equal spacing of parameters.


# we decouple the design piece and only talk about inference.
n.methods = 8 #three full data regressions, three ard regressions 1 Difference in means and 1 HT

results <- array(NA, c(n.sims, n.methods, J))

#### Simulations Start

# formula for the regression model
fmla <- formula(Y ~ 0 + deg.ratio +  deg.ratio:(H + A + H:A + frac.treated + H:frac.treated))

for(j in seq(J)){
  for(sim in seq(n.sims)){
    n = n.seq[j]
    K = K.seq[j]

    P = matrix(data = 0, nrow = K, ncol = K)

    p.in.max = p.in.seq[j]
    p.in.min = p.in.max/2


    p.out.max = p.out.seq[j]
    p.out.min = p.out.max/2


    P[upper.tri(P)] = seq(p.out.min,p.out.max, length.out = choose(K,2))
    P = P + t(P)
    diag(P) = seq(p.in.min,p.in.max, length.out = K)


    if(equal.size.clusters){
      PI = rep(1/K,K)
    } else {
      PI = seq(1,K,K)
      PI = PI/K
    }

    g.true <- generateSBM(n,P,PI)
    G.true <- g.true$G
    Z.true = g.true$groups

    H <- simH(Z.true)


    # whether to do binomial randomization or the true cluster treatment
    if(binomial.randomization){
      A <- binomialTreatment(n,p.treat)
    } else {
      A <- clusterTreatment(Z.true, p.treat)
      if(include.estimated.randomization){
        A.est  <- clusterTreatment(Z.est, p.treat)
      }
    }


    Y <- simUganderModelOutcomes(G.true,H,A,a = a, b = b, delta = delta, gamma = gamma, sigma = sigma)



    d.vec = as.numeric(G.true %*% rep(1,n))
    d.mean = mean(d.vec)
    f <- fracTreated(G.true,A)
    data <- UganderCovariates(G,H,A)
    data <- data.frame("Y" = Y, "A" = A, "frac.treated" = f, "H" = H, "deg.ratio" = d.vec/d.mean)



    model <- lm(fmla, data = data)


    # full data version.
    gate.est1 <- model$coefficients[3]/model$coefficients[1] + model$coefficients[4]/model$coefficients[2]
    gate.est2 <- model$coefficients[5]/model$coefficients[1] + model$coefficients[6]/model$coefficients[2]
    gate.est3 <- (gate.est1 + gate.est2)/2



    # Here we are asking about the traits directly for simplicity in the simulations.
    #
    traits <- g.true$groups

    X <-computeARD(traits, G.true)


    # Assuming For simplicity, K is known.

    k.means.model <- clusterARDKmeans(X,K)


    Z.hat = k.means.model$cluster
    #
    Z.hat = labelSwitching(Z.true,Z.hat)


    P.hat <- estimatePmat(Z.hat, G.true)

    res.ard <- ARDSBMLinearRegressionSim(Y, fmla, UganderCovariates, A, H, P.hat,Z.hat, B.boot = 200, verbose = T)


    ard.avg.coef <- meanOverList(res.ard$coef)



    ard.gate.est1 <- ard.avg.coef[3]/ard.avg.coef[1] + ard.avg.coef[4]/ard.avg.coef[2]
    ard.gate.est2 <- ard.avg.coef[5]/ard.avg.coef[1] + ard.avg.coef[6]/ard.avg.coef[2]
    ard.gate.est3 <- (ard.gate.est1 + ard.gate.est2)/2


    gate.HT <- HTEstimatorCluster(Y,G.true,A,c.vec = Z.true,frac = p.treat)

    gate.DM <- diffMeansEstimator(Y,A)

    res.vec <- c(ard.gate.est1, ard.gate.est2, ard.gate.est3,
                 gate.est1, gate.est2, gate.est3,
                 gate.HT, gate.DM)

    results[sim,,j] <- res.vec
  }
}


colnames(results) = c("ard1", "ard2", "ard3",
                      "reg1", "reg2", "reg3",
                      "HT", "DM")


filename <- paste0("data/UganderGATE_cluster_growth",cluster.growth,
                   "mutual_benefit_", mutual.benefit,
                   "sparse_graph_",sparse.graph,
                   "block",block,".rds")

saveRDS(results, filename)









