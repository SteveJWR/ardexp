
## Simulations to compare the full data and spillover parameters

rm(list = ls())
source("R/ardexp.R") # functions for the main
source("R/simulation_functions.R") # has the additional simulations pieces



library(reshape2)
library(ggplot2)
library(igraph)
library(blockmodels)
library(sandwich)


slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
print(Sys.getenv('SLURM_ARRAY_TASK_ID'))
if(slurm_arrayid == ""){
  id = 1
} else {
  # coerce the value to an integer
  id <- as.numeric(slurm_arrayid) # for zero indexing
}

set.seed(id)
block = id %% 77  + 1 # repeat every 77


# Maximum number of possible clusters to search for.
K.max = 22
K_set = seq(4,K.max,2)


#K_set = seq(120,300,20)
K_idx = id %/% 77 + 1
K = K_set[K_idx]

#village number
# Total 77 villages

# True model parameters
a = 1
b = -1
delta = 1

# simulation noise
sigma = 0.5 # how much noise are we willing to permit here?


mutual.benefit = T
if(mutual.benefit){
  gamma = 1
} else {
  gamma = -1
}

# public data
# data.file <- paste0("data/JPAL/Data/1. Network Data/Adjacency Matrices/adj_allVillageRelationships_vilno_", block, ".csv")

# real data
data.file <- paste0("DoNotUpload/Network Data/MicroFinance Wave 2/Graphs.mat")


graphs = R.matlab::readMat(data.file)

graphs = graphs$Z[[block]][[1]]
views = length(graphs)
G = NULL
for (view in seq(views)){
  if( is.null(G) ) {
    G = as.matrix(graphs[view][[1]][[1]])
  }

  else {
    if(! is.null(graphs[view][[1]][[1]])){
      G = G + as.matrix(graphs[view][[1]][[1]])
    }
  }
}



#G <- read.csv(data.file, header = FALSE)
#G <- as.matrix(G)
colnames(G) = NULL
rownames(G) = NULL
diag(G) = 0 # Remove self edges from the noisy sample
G <- G + t(G)
G[G > 0] = 1
G.true <- G

g <- graph_from_adjacency_matrix(G.true, mode = "undirected")
# remove unconnected individuals
idx = which(colSums(G.true) > 0)
G.true = G.true[idx,idx]
g <- graph_from_adjacency_matrix(G.true, mode = "undirected")


## Stochastic block-model clustering,
b.model <- BM_bernoulli(membership_type = 'SBM_sym', adj = G.true, explore_min = K.max, explore_max = K.max)
b.model$estimate()
clust_bm_K = apply(b.model$memberships[[K]]$Z, 1,which.max)

P.true <- b.model$model_parameters[[K]]$pi
P.true = (P.true + t(P.true))/2

estimate.sbm <- T
if(estimate.sbm){
  #sbm.true <- sbm::estimateSimpleSBM(G.true)

  # all we need are the "true clusters"
  Z.true <- clust_bm_K
  P.true <- estimatePmat(Z.true,G.true)
  P.true = (P.true + t(P.true))/2

  traits <- Z.true
  X <- computeARD(traits, G.true)
}


# Think of a second version of the simulation where there is a competitive advantage to treatment
# but that wears off if everyone else around you has it



n.sims =  50 # number of simulations per network

# Only compare the full data and partial data versions
n.methods = 5
results <- array(NA, c(n.sims, n.methods))
results_cover <- array(NA, c(n.sims, 2))
# Methods tuning parameters
B.boot = 1000

#### Simulations Start

# formula for the regression model
fmla <- formula(Y ~ 0 + deg.ratio +  deg.ratio:(H + A + H:A + frac.treated + H:frac.treated))


# SaturationRandomization Design.
sat.frac = 1/2 #treat half of the clusters at 0.9 and the other half at 0.1
p.high.level = 0.9
p.low.level = 0.1

n = nrow(G.true)
K = max(Z.true) # number of true clusters
H <- simH(Z.true)
for(sim in seq(n.sims)){
  set.seed(sim)
  cat(paste0("Sim: ", sim, "/", n.sims), end = "\r")

  # whether to do binomial randomization or the true cluster treatment
  K.high.level <- round(sat.frac*K)
  K.low.level <- K - K.high.level
  #clust.levels <- c(rep(p.high.level, K.high.level), rep(p.low.level, K.low.level))
  diag.clusters = diag(P.true)
  clust.levels <- rep(c(p.high.level,p.low.level),times = K)
  clust.levels = clust.levels[1:K]
  clust.levels <- clust.levels[order(diag.clusters)]

  A.sat <- saturationRandomizationTreatment(Z.true, levels = clust.levels)

  #outcome.sat <- simSpilloverModelOutcomes(G.true,A.sat,a = a, b = b, delta = delta, gamma0 = gamma0, gamma1 = gamma1, sigma = sigma)
  outcome.sat <- simUganderModelOutcomes(G.true, H, A.sat, a = a, b = b, delta = delta, gamma = gamma, sigma = sigma)

  Y.sat <- outcome.sat$Y
  true.gate = outcome.sat$GATE
  # print(true.gate)

  d.vec = as.numeric(G.true %*% rep(1,n))
  d.mean = mean(d.vec)
  f <- fracTreated(G.true,A.sat)
  data <- UganderCovariates(G.true,H,A.sat)
  data <- data.frame("Y" = Y.sat, "A" = A.sat, "frac.treated" = f, "deg.ratio" = d.vec/d.mean)


  model <- lm(fmla, data = data)

  X.data.0 <- data.frame("A" = rep(0,length(A.sat)), "frac.treated" = rep(0,length(A.sat)), "deg.ratio" = d.vec/d.mean)
  X.data.1 <- data.frame("A" = rep(1,length(A.sat)), "frac.treated" = rep(1,length(A.sat)), "deg.ratio" = d.vec/d.mean)

  X.data.0$Y = Y.sat
  X.data.1$Y = Y.sat

  y.0 <- predict(model, X.data.0)
  y.1 <- predict(model, X.data.1)

  # full data regression
  gate.est.regression <- mean((y.1 - y.0))

  # Inference for GATE
  phi = colMeans(model.matrix(fmla, X.data.1)) - colMeans(model.matrix(fmla, X.data.0))

  Sigma <- sandwich::sandwich(model)
  gate.var <- t(phi) %*% Sigma %*% phi
  sd.gate.full.data = sqrt(gate.var)
  cover.full.data <- abs(true.gate - gate.est.regression)/sd.gate.full.data <= 1.96

  # Here we are asking about the traits directly for simplicity in the simulations.


  # Assuming For simplicity, K is known.
  k.means.model <- clusterARDKmeans(X,K)

  Z.hat = k.means.model$cluster
  #
  Z.hat = labelSwitching(Z.true,Z.hat)

  # Z.hat = Z.true # true clusters known
  # estimate
  P.hat <- estimatePmatARD(Z.hat, X)


  # Estimated MODEL GATE
  res.ard <- ARDUganderSBMLinearRegressionSim(Y.sat, fmla,  A.sat, H, P.hat,Z.hat, B.boot = B.boot, verbose = T)

  mean.data = meanOverList(res.ard$data)

  ard.mean.model <- lm(fmla, data = mean.data)


  mean.data.0 = mean.data
  mean.data.1 = mean.data

  mean.data.0$A = 0
  mean.data.0$frac.treated = 0

  mean.data.1$A = 1
  mean.data.1$frac.treated = 1

  gate.ard.est <- mean(predict(ard.mean.model, mean.data.1) - predict(ard.mean.model, mean.data.0))

  # Inference for GATE
  phi = colMeans(model.matrix(fmla, mean.data.1)) - colMeans(model.matrix(fmla, mean.data.0))

  Sigma <- sandwich::sandwich(ard.mean.model)

  gate.var <- t(phi) %*% Sigma %*% phi
  sd.gate = sqrt(gate.var)
  # graph.sample.sd <- sqrt(meanOverList((unlist(res.ard$gate)- mean(unlist(res.ard$gate)))**2))

  cover <- abs(true.gate - gate.ard.est)/sd.gate <= 1.96
  # cover.model.sample <- abs(true.gate - gate.ard.est)/sqrt(sd.gate**2 + graph.sample.sd**2) <= 1.96

  # 'True' MODEL GATE
  res.ard <- ARDUganderSBMLinearRegressionSim(Y.sat, fmla,  A.sat, H, P.hat,Z.hat, B.boot = B.boot, verbose = T)
  mean.data = meanOverList(res.ard$data)
  ard.mean.model <- lm(fmla, data = mean.data)

  mean.data.0 = mean.data
  mean.data.1 = mean.data
  mean.data.0$A = 0
  mean.data.0$frac.treated = 0
  mean.data.1$A = 1
  mean.data.1$frac.treated = 1

  gate.ard.true.model.est <- mean(predict(ard.mean.model, mean.data.1) - predict(ard.mean.model, mean.data.0))




  A.cluster.treat <- clusterTreatment(Z.true, p.treat = 1/2)

  outcome.cluster.treat <- simUganderModelOutcomes(G.true, H, A.cluster.treat, a = a, b = b, delta = delta, gamma = gamma, sigma = sigma)
  Y.cluster.treat <- outcome.cluster.treat$Y

  DM.est <- diffMeansEstimator(Y.cluster.treat,A.cluster.treat)


  HT.est <- HTEstimatorCluster(Y.cluster.treat, G.true, A.cluster.treat, Z.true, p.treat = 1/2)

  res.vec <- c(gate.ard.est - true.gate,
               gate.ard.true.model.est - true.gate,
               gate.est.regression - true.gate,
               DM.est - outcome.cluster.treat$GATE,
               HT.est - outcome.cluster.treat$GATE)

  cover.vec <- c(cover[1,1], cover.full.data[1,1])




  names(res.vec) = NULL
  names(cover.vec) = NULL
  results[sim,] <- res.vec
  results_cover[sim,] <- cover.vec
}
# TODO: Problems seem to arize in the extreme values, there are seemingly sometimes occurances when the ARD estimate is way off =
# results_cover_subset = results_cover[abs(results[,1]) < 3,]
# colMeans(results_cover_subset)

colnames(results) = c("ard",
                      "ard.tm",
                      "reg",
                      "DM",
                      "HT")

colnames(results_cover) = c("cover asymptotic",
                            "cover full data")

print(colMeans(round(results,6)))

filename <- paste0("data/JPAL_sim_results/GATE_JPAL_village_",block,"_K_", K,".rds")
print(colMeans(round(results,6)))
saveRDS(results, filename)

filename_coverage <- paste0("data/JPAL_sim_results/GATE_JPAL_village_",block,"_K_", K,"_coverage.rds")

saveRDS(results_cover, filename_coverage)


### Plotting for debugging
# plot(g, vertex.color=Z.true, vertex.label = NA)
# n = nrow(G.true)
# g.sim <- generateSBM(n,P = P.true,Z = Z.true)
# G = g.sim$G
# plot(graph_from_adjacency_matrix(g.sim$G, mode = "undirected"), vertex.color=Z.true, vertex.label = NA)
# #
# # Yes, that is the problem right now!
# P.hat1 <- estimatePmat(Z.true,G)
# P.hat1 = (P.hat1 + t(P.hat1))/2
# print(P.hat1)
# print(P.true )
# print(P.hat1 - P.true)
# plot(g.sim, vertex.color=Z.true, vertex.label = NA)
# table(clust_greedy_K)
