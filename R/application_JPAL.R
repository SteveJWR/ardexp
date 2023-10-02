

## Simulations to compare the full data and spillover parameters

rm(list = ls())
source("R/ardexp.R") # functions for the main
source("R/simulation_functions.R") # has the additional simulations pieces



library(reshape2)
library(ggplot2)
library(igraph)
library(blockmodels)
# library(lsa) # for cosine similarity


slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
print(Sys.getenv('SLURM_ARRAY_TASK_ID'))
if(slurm_arrayid == ""){
  id = 1
} else {
  # coerce the value to an integer
  id <- as.numeric(slurm_arrayid) + 1
}


block = id %% 77  # repeat every 77


K.max = 22
K_set = seq(4,K.max,2)


#K_set = seq(120,300,20)
K_idx = id %/% 77 + 1
K = K_set[K_idx]

#village number
# Total 77 villages


# True model parameters
a = 1
b = 1
delta = 1

# simulation noise (variance)
sigma = 0.05 # how much noise are we willing to permit here?


mutual.benefit = T

if(mutual.benefit){
  gamma0 = 1
  gamma1 = 0.2
} else {
  gamma0 = -1
  gamma1 = -0.2
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


## Stochastic block-model clustering
b.model <- BM_bernoulli(membership_type = 'SBM_sym', adj = G.true, explore_max = K.max)
b.model$estimate()
clust_bm_K = apply(b.model$memberships[[K]]$Z, 1,which.max)


# table(clust_bm_K)
P.true <- b.model$model_parameters[[K]]$pi
P.true = (P.true + t(P.true))/2

estimate.sbm <- T
if(estimate.sbm){
  #sbm.true <- sbm::estimateSimpleSBM(G)

  # all we need are the "true clusters"
  Z.true <- clust_bm_K
  P.true <- estimatePmat(Z.true,G.true)
  P.true = (P.true + t(P.true))/2

  traits <- Z.true
  X <- computeARD(traits, G.true)
}




# Think of a second version of the simulation where there is a competitive advantage to treatment
# but that wears off if everyone else around you has it

#Global average treatment effect
theta = gamma1 - gamma0

n.sims =  100#

# Only compare the full data and partial data versions
n.methods = 3
results <- array(NA, c(n.sims, n.methods))

# Methods tuning parameters
B.boot = 50
#### Simulations Start

# formula for the regression model
fmla <- formula(Y ~ A*frac.treated + deg.ratio)
#fmla <- formula(Y ~ A*frac.treated + deg.ratio)

# SaturationRandomization Design.
sat.frac = 1/2 #treat half of the clusters at 0.9 and the other half at 0.1
p.high.level = 0.9
p.low.level = 0.1


for(sim in seq(n.sims)){
  cat(paste0("Sim: ", sim, "/", n.sims), end = "\r")
  n = nrow(G.true)
  K = max(Z.true) # number of true clusters



  #H <- simH(Z.true)
  H = rep(NA,n)

  # whether to do binomial randomization or the true cluster treatment
  K.high.level <- round(sat.frac*K)
  K.low.level <- K - K.high.level
  clust.levels <- c(rep(p.high.level, K.high.level), rep(p.low.level, K.low.level))
  A.sat <- saturationRandomizationTreatment(Z.true, levels = clust.levels)

  outcome.sat <- simSpilloverModelOutcomes(G.true,A.sat,a = a, b = b, delta = delta, gamma0 = gamma0, gamma1 = gamma1, sigma = sigma)

  Y.sat <- outcome.sat$Y


  d.vec = as.numeric(G.true %*% rep(1,n))
  d.mean = mean(d.vec)
  f <- fracTreated(G.true,A.sat)
  data <- SpilloverCovariates(G.true,A.sat)
  data <- data.frame("Y" = Y.sat, "A" = A.sat, "frac.treated" = f, "deg.ratio" = d.vec/d.mean)


  model <- lm(fmla, data = data)



  # full data version.
  theta.est <- model$coefficients[5]


  # Here we are asking about the traits directly for simplicity in the simulations.


  # Assuming For simplicity, K is known.
  k.means.model <- clusterARDKmeans(X,K)


  Z.hat = k.means.model$cluster
  #
  Z.hat = labelSwitching(Z.true,Z.hat)

  Z.hat = Z.true # true clusters known
  # estimate
  P.hat <- estimatePmatARD(Z.hat, X)


  res.ard <- ARDSBMSpilloverLinearRegressionSim(Y.sat, fmla, SpilloverCovariates, A.sat, P.hat,Z.hat, B.boot = B.boot, verbose = T)

  ard.mean.model <- lm(fmla, data = meanOverList(res.ard$data))
  ard.model.avg.coef <- ard.mean.model$coef
  ard.theta.est <- ard.model.avg.coef[5]

  res.ard.true <- ARDSBMSpilloverLinearRegressionSim(Y.sat, fmla, SpilloverCovariates, A.sat, P.true, Z.true, B.boot = B.boot, verbose = T)

  ard.true.mean.model <- lm(fmla, data = meanOverList(res.ard.true$data))
  ard.true.model.avg.coef <- ard.true.mean.model$coef
  ard.true.model.theta.est <- ard.true.model.avg.coef[5]

  res.vec <- c(ard.theta.est - theta,
               ard.true.model.theta.est - theta,
               theta.est - theta)


  names(res.vec) = NULL
  results[sim,] <- res.vec
}

colnames(results) = c("ard",
                      "ard.tm",
                      "reg")



print(colMeans(round(results,6)))

filename <- paste0("data/JPAL_sim_results/JPAL_village_",block,"_K_", K,".rds")

saveRDS(results, filename)

# TODO: Add the comment

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

