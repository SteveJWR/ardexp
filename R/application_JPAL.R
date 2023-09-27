

## Simulations to compare the full data and spillover parameters

rm(list = ls())
source("R/ardexp.R") # functions for the main
source("R/simulation_functions.R") # has the additional simulations pieces



library(reshape2)
library(ggplot2)
library(igraph)
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

K_set = seq(70,100,5) # There are 20 on this scale.
K_set = seq(120,300,20)
K_idx = id %/% 77 + 1
K = K_set[K_idx]
#village number
# Total 77 villages


# True model parameters
a = 1
b = 1
delta = 1

# simulation noise
sigma = 0.5 # how much noise are we willing to permit here?


mutual.benefit = T

if(mutual.benefit){
  gamma0 = 1
  gamma1 = 0.2
} else {
  gamma0 = -1
  gamma1 = -0.2
}

data.file <- paste0("data/JPAL/Data/1. Network Data/Adjacency Matrices/adj_allVillageRelationships_vilno_", block, ".csv")



G <- read.csv(data.file, header = FALSE)
G <- as.matrix(G)
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

clust_greedy = cluster_fast_greedy(g)

clust_greedy_K = clust_greedy %>% as.hclust() %>% cutree(K)

estimate.sbm <- T
if(estimate.sbm){
  #sbm.true <- sbm::estimateSimpleSBM(G)

  # all we need are the "true clusters"
  Z.true <- clust_greedy_K
  P.true <- estimatePmat(Z.true,G.true)
  P.true = (P.true + t(P.true))/2

  traits <- Z.true
  X <- computeARD(traits, G.true)
}


# Think of a second version of the simulation where there is a competative advantage to treatment
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
fmla <- formula(Y ~ A*frac.treated + deg.ratio) #TODO: fix the name
#fmla <- formula(Y ~ (H + A + H:A + frac.treated + H:frac.treated))


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
  # A.clust <- clusterTreatment(Z.true, p.treat) #TODO: remove this



  K.high.level <- round(sat.frac*K)
  K.low.level <- K - K.high.level
  clust.levels <- c(rep(p.high.level, K.high.level), rep(p.low.level, K.low.level))
  A.sat <- saturationRandomizationTreatment(Z.true, levels = clust.levels)

  outcome.sat <- simSpilloverModelOutcomes(G.true,A.sat,a = a, b = b, delta = delta, gamma0 = gamma0, gamma1 = gamma1, sigma = sigma)
  #Y.clust <- outcome.clust$Y #TODO: delete this
  Y.sat <- outcome.sat$Y


  d.vec = as.numeric(G.true %*% rep(1,n))
  d.mean = mean(d.vec)
  f <- fracTreated(G.true,A.sat)
  data <- SpilloverCovariates(G.true,A.sat)
  data <- data.frame("Y" = Y.sat, "A" = A.sat, "frac.treated" = f, "deg.ratio" = d.vec/d.mean)


  model <- lm(fmla, data = data)



  # full data version.
  # (Note I had e)
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

  ard.avg.coef <- meanOverList(res.ard$coef)
  ard.avg.data <- meanOverList(res.ard$data)

  ard.theta.est <- ard.avg.coef[5]

  res.ard.true.model <- ARDSBMSpilloverLinearRegressionSim(Y.sat, fmla, SpilloverCovariates, A.sat, P.true, Z.true, B.boot = B.boot, verbose = T)
  ard.true.model.avg.coef <- meanOverList(res.ard.true.model$coef)

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

filename <- paste0("data/JPAL_sim_results/JPAL_village_",block,"_K_", K,".rds")

saveRDS(results, filename)
