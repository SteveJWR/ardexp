rm(list = ls())

source("R/ardexp.R") # functions for the main
source("R/simulation_functions.R") # has the additional simulations pieces



library(reshape2)
library(ggplot2)
library(igraph)
library(blockmodels)
library(sandwich)

library(haven)
library(R.matlab)



slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
print(Sys.getenv('SLURM_ARRAY_TASK_ID'))
if(slurm_arrayid == ""){
  id = 68
} else {
  # coerce the value to an integer
  id <- as.numeric(slurm_arrayid) # for zero indexing
}

set.seed(id) #

block = id %% 77  + 1 # repeat every 77

# Maximum number of possible clusters to search for.
K.max = 6 # when setting this manually to 6 we can keep this simple
K_set = seq(4,K.max,2)


#K_set = seq(120,300,20)
K_idx = id %/% 77 + 1
K = K_set[K_idx]

# For the sake of simplicity, we assume that we have a 6 block SBM estimate in each of the networks.
# Fixing this in advance will be necessary for the optimal design step.
K = 6



#village number
# Total 77 villages

# True model parameters
alpha0 = -1
alpha1 = 1
T.max = 3
q.vec = c(0.05,0.005,0.0005) # fading diffusion faster than simple diffusion



# public data

on.cluster = T
if(on.cluster){
  hh_dat <- readRDS('data/hh_dat.rds')
  network_set <- readRDS('data/network_set.rds')
  cell_dat <- readRDS('data/cell_dat.rds')
  vertex_key <- readRDS('data/vertex_key.rds')

} else {

  data.file.networks <- paste0("DoNotUpload/Network Data/Gossip Data/RFENetwork.mat")
  data.file.vertex_key <- paste0("DoNotUpload/Network Data/Gossip Data/vertex_key.mat")
  # seed info
  data.file.seed <- paste0("DoNotUpload/Network Data/Gossip Data/all_hh_data_18_Nov_2016_v9.dta")
  data.file.cell <- paste0("DoNotUpload/Network Data/Gossip Data/karnataka_cell_rct.dta")
  # All data in this case is from comparing to a randomly seeded example in the treatments.


  hh_dat <- read_dta(data.file.seed) # no number 27
  network_set <- readMat(data.file.networks)
  vertex_key <- readMat(data.file.vertex_key)
  cell_dat <- read_dta(data.file.cell) # no number 27

  saveRDS(hh_dat, 'data/hh_dat.rds')
  saveRDS(network_set, 'data/network_set.rds')
  saveRDS(vertex_key, 'data/vertex_key.rds')
  saveRDS(cell_dat, 'data/cell_dat.rds')
}
# real data networks





print(colnames(hh_dat))


village_dat = hh_dat[hh_dat$villageid == block,]
seeds_village_hhid = village_dat$hhid[village_dat$seed_dummy == 1]


# which seeds were actually treated
village_key <- vertex_key$vertex.to.hhid.key[[block]][[1]]
seeds <- rep(0,nrow(village_key))

seed_idx = village_key[village_key[,2] %in% seeds_village_hhid,1]
seeds[seed_idx] = 1
## Simulations to compare the full data and spillover parameters



graphs = network_set$Z[[block]][[1]]
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

# g <- graph_from_adjacency_matrix(G.true, mode = "undirected")
# remove unconnected individuals
idx = which(colSums(G.true) > 0)
G.true = G.true[idx,idx]
seeds.true = seeds[idx]
g <- graph_from_adjacency_matrix(G.true, mode = "undirected")


# TODO: Updated using the Leiden Method
## Stochastic block-model clustering,
# b.model <- BM_bernoulli(membership_type = 'SBM_sym', adj = G.true, explore_min = K.max, explore_max = K.max)
# b.model$estimate()
# clust_bm_K = apply(b.model$memberships[[K]]$Z, 1,which.max)
clust_bm_K = cluster_leiden_K(G.true, K = K)

# table(clust_bm_K)
# P.true <- b.model$model_parameters[[K]]$pi
# P.true = (P.true + t(P.true))/2

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

n.sims =  500 # number of simulations per network

# Only compare the full data and partial data versions
n.methods = 3 # compare the alternative under SBM

# we only care about a single parameter of interest
results <- array(NA, c(n.sims, n.methods))
results_cover <- array(NA, c(n.sims, n.methods))
results_sd <- array(NA, c(n.sims, n.methods))
results_all_params <- array(NA, c(n.sims, n.methods))

# Methods tuning parameters
B.boot = 1000 # 200 # number of bootstrap samples

#### Simulations Start

# formula for the regression model
fmla <- formula(Y ~ n0 + n1 + n2 + n3) # treatment fractions of each of the neighbors up to 3 steps


#
true.coefficients <- c(-alpha0,alpha1,alpha1*c(q.vec[1], q.vec[1]*q.vec[2], q.vec[1]*q.vec[2]*q.vec[3]))


n = nrow(G.true)
K = max(Z.true) # number of true clusters

T.max = length(q.vec)



A.seed <- as.numeric(seeds.true)
# random walk coefficients
X_walk = DiffusionExampleCovariates(G.true,A.seed, T.max = T.max)
X_walk <- as.data.frame(X_walk)
colnames(X_walk) = paste0("n", 0:3)

#optimal.design.saturations <- matrix(nrow = 72, ncol = 6)
# Precomputed saturations levels for each of the relevant clusters
optimal.design.saturations <- readRDS(file = 'data/optimal_design_saturations.rds')
missing.row = sum(is.na(optimal.design.saturations[block,])) > 0
missing.row = F
if(missing.row){
  # optimal design seed example:
  library(rBayesianOptimization)

  # Assuming For simplicity, K is known.
  k.means.model <- clusterARDKmeans(X,K)

  Z.hat = k.means.model$cluster
  #
  Z.hat = labelSwitching(Z.true,Z.hat)

  # Z.hat = Z.true # true clusters known
  # estimate
  P.hat <- estimatePmatARD(Z.hat, X)

  G_set <- Generate_G_set(P.hat, Z.hat, L = 100)

  n = nrow(G.true)
  Sigma_naive = diag(rep(1,n)) # working variance, will not be exactly the same in the logistic example, but can still be useful

  #A = rbinom(n, size = 1, prob = 0.5)
  clusters = Z.true
  regression_features = DiffusionExampleCovariates
  phi = c(0,1,0,0,0)



  V_tau <- variance_function_factory(phi,G_set,clusters,regression_features, Sigma_naive)

  # package framework requires specifying the exact number of parameters
  Bayesopt_obj <- function(tau1, tau2, tau3, tau4, tau5, tau6){
    tau = c(tau1, tau2, tau3, tau4, tau5, tau6)
    result = list(Score = -V_tau(tau), pred = 0) # since the default is to maximize
    return(result)
  }
  search_bound <- list(tau1 = c(0,1), tau2 = c(0,1), tau3 = c(0,1), tau4 = c(0,1), tau5 = c(0,1), tau6 = c(0,1))


  search_grid <- data.frame(tau1 = runif(20,0,1),
                            tau2 = runif(20,0,1),
                            tau3 = runif(20,0,1),
                            tau4 = runif(20,0,1),
                            tau5 = runif(20,0,1),
                            tau6 = runif(20,0,1)
  )


  an.error.occured <- c()
  tryCatch( { bayes_opt_variance <- BayesianOptimization(FUN = Bayesopt_obj, bounds = search_bound,
                                                         init_grid_dt = search_grid, init_points = 0,
                                                         n_iter = 25, acq = 'ei'); }
            , error = function(e) {an.error.occured <<- c(an.error.occured, block)})
  print(an.error.occured)

  optimal.design.saturations[block,] = as.numeric(bayes_opt_variance$Best_Par)
  saveRDS(optimal.design.saturations, file = 'data/optimal_design_saturations.rds')
  print(paste0("Block ", block, ' of 72 complete'))
}






# optimal design saturation levels from the precomputed levels.
tau.opt = optimal.design.saturations[block,]
A.opt.sat <- saturation_random_sample(tau.opt, Z.true)



for(sim in seq(n.sims)){
  cat(paste0("Sim: ", sim, "/", n.sims), end = "\r")

  # whether to do binomial randomization or the true cluster treatment
  #K.high.level <- round(sat.frac*K)
  #K.low.level <- K - K.high.level
  #clust.levels <- c(rep(p.high.level, K.high.level), rep(p.low.level, K.low.level))


  #outcome.sat <- simSpilloverModelOutcomes(G.true,A.sat,a = a, b = b, delta = delta, gamma0 = gamma0, gamma1 = gamma1, sigma = sigma)
  #outcome.sat <- simUganderModelOutcomes(G.true, H, A.sat, a = a, b = b, delta = delta, gamma = gamma, sigma = sigma)
  df.sim = DiffusionExample(G.true,A.seed,q.vec, alpha0 = alpha0, alpha1 = alpha1)

  data = X_walk
  data$Y = df.sim$Y

  # Create a logistic regression based-diffusion model
  diff.model <- glm(fmla, data = data, family = "binomial")


  coef.est <- diff.model$coefficients

  # Inference for alpha 1
  phi = c(0,1,0,0,0)

  Sigma <- sandwich::sandwich(diff.model)
  alpha1.var <- t(phi) %*% Sigma %*% phi
  sd.alpha1.full.data <- sqrt(alpha1.var)
  cover.full.data <- abs(alpha1 - coef.est[2])/sd.alpha1.full.data <= 1.96

  results.diff.full.data  <- coef.est - true.coefficients


  # Assuming For simplicity, K is known.
  k.means.model <- clusterARDKmeans(X,K)

  Z.hat = k.means.model$cluster
  #
  Z.hat = labelSwitching(Z.true,Z.hat)

  Z.hat = Z.true # true clusters known
  # estimate
  P.hat <- estimatePmatARD(Z.hat, X)


  # Estimated MODEL GATE
  X.samp.ard = ARDDiffusionRegressionSim(A.seed, P.hat, Z.hat, T.max = T.max, B.boot = B.boot, verbose = F)


  X.mean.data = data.frame(meanOverList(X.samp.ard$data))

  colnames(X.mean.data) = paste0("n", 0:3)
  X.mean.data$Y =  df.sim$Y

  ard.mean.model <- glm(fmla, data = X.mean.data, family = "binomial")

  ard.coef.est <- ard.mean.model$coefficients

  # Inference for alpha 1
  phi = c(0,1,0,0,0)

  Sigma <- sandwich::sandwich(ard.mean.model)
  alpha1.var <- t(phi) %*% Sigma %*% phi
  sd.alpha1.ard <- sqrt(alpha1.var)
  cover.ard <- abs(alpha1 - ard.coef.est[2])/sd.alpha1.ard <= 1.96

  sd.length.ard = sd.alpha1.ard
  results.diff.ard  <- ard.coef.est - true.coefficients


  df.sim.opt = DiffusionExample(G.true,A.opt.sat,q.vec, alpha0 = alpha0, alpha1 = alpha1)


  X.mean.data.opt = ARDDiffusionRegressionSim(A.opt.sat, P.hat, Z.hat, T.max = T.max, B.boot = B.boot, verbose = F)


  X.mean.data.opt = data.frame(meanOverList(X.mean.data.opt$data))

  colnames(X.mean.data.opt) = paste0("n", 0:3)
  X.mean.data.opt$Y =  df.sim.opt$Y

  ard.mean.model.opt <- glm(fmla, data = X.mean.data.opt, family = "binomial")


  ard.coef.est.opt <- ard.mean.model.opt$coefficients

  # Inference for alpha 1
  phi = c(0,1,0,0,0)

  Sigma <- sandwich::sandwich(ard.mean.model.opt)
  alpha1.var <- t(phi) %*% Sigma %*% phi
  sd.alpha1.ard.opt <- sqrt(alpha1.var)
  cover.ard.opt <- abs(alpha1 - ard.coef.est.opt[2])/sd.alpha1.ard.opt <= 1.96

  sd.length.ard = sd.alpha1.ard
  results.diff.ard.opt  <- ard.coef.est.opt - true.coefficients



  # optimal performance of the estimate for

  res.vec <- c(ard.coef.est[2] - alpha1,
               ard.coef.est.opt[2] - alpha1,
               coef.est[2] - alpha1)

  cover.vec <- c(cover.ard.opt[1,1],cover.full.data[1,1],cover.ard[1,1])

  sd.vec <-  c(sd.alpha1.ard,sd.alpha1.ard.opt,sd.alpha1.full.data)

  res.vec.2norm <- c(sqrt(sum((ard.coef.est - true.coefficients)**2)),
                     sqrt(sum((ard.coef.est.opt - true.coefficients)**2)),
                     sqrt(sum((coef.est - true.coefficients)**2)))


  names(res.vec) = NULL
  names(cover.vec) = NULL
  names(sd.vec) = NULL
  results[sim,] <- res.vec
  results_cover[sim,] <- cover.vec
  results_sd[sim, ] <- sd.vec
  results_all_params[sim,] <- res.vec.2norm
}

colnames(results) = c("ard",
                      "ard.opt",
                      "reg.full.data")

colnames(results_cover) = c("ard",
                            "ard.opt",
                            "reg.full.data")

colnames(results_sd) = c("ard",
                         "ard.opt",
                         "reg.full.data")

colnames(results_all_params) = c("ard",
                         "ard.opt",
                         "reg.full.data")



filename <- paste0("data/Gossip_sim_results/Gossip_village_",block,"_diff.rds")
saveRDS(results, filename)

filename_coverage <- paste0("data/Gossip_sim_results/Gossip_village_",block,"_coverage.rds")
saveRDS(results_cover, filename_coverage)

filename_sd <- paste0("data/Gossip_sim_results/Gossip_village_",block,"_sd.rds")
saveRDS(results_sd, filename_sd)

filename_all_params <- paste0("data/Gossip_sim_results/Gossip_village_",block,"_all_params.rds")
saveRDS(results_all_params, filename_all_params)






