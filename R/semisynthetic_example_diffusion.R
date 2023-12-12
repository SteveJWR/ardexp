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

# real data networks
data.file.networks <- paste0("DoNotUpload/Network Data/Gossip Data/RFENetwork.mat")

# seed info
data.file.seed <- paste0("DoNotUpload/Network Data/Gossip Data/all_hh_data_18_Nov_2016_v9.dta")

data.file.cell <- paste0("DoNotUpload/Network Data/Gossip Data/karnataka_cell_rct.dta")
# All data in this case is from comparing to a randomly seeded example in the treatments.


hh_dat <- haven::read_dta(data.file.seed) # no number 27
network_set <- R.matlab::readMat(data.file.networks)
cell_dat <- haven::read_dta(data.file.cell) # no number 27



village_dat = hh_dat %>% filter(villageid == block)
seeds  <- village_dat$seed_dummy

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

g <- graph_from_adjacency_matrix(G.true, mode = "undirected")
# remove unconnected individuals
idx = which(colSums(G.true) > 0)
G.true = G.true[idx,idx]
seeds.true = seeds[idx]
g <- graph_from_adjacency_matrix(G.true, mode = "undirected")


## Stochastic block-model clustering,
b.model <- BM_bernoulli(membership_type = 'SBM_sym', adj = G.true, explore_min = K.max, explore_max = K.max)
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



n.sims =  500 # number of simulations per network

# Only compare the full data and partial data versions
n.methods = 3 # compare the alternative under SBM
n.coefs = T.max + 2
results <- array(NA, c(n.sims, n.methods, n.coefs))
results_cover <- array(NA, c(n.sims, 3))
results_se_length <- array(NA, c(n.sims, n.methods))
# Methods tuning parameters
B.boot = 200

#### Simulations Start

# formula for the regression model
fmla <- formula(Y ~ n0 + n1 + n2 + n3) # treatment fractions of each of the neighbors up to 3 steps


#
true.coefficients <- c(-alpha0,alpha1,alpha1*c(q.vec[1], q.vec[1]*q.vec[2], q.vec[1]*q.vec[2]*q.vec[3]))


n = nrow(G.true)
K = max(Z.true) # number of true clusters


A.seed <- as.numeric(seeds.true)
# random walk coefficients
X_walk = DiffusionExampleCovariates(G.true,A.seed, T.max = T.max)
X_walk <- as.data.frame(X_walk)
colnames(X_walk) = paste0("n", 0:3)

#optimal.design.saturations <- matrix(nrow = 72, ncol = 6)
optimal.design.saturations <- readRDS(file = 'data/optimal_design_saturations.rds')
missing.row = sum(is.na(optimal.design.saturations[block,])) > 0
if(missing.row){
  # optimal design seed example:
  library(rBayesianOptimization)

  # Assuming For simplicity, K is known.
  k.means.model <- clusterARDKmeans(X,K)

  Z.hat = k.means.model$cluster
  #
  Z.hat = labelSwitching(Z.true,Z.hat)

  Z.hat = Z.true # true clusters known
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







tau.opt = optimal.design.saturations[block,]
A.opt.sat <- saturation_random_sample(tau, clusters)



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
  data$y = df.sim$Y

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


  X.samp.ard.opt = ARDDiffusionRegressionSim(A.opt.sat, P.hat, Z.hat, T.max = T.max, B.boot = B.boot, verbose = F)


  X.mean.data.opt = data.frame(meanOverList(X.samp.ard.opt$data))

  colnames(X.mean.data.opt) = paste0("n", 0:3)
  X.mean.data.opt$Y =  df.sim.opt$Y

  ard.mean.model.opt <- glm(fmla, data = X.mean.data.opt, family = "binomial")
  print(ard.mean.model.opt)

  ard.coef.est.opt <- ard.mean.model.opt$coefficients

  # Inference for alpha 1
  phi = c(0,1,0,0,0)

  Sigma <- sandwich::sandwich(ard.mean.model.opt)
  alpha1.var <- t(phi) %*% Sigma %*% phi
  sd.alpha1.ard.opt <- sqrt(alpha1.var)
  cover.ard.opt <- abs(alpha1 - ard.coef.est.opt[2])/sd.alpha1.ard.opt <= 1.96

  sd.length.ard = sd.alpha1.ard
  results.diff.ard  <- ard.coef.est - true.coefficients



  # optimal performance of the estimate for

  res.vec <- c(ard.coef.est[2] - alpha1,
               ard.coef.est.opt[2] - alpha1,
               coef.est[2] - alpha1)

  cover.vec <- c( cover.ard.opt[1,1],cover.full.data[1,1],cover.ard[1,1])

  sd.vec <-  c(sd.alpha1.ard,sd.alpha1.ard.opt,sd.alpha1.full.data)



  names(res.vec) = NULL
  names(cover.vec) = NULL
  names(sd.vec) = NULL
  results[sim,] <- res.vec
  results_cover[sim,] <- cover.vec
  results_sd[sim, ] <- sd.vec
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


print(colMeans(round(results,6)))

filename <- paste0("data/Gossip_sim_results/Gossip_village_",block,"_diff.rds")
print(colMeans(round(results,6)))

saveRDS(results, filename)

filename_coverage <- paste0("data/Gossip_sim_results/Gossip_village_",block,"_coverage.rds")

saveRDS(results_cover, filename_coverage)

filename_sd <- paste0("data/Gossip_sim_results/Gossip_village_",block,"_sd.rds")

saveRDS(results_sd, filename_sd)

# filename_sd <- paste0("data/JPAL_sim_results/GATE_JPAL_village_",block,"_K_", K,"_sd.rds")
#
# saveRDS(results_sd, filename_sd)
#

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


# TODO: remove the unnecessary statement here when done testing.
#clust_bm_K <- c(1,rep(seq(1,nrow(G.true)/2), 2))
#clust_bm_K[nrow(G.true)] = nrow(G.true) - 1

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



n.sims =  20 # number of simulations per network

# Only compare the full data and partial data versions
n.methods = 5
results <- array(NA, c(n.sims, n.methods))
results_cover <- array(NA, c(n.sims, 3))
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
  cat(paste0("Sim: ", sim, "/", n.sims), end = "\r")

  # whether to do binomial randomization or the true cluster treatment
  K.high.level <- round(sat.frac*K)
  K.low.level <- K - K.high.level
  clust.levels <- c(rep(p.high.level, K.high.level), rep(p.low.level, K.low.level))
  A.sat <- saturationRandomizationTreatment(Z.true, levels = clust.levels)

  #outcome.sat <- simSpilloverModelOutcomes(G.true,A.sat,a = a, b = b, delta = delta, gamma0 = gamma0, gamma1 = gamma1, sigma = sigma)
  outcome.sat <- simUganderModelOutcomes(G.true, H, A.sat, a = a, b = b, delta = delta, gamma = gamma, sigma = sigma)

  Y.sat <- outcome.sat$Y
  true.gate = outcome.sat$GATE

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

  Z.hat = Z.true # true clusters known
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
  graph.sample.sd <- sqrt(meanOverList((unlist(res.ard$gate)- mean(unlist(res.ard$gate)))**2))

  cover <- abs(true.gate - gate.ard.est)/sd.gate <= 1.96
  cover.model.sample <- abs(true.gate - gate.ard.est)/sqrt(sd.gate**2 + graph.sample.sd**2) <= 1.96

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

  cover.vec <- c(cover[1,1],cover.model.sample[1,1], cover.full.data[1,1])




  names(res.vec) = NULL
  names(cover.vec) = NULL
  results[sim,] <- res.vec
  results_cover[sim,] <- cover.vec
}

colnames(results) = c("ard",
                      "ard.tm",
                      "reg",
                      "DM",
                      "HT")
colnames(results_cover) = c("cover asymptotic",
                            "cover asymptotic + model sample")

print(colMeans(round(results,6)))

filename <- paste0("data/JPAL_sim_results/GATE_JPAL_village_",block,"_K_", K,".rds")
print(colMeans(round(results,6)))
saveRDS(results, filename)

filename_coverage <- paste0("data/JPAL_sim_results/GATE_JPAL_village_",block,"_K_", K,"_coverage.rds")

saveRDS(results_cover, filename_coverage)







