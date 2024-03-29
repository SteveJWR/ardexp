
# Here we have examples of the various outcome models for each of the models.



### replication of the simple Ugander response model
simUganderModelOutcomes <- function(G,H,A, a = 1, b = 2, delta = 1, gamma = 1, sigma = 1) {
  n = nrow(G)
  d.vec = as.numeric(G %*% rep(1,n))
  d.mean = mean(d.vec)

  eps <- rnorm(n)
  nonoise.Y0 <- a + b*H
  Y0 <- (a + b*H + sigma*eps)*(d.vec/d.mean)

  num.treated = as.numeric(G %*% A)
  f = fracTreated(G,A)
  f[is.na(f)] = 0 # no neighbors means none are treated
  YA <-Y0*(1 + delta*A + gamma*f)
  gate = mean(Y0)*(delta + gamma)
  return(list("Y" = YA, "GATE" = gate))
}


simSpilloverModelOutcomes <- function(G,A, a = 1, b = 0, delta = 1, gamma0 = 1, gamma1 = 0.2, sigma = 1) {
  n = nrow(G)
  d.vec = as.numeric(G %*% rep(1,n))
  d.mean = mean(d.vec)

  eps <- rnorm(n)
  Y0 <- (a  + b*d.vec/d.mean + sigma*eps)

  num.treated = as.numeric(G %*% A)
  f = fracTreated(G,A)
  f[is.na(f)] = 0 # no neighbors means none are treated

  YA <- Y0 + delta*A + gamma0*(1 - A)*f + gamma1*(A)*f
  return(list("Y" = YA))
}


# true only
fracTreated <- function(G, A){
  num.treated = as.numeric(G %*% A)
  d.vec = colSums(G)
  return(num.treated/d.vec)
}


# if we want to do graph cluster randomization, we can keep things simple with the natural clusters of the graph
oneHopClustering <- function(G,w = rep(nrow(G),nrow(G)), verbose = F){
  n = nrow(G)
  X = rbeta(n,w,1)
  z.vec <- rep(NA,n)
  for(i in seq(n)){
    if(i %% 1000 == 0 & verbose){
      cat(paste0("Cluster: ", i,"/",n),end = "\r")
    }
    zi = max(X[G[i,] == 1])
    z.vec[i] = zi
  }
}


simH <- function(z.vec, sigma.within = 0.5, sigma.between = 2){
  unique.clusters <- unique(z.vec)
  K <- length(unique.clusters)
  H.block <- rnorm(K, sd = sigma.between)
  H = rep(NA, length(z.vec))
  for(i in seq(K)){
    z.true = unique.clusters[i]
    z.idx = which(z.vec == z.true)
    n.block = length(z.idx)
    H[z.idx] = rnorm(n.block, mean = H.block[i], sd = sigma.within)
  }
  return(H)
}

binomialTreatment <- function(n,p){
  A <- rbinom(n,size = 1, prob = p)
  return(A)
}

clusterTreatment <- function(Z, p.treat = 1/2){
  unique.clusters <- unique(Z)
  K <- length(unique.clusters)
  K.treat <- round(K*p.treat)
  treatments <- c(rep(1,K.treat),rep(0,K - K.treat))
  A.block <- sample(treatments, replace = F)
  A = rep(NA, length(Z))
  for(i in seq(K)){
    z.true = unique.clusters[i]
    z.idx = which(Z == z.true)
    A[z.idx] = A.block[i]
  }
  return(A)
}

saturationRandomizationTreatment <- function(Z, levels){
  K = max(Z)
  # Z must be an integer between 1 and K
  if(length(levels) != K){
    stop("saturation levels must be the same as the number of clusters")
  }
  n = length(Z)
  A = rep(0,n)
  for(i in seq(K)){
   idx.i = which(Z == i)
   ni = length(idx.i)
   ni.treat = round(ni*levels[i])
   idx.i.treat <- sample(idx.i, ni.treat, replace = F)
   A[idx.i.treat] = 1
  }
  return(A)
}








neighborVectors <- function(G){
  nbr.list <- list()
  g = igraph::graph_from_adjacency_matrix(G, mode = "undirected")
  n = nrow(G)
  nbs <- igraph::neighborhood(g)
  return(nbs)
}

#### Competing methods

# frac is the fraction of blocks treated
HTEstimatorCluster <- function(Y,G,A,Z,p.treat = 1/2){
  n = nrow(G)
  pE0 <- rep(NA, n)
  pE1 <- rep(NA, n)
  IE0 <- rep(0, n)
  IE1 <- rep(0, n)
  nbrs <- neighborVectors(G)
  K = length(unique(Z))
  K.treat = round(K*p.treat)
  # rate limiting step
  neighbors <- neighborVectors(G)

  for(i in seq(n)){
    #cat(paste0("Node: ", i, "/",n), end = "\r")
    neighborhood <- c(i,neighbors[[i]])

    if(all(A[neighborhood] == 1)){
      IE1[i] = 1
      num.blocks <- length(unique(Z[neighborhood]))
      pE1[i] = choose(K.treat,num.blocks)/choose(K,num.blocks)
    }
    if(all(A[neighborhood] == 0)){
      IE0[i] = 1
      num.blocks <- length(unique(Z[neighborhood]))
      pE0[i] = choose(K - K.treat,num.blocks)/choose(K,num.blocks)
    }
  }

  out.1 <- (1/n)*Y*IE1/pE1
  out.2 <- (1/n)*Y*IE0/pE0

  not.na.idx1 <- which(!is.na(out.1))
  not.na.idx2 <- which(!is.na(out.2))

  if(length(not.na.idx1) > 0 & length(not.na.idx1) > 0){
    HT = sum(out.1[not.na.idx1]) - sum(out.2[not.na.idx2])
  } else {
    HT = NA
  }
  return(HT)
}


# HTEstimatorBernoulli <- function(Y,G,A,Z,p = 1/2){
#   n = nrow(G)
#   pE0 <- rep(NA, n)
#   pE1 <- rep(NA, n)
#   IE0 <- rep(0, n)
#   IE1 <- rep(0, n)
#   nbrs <- neighborVectors(G)
#   K = length(unique(Z))
#   K.treat = round(K*frac)
#   neighbors <- neighborVectors(G)
#
#   for(i in seq(n)){
#     neighborhood <- c(i,neighbors[[i]])
#
#     if(all(A[neighborhood] == 1)){
#       IE1[i] = 1
#       num.blocks <- length(neighborhood)
#       pE1[i] = p^(num.blocks)
#     }
#     if(all(A[neighborhood] == 0)){
#       IE0[i] = 1
#       num.blocks <- length(neighborhood)
#       pE0[i] = (1 - p)^(num.blocks)
#     }
#   }
#
#   out.1 <- (1/n)*Y*IE1/pE1
#   out.2 <- (1/n)*Y*IE0/pE0
#
#
#   not.na.idx1 <- which(!is.na(out.1))
#   not.na.idx2 <- which(!is.na(out.2))
#
#   if(length(not.na.idx1) > 0 & length(not.na.idx1) > 0){
#     HT = sum(out.1[not.na.idx1]) - sum(out.2[not.na.idx2])
#   } else {
#     HT = NA
#   }
#   return(HT)
# }


diffMeansEstimator <- function(Y,A){
  tau = mean(Y[A == 1]) - mean(Y[A == 0])
  return(tau)
}



meanOverList <- function(lst){
  obj <- lst[[1]]
  B = length(lst)
  for(b in seq(2,B)){
    obj <- obj + lst[[b]]
  }
  return((1/B)*obj)
}

concatenateOverList <- function(lst){
  obj <- lst[[1]]
  B = length(lst)
  for(b in seq(2,B)){
    obj <- cbind(obj,lst[[b]])
  }
  return(obj)
}

varOverList <- function(lst){
  mean.obj <- meanOverList(lst)
  B = length(lst)
  obj <- (lst[[1]] - mean.obj)^2
  for(b in seq(2,B)){
    obj <- obj + (lst[[b]] - mean.obj)^2
  }
  return(1/(B - 1)*obj)
}

# no intercept must be suggested for the linear model
# or otherwise we let it just be 0
#


# The user should be able to define their own functions
# H is a variable for covariates here.
UganderCovariates <- function(G,H,A){
  n = nrow(G)
  d.vec = as.numeric(G %*% rep(1,n))
  d.mean = mean(d.vec)
  f <- fracTreated(G,A)
  cov.data <- data.frame("A" = A,
                         "frac.treated" = f,
                         "H" = H,
                         "deg.ratio" = d.vec/d.mean)
  return(cov.data)
}

SpilloverCovariates <- function(G,A){
  n = nrow(G)
  d.vec = as.numeric(G %*% rep(1,n))
  d.mean = mean(d.vec)
  f <- fracTreated(G,A)

  # Exception for the unconnected component
  isolated_points = d.vec == 0
  f[isolated_points] = 0

  cov.data <- data.frame("A" = A,
                         "frac.treated" = f,
                         "deg.ratio" = d.vec/d.mean)
  return(cov.data)
}


LinearDiffusionFeatures <- function(G,A,T_steps = 5){
  n = nrow(G)
  cov.mat = matrix(ncol = T_steps + 1, nrow = n)
  cov.mat[,1] = A
  G_power = diag(1, n)
  for(t_step in seq(T_steps)){
    G_power = G_power %*% G
    t_neighbors = G_power %*% A
    cov.mat[,1 + t_step] = t_neighbors
  }
  cov.data <- as.data.frame(cov.mat)
  names(cov.data) = c('A', paste('treated_walk_',as.character(seq(T_steps)), sep=""))
  return(cov.data)
}


# H is a set of covariates
# A are the treatments
# ,... can pass other parameters into the linear model
# This is just a parametric Bootstrap
ARDUganderSBMLinearRegressionSim <- function(Y,fmla, A, H, P, Z, B.boot = 200, verbose = F){
  n = length(Y)
  A.0 <- rep(0,n)
  A.1 <- rep(1,n)
  coef.list <- list()
  var.list <- list()
  data.list <- list()
  gate.list <- list()
  for(b in seq(B.boot)){
    if (verbose) {
      m1 = (round(20 * b/B.boot))
      m2 = 20 - m1
      progress.bar = paste0("|", strrep("=", m1), strrep("-",
                                                         m2), "|")
      cat(paste0("Resample:", b, "/", B.boot, "  ", progress.bar),
          end = "\r")
    }
    g.sim <- generateSBM(n,P,Z = Z)
    G = g.sim$G
    data <- UganderCovariates(G,H,A)
    data$Y <- Y

    data.0 <- UganderCovariates(G,H,A.0)
    data.1 <- UganderCovariates(G,H,A.1)

    ## manual fix
    data.0$frac.treated = 0
    data.1$frac.treated = 1


    model.sim <- lm(fmla, data = data)
    Y.0 <- predict(model.sim, data.0)
    Y.1 <- predict(model.sim, data.1)
    gate = mean(Y.1) - mean(Y.0)
    coef.list[[b]] <- model.sim$coefficients
    var.list[[b]] <- sandwich::sandwich(model.sim)
    data.list[[b]] <- data
    gate.list[[b]] <- gate
  }

  return(list("coef" = coef.list, "var" = var.list, "data" = data.list, "gate" = gate.list))
}


ARDDiffusionRegressionSim <- function(A, P, Z, T.max = 3, B.boot = 200, verbose = F){

  coef.list <- list()
  var.list <- list()
  data.list <- list()
  gate.list <- list()
  for(b in seq(B.boot)){
    if (verbose) {
      m1 = (round(20 * b/B.boot))
      m2 = 20 - m1
      progress.bar = paste0("|", strrep("=", m1), strrep("-",
                                                         m2), "|")
      cat(paste0("Resample:", b, "/", B.boot, "  ", progress.bar),
          end = "\r")
    }
    g.sim <- generateSBM(n,P,Z = Z)
    G = g.sim$G
    data <- DiffusionExampleCovariates(G,A, T.max = T.max)
    data.list[[b]] <- data
  }

  return(list( "data" = data.list))
}

# Here we use a simple linear model for maximizing over the treatment assignments.
ExampleLinearModelCovariates <- function(G,A) {
  n = nrow(G)
  d.vec = as.numeric(G %*% rep(1,n))
  f <- fracTreated(G,A)
  G2 = G %*% G
  f2 = fracTreated(G2,A)

  f[is.na(f)] = 0 # no neighbors means none are treated
  f2[is.na(f2)] = 0 # no neighbors means none are treated
  cov.data <- data.frame("A" = A,
                         "frac.treated" = f,
                         "frac.treated.A0" = f*(1 - A),
                         "frac.treated.second" = f2,
                         "degree" = d.vec)
}

simExampleLinearModelOutcomes <- function(G,A,a = 1,delta = 1,gamma = 0.5, gamma2 = -0.5, sigma = 1) {
  n = nrow(G)
  eps <- rnorm(n)
  f = fracTreated(G,A)
  G2 = G %*% G
  f2 = fracTreated(G2,A)

  f[is.na(f)] = 0 # no neighbors means none are treated
  f2[is.na(f2)] = 0 # no neighbors means none are treated
  Y <- a + delta*A + gamma*f + gamma2*f2 + sigma*eps
  meanY <- a + delta*A + gamma*f + gamma2*f2
  return(list("Y" = Y, "meanY" = meanY))
}

ARDExampleLinearModelSim <- function(Y,fmla, A, P, Z, B.boot = 200, verbose = F){
  n = length(Y)
  A.0 <- rep(0,n)
  A.1 <- rep(1,n)
  coef.list <- list()
  var.list <- list()
  data.list <- list()
  gate.list <- list()
  for(b in seq(B.boot)){
    if (verbose) {
      m1 = (round(20 * b/B.boot))
      m2 = 20 - m1
      progress.bar = paste0("|", strrep("=", m1), strrep("-",
                                                         m2), "|")
      cat(paste0("Resample:", b, "/", B.boot, "  ", progress.bar),
          end = "\r")
    }
    g.sim <- generateSBM(n,P,Z = Z)
    G = g.sim$G
    data <- ExampleLinearModelCovariates(G,A)
    data$Y <- Y # adding the outcome


    model.sim <- lm(fmla, data = data)
    coef.list[[b]] <- model.sim$coefficients
    var.list[[b]] <- sandwich::sandwich(model.sim) # Sandwich variances for the examples
    data.list[[b]] <- data
  }
  return(list("coef" = coef.list, "var" = var.list, "data" = data.list))
}



ARDSBMSpilloverLinearRegressionSim <- function(Y,fmla, graphMapping, A, P, Z, B.boot = 200, verbose = F){
  n = length(Y)
  A.0 <- rep(0,n)
  A.1 <- rep(1,n)
  coef.list <- list()

  var.list <- list()
  data.list <- list()
  gate.list <- list()
  for(b in seq(B.boot)){
    if (verbose) {
      m1 = (round(20 * b/B.boot))
      m2 = 20 - m1
      progress.bar = paste0("|", strrep("=", m1), strrep("-",
                                                         m2), "|")
      cat(paste0("Resample:", b, "/", B.boot, "  ", progress.bar),
          end = "\r")
    }
    g.sim <- generateSBM(n,P,Z = Z)
    G = g.sim$G
    data <- graphMapping(G,A)
    data$Y <- Y

    ## manual fix

    model.sim <- lm(fmla, data = data)
    coef.list[[b]] <- model.sim$coefficients
    var.list[[b]] <- sandwich::sandwich(model.sim)
    data.list[[b]] <- data
  }

  return(list("coef" = coef.list, "var" = var.list, "data" = data.list))
}

# only randomizes over Z blocks
optimalDesignSpilloverSBMGreedy <- function(phi,graphMapping,P,Z,grid.steps = 10, B.boot = 200, verbose = F, max.cycles = 10){
  K = length(unique(Z))
  tau.grid = seq(0,1,length.out = grid.steps)
  tau.list =list()
  G.list <- list()
  for(b in seq(B.boot)){
    if (verbose) {
      m1 = (round(20 * b/B.boot))
      m2 = 20 - m1
      progress.bar = paste0("|", strrep("=", m1), strrep("-",
                                                         m2), "|")
      cat(paste0("Resample:", b, "/", B.boot, "  ", progress.bar),
          end = "\r")
    }
    g.sim <- generateSBM(n,P,Z = Z)
    G = g.sim$G
  }

  for(j in seq(grid.steps)){
    tau.list[[j]] = tau.grid
  }
  tau.levels <- expand.grid(tau.list)
}


# Helper function to clean the graph data
clean_diffusion_graphs <- function(graphs){
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
  colnames(G) = NULL
  rownames(G) = NULL
  diag(G) = 0 # Remove self edges from the noisy sample
  G <- G + t(G)
  G[G > 0] = 1
  G.true <- G
  idx = which(colSums(G.true) > 0)
  G.true = G.true[idx,idx]
  return(G.true)
}



assign_balls_to_bins <- function(N, bin_sizes, decreasing = TRUE) {
  # Sort bin sizes from largest to smallest
  sorted_bin_sizes <- sort(bin_sizes, decreasing = decreasing)

  bin_order = order(bin_sizes, decreasing = decreasing)
  # Initialize vector to store number of balls in each bin
  balls_in_bins <- rep(0, length(bin_sizes))

  remaining_balls = N
  k = 1
  while(remaining_balls > 0){
    bin_index = bin_order[k]
    open_space = bin_sizes[bin_index]
    if(open_space <= remaining_balls) {
      balls_in_bins[bin_index] = bin_sizes[bin_index] # we fill the bin
    } else {
      balls_in_bins[bin_index] = remaining_balls # we place the remaining balls
    }
    remaining_balls = remaining_balls - balls_in_bins[bin_index]
    k = k + 1
  }
  return(balls_in_bins)
}



cluster_leiden_K <- function(G, K = 4, adj_mode = "undirected", max_resolution_parameter = 10, n_iterations = 5){
  g <- graph_from_adjacency_matrix(G, mode = adj_mode)
  tmp_cluster = cluster_leiden(g, objective_function = "modularity",
                               n_iterations = n_iterations, resolution_parameter = max_resolution_parameter)
  resolution_parameter_upper = max_resolution_parameter
  resolution_parameter_lower = 0
  iter = 0
  if(tmp_cluster$nb_clusters == K){
    return(tmp$membership) # we are done
  }
  if(tmp_cluster$nb_clusters <  K) {
    stop("Not enough clusters")
  } else {
    while(tmp_cluster$nb_clusters != K){
      iter = iter + 1
      resolution_parameter_tmp = (resolution_parameter_upper + resolution_parameter_lower)/2
      tmp_cluster = cluster_leiden(g, objective_function = "modularity",
                                   n_iterations = n_iterations, resolution_parameter = resolution_parameter_tmp)
      if(tmp_cluster$nb_clusters > K){
        resolution_parameter_upper = resolution_parameter_tmp
      } else {
        resolution_parameter_lower = resolution_parameter_tmp
      }
      if( iter > 100){
        stop("Too many iterations")
      }

    }
    return(tmp_cluster$membership)
  }
}


