

### replication of the simple Ugander response model
# generates an SBM with fixed proportions PI.
# or if there is no PI, then use Z
generateSBM <- function(n,P,PI,Z){
  if(missing(Z)){
    if(nrow(P) != ncol(P)){
      stop("P must be symmetric")
    }
    if(sum(abs(P  - t(P))) > 0.01){
      stop("P must be symmetric")
    }
    if(length(PI) != nrow(P)){
      stop("proportions and P must be of the same dimension")
    }
    n.groups = round(n*PI/(sum(PI)))
    K = length(PI)
    while(sum(n.groups) < n){
      i = sample(seq(K), size = 1)
      n.groups[i] = n.groups[i] + 1
    }
    while(sum(n.groups) > n){
      i = sample(seq(K), size = 1)
      n.groups[i] = n.groups[i] - 1
    }
    K = length(PI)
    A = rep(seq(K),times = n.groups)
    g = igraph::sample_sbm(n,pref.matrix = P, block.sizes = n.groups)
    G = igraph::as_adjacency_matrix(g)
  } else {
    K = max(Z)
    n.groups <- rep(NA,K)
    for(i in seq(K)){
      n.groups[i] = sum(Z == i)
    }

    # permute the ordering at the end so that it will agree with Z
    g = igraph::sample_sbm(n,pref.matrix = P, block.sizes = n.groups)
    G = igraph::as_adjacency_matrix(g)

    # marks the end of the groups
    group.end.idx <- cumsum(n.groups)

    rev.perm <- rep(NA,n)
    for(i in seq(n)){
      rev.perm[i] = group.end.idx[Z[i]]
      group.end.idx[Z[i]] = group.end.idx[Z[i]] - 1
    }
    perm = rev.perm[invPerm(rev.perm)]

    # permute the appropriate Z
    G <- G[perm,perm]
    A = Z

  }
  return(list('G' = G, "groups" = A))
}

simUganderModelOutcomes <- function(G,H,A, a = 1, b = 2, delta = 1, gamma = 1, sigma = 1, scale.deg = T) {
  n = nrow(G)
  d.vec = as.numeric(G %*% rep(1,n))
  d.mean = mean(d.vec)

  eps <- rnorm(n)
  nonoise.Y0 <- a + b*H
  if(scale.deg){
    Y0 <- (a + b*H + sigma*eps)*(d.vec/d.mean)
  } else {
    Y0 <- (a + b*H + sigma*eps)
  }

  num.treated = as.numeric(G %*% A)
  f = fracTreated(G,A)
  f[is.na(f)] = 0 # no neighbors means none are treated
  YA <-Y0*(1 + delta*A + gamma*f)
  gate = mean(nonoise.Y0)*(delta + gamma)
  return(list("Y" = YA, "GATE" = gate))
}

simSpilloverModelOutcomes <- function(G,A, a = 1, b = 1, delta = 1, gamma0 = 1, gamma1 = 0.2, sigma = 1) {
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

# z.vec is the cluster vector
simH <- function(z.vec, sigma.within = 1, sigma.between = 5){
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

clusterTreatment <- function(z.vec, frac = 1/2){
  unique.clusters <- unique(z.vec)
  K <- length(unique.clusters)
  K.treat <- round(K*frac)
  treatments <- c(rep(1,K.treat),rep(0,K - K.treat))
  A.block <- sample(treatments, replace = F)
  A = rep(NA, length(z.vec))
  for(i in seq(K)){
    z.true = unique.clusters[i]
    z.idx = which(z.vec == z.true)
    A[z.idx] = A.block[i]
  }
  return(A)
}

saturationRandomizationTreatment <- function(c.vec, levels){
  K = max(c.vec)
  # c.vec must be an integer between 1 and K
  if(length(levels) != K){
    stop("saturation levels must be the same as the number of clusters")
  }
  n = length(c.vec)
  A = rep(0,n)
  for(i in seq(K)){
   idx.i = which(c.vec == i)
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

#TODO: think about where the c.vec should be denoted and where it should be z for a true clustering
# frac is the fraction of blocks treated
HTEstimatorCluster <- function(Y,G,A,c.vec,p.treat = 1/2){
  n = nrow(G)
  pE0 <- rep(NA, n)
  pE1 <- rep(NA, n)
  IE0 <- rep(0, n)
  IE1 <- rep(0, n)
  nbrs <- neighborVectors(G)
  K = length(unique(c.vec))
  K.treat = round(K*p.treat)
  # rate limiting step
  neighbors <- neighborVectors(G)

  for(i in seq(n)){
    #cat(paste0("Node: ", i, "/",n), end = "\r")
    neighborhood <- c(i,neighbors[[i]])

    if(all(A[neighborhood] == 1)){
      IE1[i] = 1
      num.blocks <- length(unique(c.vec[neighborhood]))
      pE1[i] = choose(K.treat,num.blocks)/choose(K,num.blocks)
    }
    if(all(A[neighborhood] == 0)){
      IE0[i] = 1
      num.blocks <- length(unique(c.vec[neighborhood]))
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


HTEstimatorBernoulli <- function(Y,G,A,p = 1/2){
  n = nrow(G)
  pE0 <- rep(NA, n)
  pE1 <- rep(NA, n)
  IE0 <- rep(0, n)
  IE1 <- rep(0, n)
  nbrs <- neighborVectors(G)
  K = length(unique(c.vec))
  K.treat = round(K*frac)
  neighbors <- neighborVectors(G)

  for(i in seq(n)){
    neighborhood <- c(i,neighbors[[i]])

    if(all(A[neighborhood] == 1)){
      IE1[i] = 1
      num.blocks <- length(neighborhood)
      pE1[i] = p^(num.blocks)
    }
    if(all(A[neighborhood] == 0)){
      IE0[i] = 1
      num.blocks <- length(neighborhood)
      pE0[i] = (1 - p)^(num.blocks)
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
# H is a good variable for covariates here.
UganderCovariates <- function(G,H,A){
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
  d.vec = as.numeric(G %*% rep(1,n))
  d.mean = mean(d.vec)
  f <- fracTreated(G,A)
  cov.data <- data.frame("A" = A,
                         "frac.treated" = f,
                         "deg.ratio" = d.vec/d.mean)
  return(cov.data)
}

# H is a set of covariates
# A are the treatments
# ,... can pass other parameters into the linear model
# This is just a parametric Bootstrap
ARDSBMLinearRegressionSim <- function(Y,fmla, graphMapping, A, H, P, Z, B.boot = 200, verbose = F){
  n = length(Y)
  A.0 <- rep(0,n)
  A.1 <- rep(1,n)
  coef.list <- list()
  # almost certainly we should use a robust variance here
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
    data <- graphMapping(G,H,A)
    data$Y <- Y

    data.0 <- graphMapping(G,H,A.0)
    data.1 <- graphMapping(G,H,A.1)
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


ARDSBMSpilloverLinearRegressionSim <- function(Y,fmla, graphMapping, A, P, Z, B.boot = 200, verbose = F){
  n = length(Y)
  A.0 <- rep(0,n)
  A.1 <- rep(1,n)
  coef.list <- list()
  # almost certainly we should use a robust variance here
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


optimalDesignSpilloverSBM <- function(phi,graphMapping,P,Z,C,grid.steps = 10, B.boot = 200, M = 100, verbose = F){
  length(unique())
}








