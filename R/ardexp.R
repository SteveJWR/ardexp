
library(igraph)
library(ClusterR)
library(caret)
library(Matrix)
library(sandwich)
#library(mltools)
#library(data.table)



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
    groups = rep(seq(K),times = n.groups)
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

    Z.original = rep(seq(K),times = n.groups)
    perm = invPerm(order(Z))
    # marks the end of the groups

    # permute the appropriate Z
    G <- G[perm,perm]
    groups = Z

  }
  return(list('G' = G, "groups" = groups))
}


traitsSample <- function(Z,Q){
  n = length(Z)
  traits.vec = rep(NA, n)
  Tr = nrow(Q)
  for(i in seq(n)){
    trait = sample(seq(Tr), size = 1, prob = Q[,Z[i]])
    traits.vec[i] = trait
  }
  return(traits.vec)
}


computeARD <- function(traits.vec, G){
  n = nrow(G)
  Tr = max(traits.vec)
  X = matrix( nrow = n, ncol = Tr)

  for(k in seq(Tr)){
    ind = 1*(traits.vec == k)
    X[,k] = as.numeric(G %*% ind)
  }

  return(X)
}


# clusterARD <- function(X){
#   n = nrow(X)
#
#   W = 1:n
#   cluster.means = list()
#   while(length(W) > 0){
#     i = sample(W,size = 1)
#   }
#
# }

clusterARDKmeans <- function(X, K){
    res = stats::kmeans(X, centers = K)
    return(res)
}



estimatePmat <- function(Z.hat, G){
  K <- length(unique(Z.hat))
  n <- nrow(G)
  P <- matrix(NA, nrow = K, ncol = K)
  for(i in seq(K)){
    for(j in seq(K)){
      if(i <= j){
        idx.1 = which(Z.hat == i)
        idx.2 = which(Z.hat == j)
        if(i == j){
          G.sub = G[idx.1, idx.2]
          if(length(G.sub) == 1){
            p.hat = 0
          } else {
            n.sub <- nrow(G.sub)
            p.hat = sum(triu(G.sub))/choose(n.sub,2) # since it keeps the lower part.
          }

          P[i,j] = p.hat
        } else {
          G.sub = G[idx.1, idx.2]
          p.hat = mean(G.sub)
          P[i,j] = p.hat
        }
      }
    }
  }
  P[lower.tri(P)] = t(P)[lower.tri(t(P))]
  return(P)
}


estimatePmatARD <- function(Z.hat, X){
  K <- length(unique(Z.hat))
  n <- nrow(X)
  P <- matrix(NA, nrow = K, ncol = K)

  for(i in seq(K)){
    for(j in seq(K)){
      if(i <= j){
        idx.1 = which(Z.hat == i)

        X.sub = X[idx.1, j]
        if(i != j){
          n.total <- length(idx.1)*(sum(Z.hat == j))
        } else {
          n.total <- length(idx.1)*(sum(Z.hat == j) - 1)
        }

        p.hat = sum(X.sub)/n.total
        P[i,j] = p.hat

      }
    }
  }
  P[lower.tri(P)] = t(P)[lower.tri(t(P))]
  return(P)
}




linearOutcomeNormalSim <- function(Z,G,a,betasZ,  beta.a, beta.dtau,  noise = 1){
  d.tau = as.numeric(G %*% a) # number of treated neighbours
  K = max(Z)
  n = length(Z)
  Z.df = data.frame("z" = factor(Z, levels = seq(K)))
  dmy <- dummyVars(" ~ .", data = Z.df)
  Z.oh <- data.frame(predict(dmy, newdata = Z.df))
  Z.oh <- as.matrix(Z.oh)
  means = Z.oh %*% betasZ + beta.a *a + beta.dtau * d.tau
  Y = as.numeric(rnorm(n = n, mean = means, sd = noise))
  return(Y)
}



# linearOutcomeNormalSimByGroup <- function(Z,G,a,V, betasZ, betaV, noise = 1){
#   means
# }


treatmentAssignment <- function(Z,prob){
  K = max(Z)
  n = length(Z)
  a = rep(0,n)
  for(k in seq(K)){
    idx = which(Z == k)
    n.group = length(idx)
    n.treat = round(n.group*prob[k])
    idx.treat = sample(idx, size = n.treat, replace = F)
    a[idx.treat] = 1
  }
  return(a)
}


estimateQ <- function(traits, Z){
  Tr = max(traits)
  K = max(Z)
  Q = matrix(nrow = Tr, ncol = K)
  for(k in seq(K)){
    nc = sum(Z == k)
    for(tr in seq(Tr)){
      Q[tr,k] = sum(traits == tr & Z == k)/nc
    }
  }
  return(Q)
}

reverseConditionalQ <- function(traits, Z){
  Tr = max(traits)
  K = max(Z)
  Q = matrix(nrow = K, ncol = Tr)
  for(tr in seq(Tr)){
    nt = sum(traits == tr)
    for(k in seq(K)){
      Q[k,tr] = sum(traits == tr & Z == k)/nt
    }
  }
  return(Q)
}


# label swapping for when classification is good
labelSwitching <- function(Z,Z.hat){
  K = max(Z.hat)
  switch.model <- label.switching::ecr(Z.hat,matrix(Z, nrow = 1), K)
  return(switch.model$permutations[Z.hat])
}


linearOutcomeNormalizedNormalSim <- function(Z,G,a,betasZ,  beta.a, beta.dtau,  noise = 1){
  d.tau = as.numeric(G %*% a) # number of treated neighbours
  d.vec = colSums(G)

  K = max(Z)
  n = length(Z)
  Z.df = data.frame("z" = factor(Z, levels = seq(K)))
  dmy <- dummyVars(" ~ .", data = Z.df)
  Z.oh <- data.frame(predict(dmy, newdata = Z.df))
  Z.oh <- as.matrix(Z.oh)
  means = Z.oh %*% betasZ + beta.a *a + beta.dtau * (d.tau/d.vec)
  Y = as.numeric(rnorm(n = n, mean = means, sd = noise))
  return(Y)
}


linearOutcomeNormalizedByGroupNormalSim <- function(Z,G,a,betasZ,  beta.a, beta.dtau,  noise = 1){

  d.tau = as.numeric(G %*% a) # number of treated neighbours

  d.vec = colSums(G)

  K = max(Z)
  n = length(Z)
  Z.df = data.frame("z" = factor(Z, levels = seq(K)))
  dmy <- dummyVars(" ~ .", data = Z.df)
  Z.oh <- data.frame(predict(dmy, newdata = Z.df))
  Z.oh <- as.matrix(Z.oh)
  means = Z.oh %*% betasZ + beta.a *a + beta.dtau * (d.tau/d.vec)
  Y = as.numeric(rnorm(n = n, mean = means, sd = noise))
  return(Y)
}
