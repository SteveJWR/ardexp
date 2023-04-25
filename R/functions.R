
library(igraph)
K = 10
P = matrix(data = rep(0.1,K^2), nrow = K, ncol = K)
diag(P) = 0.4
pi = rep(1/K,K)
# generates an SBM with fixed proportions pi.
generateSBM <- function(n,P,pi){
  if(nrow(P) != ncol(P)){
    stop("P must be symmetric")
  }
  if(P != t(P)){
    stop("P must be symmetric")
  }
  if(length(pi) != P){
    stop("proportions and P must be of the same dimension")
  }
  n.groups = round(n*pi/(sum(pi)))
  K = length(n.groups)
  Z = rep(seq(K),times = n.groups)
  g = igraph::sample_sbm(n,pref.matrix = P, block.sizes = n.groups)
  G = igraph::as_adjacency_matrix(g)
  return(list('G' = G, "groups" = z))
}

# Q = P(T = t|C = c)
Tr = 12 # number of traits
input.q = c(1,rep(0,K))
Q = matrix(rep(input.q, Tr)[1:(Tr*K)], nrow = Tr, ncol = K, byrow = T)
Q = t(t(Q)/colSums(Q))
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

#library(caret)
library(ClusterR)
cluster.res = external_validation(res$cluster,
                                  Z,
                                  method = "adjusted_rand_index",
                                  summary_stats = T)




linearOutcomeNormalSim <- function(Z,G,a,V, betasZ, betaV, noise = 1){
  means = betasZ %*%
}

