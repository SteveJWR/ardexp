

## Simuations to compare the full data and spillover parameters

rm(list = ls())
source("R/ardexp.R") # functions for the main
source("R/simulation_functions.R") # has the additional simulations pieces




slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
print(Sys.getenv('SLURM_ARRAY_TASK_ID'))
if(slurm_arrayid == ""){
  id = 1
} else {
  # coerce the value to an integer
  id <- as.numeric(slurm_arrayid)
}


if(id %% 8 == 1){
  mutual.benefit = T

} else if(id %% 8 == 2) {
  mutual.benefit = F

}

#village number
block = ceiling(id/2)

# True model parameters
a = 1
b = 1
delta = 1

# simulation noise
sigma = 0.5 # how much noise are we willing to permit here?



if(mutual.benefit){
  gamma0 = 1
  gamma1 = 0.2
} else {
  gamma0 = -1
  gamma1 = -0.2
}

data.file <- paste0("data/JPAL/Data/1. Network Data/Adjacency Matrices/adj_allVillageRelationships_vilno_", block, ".csv")

#data.file <- paste0("data/JPAL/Data/1. Network Data/Adjacency Matrices/adj_andRelationships_vilno_", block, ".csv")

#data.file.lend <- paste0("data/JPAL/Data/1. Network Data/Adjacency Matrices/adj_lendmoney_vilno_", block, ".csv")
#data.file.borrow <- paste0("data/JPAL/Data/1. Network Data/Adjacency Matrices/adj_borrowmoney_vilno_", block, ".csv")
#data.file.keroricecome <- paste0("data/JPAL/Data/1. Network Data/Adjacency Matrices/adj_keroricecome_vilno_", block, ".csv")
#data.file.keroricego <- paste0("data/JPAL/Data/1. Network Data/Adjacency Matrices/adj_keroricego_vilno_", block, ".csv")
#data.file<- paste0("data/JPAL/Data/1. Network Data/Adjacency Matrices/adj_nonrel_vilno_", block, ".csv")




# G.lend <- read.csv(data.file.lend, header = FALSE)
# G.lend <- as.matrix(G.lend)
# colnames(G.lend) = NULL
# rownames(G.lend) = NULL
#
# G.borrow <- read.csv(data.file.borrow, header = FALSE)
# G.borrow <- as.matrix(G.borrow)
# colnames(G.borrow) = NULL
# rownames(G.borrow) = NULL
#
# G.keroricecome <- read.csv(data.file.keroricecome, header = FALSE)
# G.keroricecome <- as.matrix(G.keroricecome)
# colnames(G.keroricecome) = NULL
# rownames(G.keroricecome) = NULL
#
#
# G.keroricego <- read.csv(data.file.keroricego, header = FALSE)
# G.keroricego <- as.matrix(G.keroricego)
# colnames(G.keroricego) = NULL
# rownames(G.keroricego) = NULL

# G <- (G.lend + G.borrow)
#



G <- read.csv(data.file, header = FALSE)
G <- as.matrix(G)
colnames(G) = NULL
rownames(G) = NULL
diag(G) = 0 #TODOLATER: Ask Tyler About this, why there are any self edges
G <- G + t(G)
G[G > 0] = 1

G.true <- G


g <- graph_from_adjacency_matrix(G.true, mode = "undirected")

clust <- igraph::spectrum(g)

plot(g, vertex.size=3, vertex.label=NA)

# remove dangling individuals
idx = which(colSums(G.true) > 0)
G.true = G.true[idx,idx]
g <- graph_from_adjacency_matrix(G.true, mode = "undirected")

# Off-the-shelf community detection algorithm
clust = cluster_fast_greedy(g)

# modularity measure
modularity(clust)
## [1] 0.3934089
# plot communities with shaded regions
coords = layout_with_fr(g)
plot(clust, g, layout=coords, vertex.size=3, vertex.label=NA)


estimate.sbm <- T
if(estimate.sbm){
  #sbm.true <- sbm::estimateSimpleSBM(G)

  # all we need are the "true clusters"
  Z.true <- clust$membership
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


# treatment probability for graph cluster designs
#p.treat = 0.5 #TODO: remove this



#p.in.seq = rep(0.2, J)
#p.out.seq = rep(0.02, J)

# other simulation parameters
# binomial.randomization <- F #TODO: remove this

# In order to have some degree heterogeneity, we say that the clusters have an
# equal spacing of parameters.


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


  #data.mean.sub <- data
  #data.mean.sub$frac.treated <- 0.2630463*(data.mean.sub$frac.treated < 0.5) + 0.7369646*(data.mean.sub$frac.treated > 0.5)
  #cor(data.mean.sub)
  #model.mean.sub <- lm(fmla, data.mean.sub)

  #model.mean.sub$coefficients
  #TODO: make the HT estimator even better

  res.vec <- c(ard.theta.est - theta,
               ard.true.model.theta.est - theta,
               theta.est - theta)
  names(res.vec) = NULL
  results[sim,] <- res.vec
}

colnames(results) = c("ard",
                      "ard.tm",
                      "reg")

filename <- paste0("data/JPAL_village_",block,
                   "mutual_benefit_", mutual.benefit,".rds")

saveRDS(results, filename)




library(reshape2)
library(ggplot2)

longData<-melt(P.true)
longData<-melt(P.hat - P.true)
longData<-longData[longData$value!=0,]

ggplot(longData, aes(x = Var2, y = Var1)) +
  geom_raster(aes(fill=value)) +
  scale_fill_gradient(low="grey90", high="red") +
  labs(x="letters", y="LETTERS", title="Matrix") +
  theme_bw() + theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
                     axis.text.y=element_text(size=9),
                     plot.title=element_text(size=11))



# figure out this plots thing next
make.plots = F
if(make.plots){
  for(id in seq(2)){
    if(id %% 8 == 1){
      mutual.benefit = T

    } else if(id %% 8 == 2) {
      mutual.benefit = F

    }
    library(ggpubr)
    library(abind)
    png.width = 1200
    png.height = 1000
    png.res = 200

    filename <- paste0("data/JPAL_village_",1,
                       "mutual_benefit_", mutual.benefit,".rds")
    results = readRDS(filename)
    for(i in seq(2,72)){
      filename <- paste0("data/JPAL_village_",i,
                         "mutual_benefit_", mutual.benefit,".rds")
      res.tmp <-readRDS(filename)
      results <- abind(results, res.tmp, along = 1)
    }
    n.sims = dim(results)[1]

    methods <- c("ard",
                 "ard.tm",
                 "reg")

    J = length(methods)
    n.seq <- c(100,316,1000,3162, 10000)
    sample.size.vec <- rep(n.seq, each = J)
    N = length(n.seq)

    bias.vec = c()
    rmse.vec = c()
    for(j in seq(N)){

      res.tmp = as.matrix(results[,,j])
      bias.vec = c(bias.vec, colMeans(res.tmp, na.rm = T))
      rmse.vec = c(rmse.vec, colMeans(abs(res.tmp)^2, na.rm = T)) # change to the mean absolute deviation
    }
    rmse.vec <- sqrt(rmse.vec)

    res.data <- data.frame("SampleSize" = sample.size.vec,
                           "Method" = methods,
                           "Bias" = bias.vec,
                           "RMSE" = rmse.vec)
    method.subset =  c("ard",
                       "ard.tm",
                       "reg")
    res.data <- res.data[res.data$Method %in% method.subset, ]

    plt.bias <- ggplot(res.data, aes(x = log(SampleSize), y = Bias, group = Method,color = Method)) +
      geom_line() +
      #geom_point() +
      #geom_errorbar(aes(ymin = ModelDev - 2*ModelDev_sd, ymax = ModelDev + 2*ModelDev_sd)) +
      ggtitle("RMSE of Methods") +
      xlab("log-Sample Size") +
      ylab("Bias")
    #geom_errorbar(aes(ymin = lik.mean.scaled - 2*lik.sd.scaled, ymax = lik.mean.scaled + 2*lik.sd.scaled))

    plt.bias

    plt.rmse <- ggplot(res.data, aes(x = log(SampleSize), y = RMSE, group = Method,color = Method)) +
      geom_line() +
      #geom_point() +
      #geom_errorbar(aes(ymin = ModelDev - 2*ModelDev_sd, ymax = ModelDev + 2*ModelDev_sd)) +
      ggtitle("RMSE of Methods") +
      xlab("log-Sample Size") +
      ylab("RMSE") +
      coord_cartesian(
        xlim =c(min(log(sample.size.vec)),max(log(sample.size.vec))),
        ylim = c(0,1)
      )
    plt.rmse
    #geom_errorbar(aes(ymin = lik.mean.scaled - 2*lik.sd.scaled, ymax = lik.mean.scaled + 2*lik.sd.scaled))

    #plt.rmse
    file.bias <- paste0("plot/Bias_Spillover_cluster_growth",cluster.growth,
                        "mutual_benefit_", mutual.benefit,
                        "cluster_equal_",cluster.equal.size,".png")
    file.rmse <- paste0("plot/RMSE_Spillover_cluster_growth",cluster.growth,
                        "mutual_benefit_", mutual.benefit,
                        "cluster_equal_",cluster.equal.size,".png")

    ggsave(
      filename = file.bias,
      plot = plt.bias,
      scale = 1,
      width = png.width,
      height = png.height,
      units = "px",
      dpi = png.res
    )

    ggsave(
      filename = file.rmse,
      plot = plt.rmse,
      scale = 1,
      width = png.width,
      height = png.height,
      units = "px",
      dpi = png.res
    )
  }
}







