

## Simuations to compare the full data and spillover parameters

rm(list = ls())
source("R/ardexp.R") # functions for the main
source("R/simulation_functions.R") # has the additional simulations pieces

## competing methods DM, HT, Oracle Regression, Partial Data Regression, Partial Data Regression (Estimated Clusters)


slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
print(Sys.getenv('SLURM_ARRAY_TASK_ID'))
if(slurm_arrayid == ""){
  id = 2
} else {
  # coerce the value to an integer
  id <- as.numeric(slurm_arrayid)
}


if(id %% 8 == 1){
  mutual.benefit = T
  cluster.growth = T
  cluster.equal.size = T

} else if(id %% 8 == 2) {
  mutual.benefit = F
  cluster.growth = T
  cluster.equal.size = T

} else if(id %% 8 == 3) {
  mutual.benefit = T
  cluster.growth = F
  cluster.equal.size = T

} else if(id %% 8 == 4) {
  mutual.benefit = F
  cluster.growth = F
  cluster.equal.size = T

} else if(id %% 8 == 5) {
  mutual.benefit = T
  cluster.growth = T
  cluster.equal.size = F

} else if(id %% 8 == 6) {
  mutual.benefit = F
  cluster.growth = T
  cluster.equal.size = F

} else if(id %% 8 == 7) {
  mutual.benefit = T
  cluster.growth = F
  cluster.equal.size = F

} else if(id %% 8 == 0) {
  mutual.benefit = F
  cluster.growth = F
  cluster.equal.size = F

}

block = ceiling(id/8)

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


# Think of a second version of the simulation where there is a competative advantage to treatment
# but that wears off if everyone else around you has it

#Global average treatment effect
theta = gamma1 - gamma0

n.sims =  10#

n.seq <- c(100,316,1000,3162,10000) # log scale population growth.
#n.seq <- c(100,316,1000,3162) #TODO: Do this with the larger sample size.
J = length(n.seq) # number of sample sizes to look at

# treatment probability for graph cluster designs
p.treat = 0.5


if(cluster.growth){
  K.seq = round((n.seq)**(1/2)) # square root growth in number of clusters
} else {
  K.seq = rep(10,J)
}




p.in.seq = rep(0.2, J)
p.out.seq = rep(0.02, J)

# other simulation parameters
binomial.randomization <- F

# In order to have some degree heterogeneity, we say that the clusters have an
# equal spacing of parameters.


# we decouple the design piece and only talk about inference.
n.methods = 3
results <- array(NA, c(n.sims, n.methods, J))

# Methods tuning parameters
B.boot = 50
#### Simulations Start

# formula for the regression model
fmla <- formula(Y ~ deg.ratio + A*frac.treated )
#fmla <- formula(Y ~ (H + A + H:A + frac.treated + H:frac.treated))


# SaturationRandomization Design.
sat.frac = 1/2 #treat half of the clusters at 0.9 and the other half at 0.1
p.high.level = 0.9
p.low.level = 0.1


for(j in seq(J)){
  for(sim in seq(n.sims)){
    cat(paste0("Sim: ", sim, "/", n.sims, " --- Sample Size: ", j,"/", J), end = "\r")
    n = n.seq[j]
    K = K.seq[j]

    P = matrix(data = 0, nrow = K, ncol = K)

    p.in.max = p.in.seq[j]
    p.in.min = p.in.max #p.in.max/2 TODO: put this back


    p.out.max = p.out.seq[j]
    p.out.min = p.out.max #p.out.max/2 TODO: put this back


    P[upper.tri(P)] = seq(p.out.min,p.out.max, length.out = choose(K,2))
    P = P + t(P)
    diag(P) = seq(p.in.min,p.in.max, length.out = K)


    if(cluster.equal.size){
      PI = rep(1/K,K)
    } else {
      PI = seq(1,sqrt(K),K)
      PI = PI/sum(PI)
    }

    g.true <- generateSBM(n,P,PI)
    G.true <- g.true$G
    Z.true = g.true$groups

    #H <- simH(Z.true)
    H = rep(NA,n)

    # whether to do binomial randomization or the true cluster treatment
    A.clust <- clusterTreatment(Z.true, p.treat)



    K.high.level <- round(sat.frac*K)
    K.low.level <- K - K.high.level
    clust.levels <- c(rep(p.high.level, K.high.level), rep(p.low.level, K.low.level))
    A.sat <- saturationRandomizationTreatment(Z.true, levels = clust.levels)


    outcome.clust <- simSpilloverModelOutcomes(G.true,A.clust,a = a, b = b, delta = delta, gamma0 = gamma0, gamma1 = gamma1, sigma = sigma)
    outcome.sat <- simSpilloverModelOutcomes(G.true,A.sat,a = a, b = b, delta = delta, gamma0 = gamma0, gamma1 = gamma1, sigma = sigma)
    Y.clust <- outcome.clust$Y
    Y.sat <- outcome.sat$Y


    d.vec = as.numeric(G.true %*% rep(1,n))
    d.mean = mean(d.vec)
    f <- fracTreated(G.true,A.sat)
    data <- SpilloverCovariates(G.true,A.sat)
    data <- data.frame("Y" = Y.sat, "A" = A.sat, "frac.treated" = f, "deg.ratio" = d.vec/d.mean)


    model <- lm(fmla, data = data)



    # full data version.
    # (Note I had e)
    theta.est <- model$coefficients[5] - model$coefficients[4]


    # Here we are asking about the traits directly for simplicity in the simulations.

    traits <- g.true$groups

    X <-computeARD(traits, G.true)

    # Assuming For simplicity, K is known.
    k.means.model <- clusterARDKmeans(X,K)


    Z.hat = k.means.model$cluster
    #
    Z.hat = labelSwitching(Z.true,Z.hat)


    P.hat <- estimatePmat(Z.hat, G.true)


    res.ard <- ARDSBMSpilloverLinearRegressionSim(Y.sat, fmla, SpilloverCovariates, A.sat, P.hat,Z.hat, B.boot = B.boot, verbose = T)

    ard.avg.coef <- meanOverList(res.ard$coef)
    ard.avg.data <- meanOverList(res.ard$data)

    ard.theta.est <- ard.avg.coef[5] - ard.avg.coef[4]

    res.ard.true.model <- ARDSBMSpilloverLinearRegressionSim(Y.sat, fmla, SpilloverCovariates, A.sat, P, Z.true, B.boot = B.boot, verbose = T)
    ard.true.model.avg.coef <- meanOverList(res.ard.true.model$coef)

    ard.true.model.theta.est <- ard.true.model.avg.coef[5] - ard.true.model.avg.coef[4]


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
    results[sim,,j] <- res.vec
  }
}


colnames(results) = c("ard",
                      "ard.tm",
                      "reg")

filename <- paste0("data/SpilloverGATE_cluster_growth",cluster.growth,
                   "mutual_benefit_", mutual.benefit,
                   "cluster_equal_",cluster.equal.size,
                   "block",block,".rds")

saveRDS(results, filename)




make.plots = F
if(make.plots){
  library(ggpubr)
  library(abind)
  png.width = 1200
  png.height = 1000
  png.res = 200



  filename <- paste0("data/SpilloverGATE_cluster_growth",cluster.growth,
                     "mutual_benefit_", mutual.benefit,
                     "cluster_equal_",cluster.equal.size,
                     "block",1,".rds")
  results = readRDS(filename)
  for(i in seq(2,100)){
    filename <- paste0("data/SpilloverGATE_cluster_growth",cluster.growth,
                       "mutual_benefit_", mutual.benefit,
                       "cluster_equal_",cluster.equal.size,
                       "block",i,".rds")
    res.tmp <-readRDS(filename)
    results <- abind(results, res.tmp, along = 1)
  }
  n.sims = dim(results)[1]

  methods <- c("ard1", "ard2", "ard3",
               "ard.tm1", "ard.tm2", "ard.tm3",
               "reg1", "reg2", "reg3",
               "HT", "DM")

  J = length(methods)
  n.seq <- c(100,316,1000,3162)
  sample.size.vec <- rep(n.seq, each = J)
  N = length(n.seq)

  bias.vec = c()
  rmse.vec = c()
  for(j in seq(N)){

    res.tmp = as.matrix(results[,,j])
    bias.vec = c(bias.vec, colMeans(res.tmp, na.rm = T))
    rmse.vec = c(rmse.vec, colMeans(abs(res.tmp), na.rm = T)) # change to the mean absolute deviation
  }
  #rmse.vec <- sqrt(rmse.vec)

  res.data <- data.frame("SampleSize" = sample.size.vec,
                         "Method" = methods,
                         "Bias" = bias.vec,
                         "RMSE" = rmse.vec)
  method.subset =  c("ard1", "ard2", "ard3","reg1", "reg2", "reg3", "DM", "HT")
  res.data <- res.data[res.data$Method %in% method.subset, ]

  plt.bias <- ggplot(res.data, aes(x = log(SampleSize), y = Bias, group = Method,color = Method)) +
    geom_line() +
    #geom_point() +
    #geom_errorbar(aes(ymin = ModelDev - 2*ModelDev_sd, ymax = ModelDev + 2*ModelDev_sd)) +
    ggtitle("Bias of Methods") +
    xlab("log-Sample Size") +
    ylab("Bias")
  #geom_errorbar(aes(ymin = lik.mean.scaled - 2*lik.sd.scaled, ymax = lik.mean.scaled + 2*lik.sd.scaled))

  plt.bias

  plt.rmse <- ggplot(res.data, aes(x = log(SampleSize), y = RMSE, group = Method,color = Method)) +
    geom_line() +
    #geom_point() +
    #geom_errorbar(aes(ymin = ModelDev - 2*ModelDev_sd, ymax = ModelDev + 2*ModelDev_sd)) +
    ggtitle("Mean Absolute Deviation of Methods") +
    xlab("log-Sample Size") +
    ylab("MAD") +
    coord_cartesian(
      xlim =c(min(log(sample.size.vec)),max(log(sample.size.vec))),
      ylim = c(0,10)
    )

  #geom_errorbar(aes(ymin = lik.mean.scaled - 2*lik.sd.scaled, ymax = lik.mean.scaled + 2*lik.sd.scaled))

  #plt.rmse
  file.bias <- paste0("plot/Bias_SpilloverGATE_cluster_growth",cluster.growth,
                      "mutual_benefit_", mutual.benefit,
                      "cluster_equal_",cluster.equal.size,".png")
  file.rmse <- paste0("plot/RMSE_SpilloverGATE_cluster_growth",cluster.growth,
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





