

## Simuations to compare the performance against both full data and Ugander Methods

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
b = -1
delta = 1

# simulation noise
sigma = 0.5 # how much noise are we willing to permit here?



if(mutual.benefit){
  gamma = 1
} else {
  gamma = -1
}


# Think of a second version of the simulation where there is a competative advantage to treatment
# but that wears off if everyone else around you has it

#Global average treatment effect
gate = gamma + delta

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
n.methods = 5
results <- array(NA, c(n.sims, n.methods, J))

# Methods tuning parameters
B.boot = 50
#### Simulations Start

# formula for the regression model
fmla <- formula(Y ~ 0 + deg.ratio +  deg.ratio:(H + A + H:A + frac.treated + H:frac.treated))
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
      PI = seq(1,sqrt(K),length.out = K)
      PI = PI/sum(PI)
    }

    g.true <- generateSBM(n,P,PI)
    G.true <- g.true$G
    Z.true = g.true$groups

    H <- simH(Z.true)


    # whether to do binomial randomization or the true cluster treatment
    A.clust <- clusterTreatment(Z.true, p.treat)



    K.high.level <- round(sat.frac*K)
    K.low.level <- K - K.high.level
    clust.levels <- c(rep(p.high.level, K.high.level), rep(p.low.level, K.low.level))
    A.sat <- saturationRandomizationTreatment(Z.true, levels = clust.levels)


    outcome.clust <- simUganderModelOutcomes(G.true,H,A.clust,a = a, b = b, delta = delta, gamma = gamma, sigma = sigma, scale.deg = T)
    outcome.sat <- simUganderModelOutcomes(G.true,H,A.sat,a = a, b = b, delta = delta, gamma = gamma, sigma = sigma, scale.deg = T)
    Y.clust <- outcome.clust$Y
    Y.sat <- outcome.sat$Y

    gate = outcome.clust$GATE

    d.vec = as.numeric(G.true %*% rep(1,n))
    d.mean = mean(d.vec)
    f <- fracTreated(G.true,A.sat)
    data <- UganderCovariates(G.true,H,A.sat)
    data <- data.frame("Y" = Y.sat, "A" = A.sat, "frac.treated" = f, "H" = H, "deg.ratio" = d.vec/d.mean)


    data.0 <- data
    data.1 <- data
    data.0$A = 0
    data.1$A = 1
    data.0$frac.treated = 0
    data.1$frac.treated = 1

    model <- lm(fmla, data = data)

    Y.0 <- predict(model, data.0)
    Y.1 <- predict(model, data.1)

    # full data version.
    # (Note I had e)
    gate.est1 <- model$coefficients[3]/model$coefficients[1] + model$coefficients[4]/model$coefficients[1]
    gate.est2 <- model$coefficients[5]/model$coefficients[2] + model$coefficients[6]/model$coefficients[2]
    gate.est3 <- (gate.est1 + gate.est2)/2

    gate.est <- mean(Y.1) - mean(Y.0)

    # Here we are asking about the traits directly for simplicity in the simulations.

    traits <- g.true$groups

    X <-computeARD(traits, G.true)

    # Assuming For simplicity, K is known.
    k.means.model <- clusterARDKmeans(X,K)


    Z.hat = k.means.model$cluster
    #
    Z.hat = labelSwitching(Z.true,Z.hat)


    P.hat <- estimatePmat(Z.hat, G.true)


    res.ard <- ARDSBMLinearRegressionSim(Y.sat, fmla, UganderCovariates, A.sat, H, P.hat,Z.hat, B.boot = B.boot, verbose = T)

    ard.avg.coef <- meanOverList(res.ard$coef)
    ard.avg.data <- meanOverList(res.ard$data)

    ard.gate.est1 <- ard.avg.coef[3]/ard.avg.coef[1] + ard.avg.coef[4]/ard.avg.coef[1]
    ard.gate.est2 <- ard.avg.coef[5]/ard.avg.coef[2] + ard.avg.coef[6]/ard.avg.coef[2]
    ard.gate.est3 <- (ard.gate.est1 + ard.gate.est2)/2

    ard.gate.est <- meanOverList(res.ard$gate)

    res.ard.true.model <- ARDSBMLinearRegressionSim(Y.sat, fmla, UganderCovariates, A.sat, H, P,Z.true, B.boot = B.boot, verbose = T)


    ard.true.model.avg.coef <- meanOverList(res.ard.true.model$coef)

    ard.true.model.data <- meanOverList(res.ard.true.model$data)

    ard.true.model.gate.est1 <- ard.true.model.avg.coef[3]/ard.true.model.avg.coef[1] + ard.true.model.avg.coef[4]/ard.true.model.avg.coef[1]
    ard.true.model.gate.est2 <- ard.true.model.avg.coef[5]/ard.true.model.avg.coef[2] + ard.true.model.avg.coef[6]/ard.true.model.avg.coef[2]
    ard.true.model.gate.est3 <- (ard.true.model.gate.est1 + ard.true.model.gate.est2)/2

    ard.true.model.gate.est <-  meanOverList(res.ard.true.model$gate)

    #data.mean.sub <- data
    #data.mean.sub$frac.treated <- 0.2630463*(data.mean.sub$frac.treated < 0.5) + 0.7369646*(data.mean.sub$frac.treated > 0.5)
    #cor(data.mean.sub)
    #model.mean.sub <- lm(fmla, data.mean.sub)

    #model.mean.sub$coefficients
    #TODO: make the HT estimator even better
    gate.HT <- HTEstimatorCluster(Y.clust,G.true,A.clust,c.vec = Z.true,p.treat = p.treat)


    gate.DM <- diffMeansEstimator(Y.clust,A.clust)

    res.vec <- c(ard.gate.est1 - gate, ard.gate.est2 - gate, ard.gate.est3 - gate,
                 ard.true.model.gate.est1 - gate, ard.true.model.gate.est2 - gate,
                 ard.true.model.gate.est3 - gate,
                 gate.est1 - gate, gate.est2 - gate, gate.est3 - gate,
                 gate.HT - gate, gate.DM - gate)
    res.vec <- c(ard.gate.est - gate,
                 ard.true.model.gate.est - gate,
                 gate.est - gate,
                 gate.HT - gate, gate.DM - gate)
    names(res.vec) = NULL
    results[sim,,j] <- res.vec
  }
}


colnames(results) = c("ard",
                      "ard.tm",
                      "reg",
                      "HT", "DM")

filename <- paste0("data/UganderGATE_cluster_growth",cluster.growth,
                   "mutual_benefit_", mutual.benefit,
                   "cluster_equal_",cluster.equal.size,
                   "block",block,".rds")

saveRDS(results, filename)




make.plots = F
if(make.plots){
  for(id in seq(8)){
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
    library(ggpubr)
    library(abind)
    png.width = 1200
    png.height = 1000
    png.res = 200


    filename <- paste0("data/UganderGATE_cluster_growth",cluster.growth,
                       "mutual_benefit_", mutual.benefit,
                       "cluster_equal_",cluster.equal.size,
                       "block",1,".rds")
    results = readRDS(filename)
    for(i in seq(2,100)){
      filename <- paste0("data/UganderGATE_cluster_growth",cluster.growth,
                         "mutual_benefit_", mutual.benefit,
                         "cluster_equal_",cluster.equal.size,
                         "block",i,".rds")
      res.tmp <-readRDS(filename)
      results <- abind(results, res.tmp, along = 1)
    }
    n.sims = dim(results)[1]

    methods <- c("ard",
                 "ard.tm",
                 "reg",
                 "HT", "DM")

    J = length(methods)
    n.seq <- c(100,316,1000,3162,10000)
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
    method.subset = c("ard",
                      "ard.tm",
                      "reg",
                      "HT", "DM")
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
      ggtitle("RMSE of Methods") +
      xlab("log-Sample Size") +
      ylab("RMSE") +
      coord_cartesian(
        xlim =c(min(log(sample.size.vec)),max(log(sample.size.vec))),
        ylim = c(0,10)
      )
    plt.rmse
    #geom_errorbar(aes(ymin = lik.mean.scaled - 2*lik.sd.scaled, ymax = lik.mean.scaled + 2*lik.sd.scaled))

    #plt.rmse
    file.bias <- paste0("plot/Bias_UganderGATE_cluster_growth",cluster.growth,
                        "mutual_benefit_", mutual.benefit,
                        "cluster_equal_",cluster.equal.size,".png")
    file.rmse <- paste0("plot/RMSE_UganderGATE_cluster_growth",cluster.growth,
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





