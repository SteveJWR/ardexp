
source("R/functions.R")
library(label.switching)

# Q = P(T = t|C = c)
Tr = 2 # number of traits
K = 2
# just one example of an easy to classify trait set
input.q = c(1,rep(0,K))
Q = matrix(rep(input.q, Tr)[1:(Tr*K)], nrow = Tr, ncol = K, byrow = T)
Q = t(t(Q)/colSums(Q))

#population versions
# pt = as.numeric(Q %*% pi)
# R = Q *outer(1/pt,pi,"*")
# R[is.nan(R)] = 0


pi = rep(1/K,K)

P = matrix(0.1, nrow = K, ncol = K )
P = (P + t(P))/2
diag(P) = rep(0.4,K)

n.sims = 200
m = 5
n.set = round(seq(100,500,length.out = m))


# Current block for this variance
prob = seq(0.3,0.7,length.out = K)

beta.Z = seq(0,1, length.out = K)
beta.a = 1

noise = 0.2
# Current Setup for K = 2
# only a fraction of treated neighbours.
beta.V.Z = c(1)

# beta1 is the intercept
true.param.vec <- c(beta.Z, beta.a, beta.V.Z)

p = length(true.param.vec)
results <- array(NA, c(n.sims,p,m))
results.true <- results
results.baseline <- results

set.seed(1)

for(i in seq(n.set)){
  n = n.set[i]
  print(paste0(i,"/",m))
  #beta.dtau = 50/sqrt(n)
  for(sim in seq(n.sims)){
    cat(paste0("Simulation ", sim,"/",n.sims), end = "\r")

    # TODO: can we make this faster?
    g.sbm = generateSBM(n,P,pi)


    G = g.sbm$G
    Z = g.sbm$groups


    traits = traitsSample(Z,Q)

    X = computeARD(traits, G)
    d.vec <- rowSums(X) #exhaustive ard

    # Assuming correct number of known clusters
    model.kmeans = clusterARDKmeans(X, K)
    Z.hat = model.kmeans$cluster

    #
    Z.hat = labelSwitching(Z,Z.hat)

    a = treatmentAssignment(Z.hat, prob)


    ## Simulate Outcomes

    diag.Z.blocks = array(0, c(n,n,K))
    d.vec = rowSums(X)
    for(k in seq(K)){
      idx = which(Z == k)
      block.mat = array(0,c(n,n))
      diag(block.mat)[idx] = 1
      diag.Z.blocks[,,k] = block.mat
    }

    # TODO: can we make this faster?
    V = as.numeric(G %*%  a)
    V = V/d.vec

    # new means
    means = beta.Z[Z] + beta.a*a +   beta.V.Z*V
    eps = rnorm(n, sd = noise)
    Y = means + eps
    # baselines
    Y.baseline = beta.Z[Z] + eps

    Q.hat = estimateQ(traits, Z.hat)
    Qp.hat = reverseConditionalQ(traits, Z.hat)

    #TODO: fix this in the code
    Qp.hat[is.nan(Qp.hat)] = 0

    V.ave = X %*% t(Qp.hat)
    V.ave = rowSums(t(t(V.ave) *prob)) /d.vec

    data = data.frame("Y" = Y, "A" = a, "Z" = as.factor(Z.hat), "V" = V.ave)

    fmla <- formula("Y ~ Z + A + V")
    model <- lm(fmla,data)

    data.true = data.frame("Y" = Y, "A" = a, "Z" = as.factor(Z), "V" = V)

    model.true <- lm(fmla,data.true)
    #data.sim[sim,] = model$coefficients
    results[sim,,i] = model$coefficients
    results.true[sim,,i] = model.true$coefficients

    data.baseline = data.frame("Y" = c(Y.baseline,Y), "A" = c(rep(0,n),a), "Z" = as.factor(c(Z.hat,Z.hat)), "V" = c(rep(0,n),V.ave))
    model.baseline <- lm(fmla,data.baseline)
    results.baseline[sim,,i] = model.baseline$coefficients

  }
}



#true.coefficients <- c(0,1,1,1)
#beta.dt.true = 50/sqrt(n.set)

#TODO: Fix this section
rmse.mat = matrix(NA, ncol = p, nrow = m)
bias.mat = matrix(NA, ncol = p, nrow = m)

rmse.true.mat = matrix(NA, ncol = p, nrow = m)
bias.true.mat = matrix(NA, ncol = p, nrow = m)

rmse.baseline.mat = matrix(NA, ncol = p, nrow = m)
bias.baseline.mat = matrix(NA, ncol = p, nrow = m)


for(i in seq(m)){
  for(j in seq(p)){
    rmse.mat[i,j] = sqrt(mean((results[,j,i] - true.param.vec[j])^2))
    bias.mat[i,j] = mean(results[,j,i] - true.param.vec[j])
  }
}


for(i in seq(m)){
  for(j in seq(p)){
    rmse.baseline.mat[i,j] = sqrt(mean((results.baseline[,j,i] - true.param.vec[j])^2))
    bias.baseline.mat[i,j] = mean(results.baseline[,j,i] - true.param.vec[j])
  }
}


for(i in seq(m)){
  for(j in seq(p)){
    rmse.true.mat[i,j] = sqrt(mean((results.true[,j,i] - true.param.vec[j])^2))
    bias.true.mat[i,j] = mean(results.true[,j,i] - true.param.vec[j])
  }
}

#rmse.dt <- sqrt(colMeans((t(t(data.sim.dt) - beta.dt.true))**2))

plot(log(n.set), log(rmse.mat[,1]), type = "l")
plot(log(n.set), log(rmse.mat[,2]), type = "l")
plot(log(n.set), log(rmse.mat[,3]), type = "l")
plot(log(n.set), log(rmse.mat[,4]), type = "l")
#plot(log(n.set), log(rmse.mat[,5]), type = "l")
#plot(log(n.set), log(rmse.mat[,6]), type = "l")
#plot(log(n.set), log(rmse.mat[,7]), type = "l")




plot(log(n.set), bias.mat[,1], type = "l")
plot(log(n.set), bias.mat[,2], type = "l")
plot(log(n.set), bias.mat[,3], type = "l")
plot(log(n.set), bias.mat[,4], type = "l")
#plot(log(n.set), bias.mat[,5], type = "l")
#plot(log(n.set), bias.mat[,6], type = "l")
#plot(log(n.set), bias.mat[,7], type = "l")


plot.dat.rmse= data.frame("rmse" = as.numeric(rmse.mat), "n" = rep(n.set, times = 4),"parameter" = as.factor(rep(seq(4),each = m)))
plt.rmse <- ggplot(data = plot.dat.rmse, aes(x = n, y = rmse, colour = parameter)) + ylim(0,1.2) + geom_line() + ggtitle("RMSE (treated neighbors)")
plt.rmse

plot.dat.bias= data.frame("bias" = as.numeric(bias.mat), "n" = rep(n.set, times = 4),"parameter" = as.factor(rep(seq(4),each = m)))
plt.bias <- ggplot(data = plot.dat.bias, aes(x = n, y = bias, colour = parameter)) + ylim(-1,1) + geom_line() + ggtitle("Bias (treated neighbors)")
plt.bias


#####
plot.dat.true.rmse = data.frame("rmse" = as.numeric(rmse.true.mat), "n" = rep(n.set, times = 4),"parameter" = as.factor(rep(seq(4),each = m)))
plt.true.rmse <- ggplot(data = plot.dat.true.rmse, aes(x = n, y = rmse, colour = parameter)) + ylim(0,1.2) + geom_line() + ggtitle("RMSE (treated neighbors), True Exposure")
plt.true.rmse

plot.dat.true.bias= data.frame("bias" = as.numeric(bias.true.mat), "n" = rep(n.set, times = 4),"parameter" = as.factor(rep(seq(4),each = m)))
plt.true.bias <- ggplot(data = plot.dat.true.bias, aes(x = n, y = bias, colour = parameter)) + ylim(-1,1) + geom_line() + ggtitle("Bias (treated neighbors), True Exposure")
plt.true.bias


##### .baseline
plot.dat.baseline.rmse= data.frame("rmse" = as.numeric(rmse.baseline.mat), "n" = rep(n.set, times = 4),"parameter" = as.factor(rep(seq(4),each = m)))
plt.baseline.rmse <- ggplot(data = plot.dat.baseline.rmse, aes(x = n, y = rmse, colour = parameter)) + ylim(0,0.2) + geom_line() + ggtitle("RMSE (treated neighbors)")
plt.baseline.rmse


plot.dat.baseline.bias= data.frame("bias" = as.numeric(bias.baseline.mat), "n" = rep(n.set, times = 4),"parameter" = as.factor(rep(seq(4),each = m)))
plt.baseline.bias <- ggplot(data = plot.dat.baseline.bias, aes(x = n, y = bias, colour = parameter)) + ylim(-1,1) + geom_line() + ggtitle("Bias (treated neighbors)")
plt.baseline.bias



plot(log(n.set), log(rmse.true.mat[,1]), type = "l")
plot(log(n.set), log(rmse.true.mat[,2]), type = "l")
plot(log(n.set), log(rmse.true.mat[,3]), type = "l")
plot(log(n.set), log(rmse.true.mat[,4]), type = "l")

#plot(log(n.set), log(rmse.true.mat[,5]), type = "l")
#plot(log(n.set), log(rmse.true.mat[,6]), type = "l")
#plot(log(n.set), log(rmse.true.mat[,7]), type = "l")



plot(log(n.set), bias.true.mat[,1], type = "l")
plot(log(n.set), bias.true.mat[,2], type = "l")
plot(log(n.set), bias.true.mat[,3], type = "l")
plot(log(n.set), bias.true.mat[,4], type = "l")
#plot(log(n.set), bias.true.mat[,5], type = "l")
#plot(log(n.set), bias.true.mat[,6], type = "l")
#plot(log(n.set), bias.true.mat[,7], type = "l")









# Allow for different slopes by connections to each degree
# Q = P(T = t|C = c)
Tr = 2 # number of traits
K = 2
# just one example of an easy to classify trait set
input.q = c(1,rep(0,K))
Q = matrix(rep(input.q, Tr)[1:(Tr*K)], nrow = Tr, ncol = K, byrow = T)
Q = t(t(Q)/colSums(Q))

pi = rep(1/K,K)

set.seed(1)
P = matrix(0.1, nrow = K, ncol = K )
P = (P + t(P))/2
diag(P) = rep(0.4,K)

n.sims = 200
m = 5
n.set = round(exp(seq(log(50),log(1000),length.out = m)))

# Current block for this variance
prob = seq(0.3,0.7,length.out = K)

beta.Z = seq(0,1, length.out = K)
beta.a = 1

noise = 0.1
# Current Setup for K = 2
beta.V.Z = c(1,5)

# beta1 is the intercept
true.param.vec <- c(beta.Z, beta.a, beta.V.Z)

p = length(true.param.vec)
results <- array(NA, c(n.sims,p,m))
results.true <- results



set.seed(1)

for(i in seq(n.set)){
  n = n.set[i]
  print(paste0(i,"/",m))
  #beta.dtau = 50/sqrt(n)
  for(sim in seq(n.sims)){
    cat(paste0("Simulation ", sim,"/",n.sims), end = "\r")

    # TODO: can we make this faster?
    g.sbm = generateSBM(n,P,pi)


    G = g.sbm$G
    Z = g.sbm$groups


    traits = traitsSample(Z,Q)

    X = computeARD(traits, G)
    d.vec <- rowSums(X) #exhaustive ard

    # Assuming correct number of known clusters
    model.kmeans = clusterARDKmeans(X, K)
    Z.hat = model.kmeans$cluster

    #
    Z.hat = labelSwitching(Z,Z.hat)

    a = treatmentAssignment(Z.hat, prob)


    ## Simulate Outcomes

    diag.Z.blocks = array(0, c(n,n,K))
    d.vec = rowSums(X)
    for(k in seq(K)){
      idx = which(Z == k)
      block.mat = array(0,c(n,n))
      diag(block.mat)[idx] = 1
      diag.Z.blocks[,,k] = block.mat
    }

    # TODO: can we make this faster?
    V = array(0,c(n,K))
    for(k in seq(K)){
      V[,k] = as.numeric(G %*% diag.Z.blocks[,,k] %*% a)
    }
    V = V/d.vec

    # new means
    means = beta.Z[Z] + beta.a*a +   V %*% beta.V.Z
    eps = rnorm(n, sd = noise)
    Y = means + eps

    Q.hat = estimateQ(traits, Z.hat)
    Qp.hat = reverseConditionalQ(traits, Z.hat)

    #TODO: fix this in the code
    Qp.hat[is.nan(Qp.hat)] = 0

    V.ave = X %*% t(Qp.hat)
    V.ave = t(t(V.ave) *prob) /d.vec

    data = data.frame("Y" = Y, "A" = a, "Z" = as.factor(Z.hat), "V" = V.ave)

    fmla <- formula("Y ~ Z + A + V")
    model <- lm(fmla,data)

    data.true = data.frame("Y" = Y, "A" = a, "Z" = as.factor(Z), "V" = V)

    model.true <- lm(fmla,data.true)
    #data.sim[sim,] = model$coefficients
    results[sim,,i] = model$coefficients
    results.true[sim,,i] = model.true$coefficients
  }
}



#true.coefficients <- c(0,1,1,1)
#beta.dt.true = 50/sqrt(n.set)

#TODO: Fix this section
rmse.mat = matrix(NA, ncol = p, nrow = m)
bias.mat = matrix(NA, ncol = p, nrow = m)

rmse.true.mat = matrix(NA, ncol = p, nrow = m)
bias.true.mat = matrix(NA, ncol = p, nrow = m)


for(i in seq(m)){
  for(j in seq(p)){
    rmse.mat[i,j] = sqrt(mean((results[,j,i] - true.param.vec[j])^2))
    bias.mat[i,j] = mean(results[,j,i] - true.param.vec[j])
  }
}


for(i in seq(m)){
  for(j in seq(p)){
    rmse.true.mat[i,j] = sqrt(mean((results.true[,j,i] - true.param.vec[j])^2))
    bias.true.mat[i,j] = mean(results.true[,j,i] - true.param.vec[j])
  }
}

#rmse.dt <- sqrt(colMeans((t(t(data.sim.dt) - beta.dt.true))**2))

plot(log(n.set), log(rmse.mat[,1]), type = "l")
plot(log(n.set), log(rmse.mat[,2]), type = "l")
plot(log(n.set), log(rmse.mat[,3]), type = "l")
plot(log(n.set), log(rmse.mat[,4]), type = "l")
plot(log(n.set), log(rmse.mat[,5]), type = "l")
#plot(log(n.set), log(rmse.mat[,6]), type = "l")
#plot(log(n.set), log(rmse.mat[,7]), type = "l")



plot(log(n.set), bias.mat[,1], type = "l")
plot(log(n.set), bias.mat[,2], type = "l")
plot(log(n.set), bias.mat[,3], type = "l")
plot(log(n.set), bias.mat[,4], type = "l")
plot(log(n.set), bias.mat[,5], type = "l")
#plot(log(n.set), bias.mat[,6], type = "l")
#plot(log(n.set), bias.mat[,7], type = "l")


plot(log(n.set), log(rmse.true.mat[,1]), type = "l")
plot(log(n.set), log(rmse.true.mat[,2]), type = "l")
plot(log(n.set), log(rmse.true.mat[,3]), type = "l")
plot(log(n.set), log(rmse.true.mat[,4]), type = "l")
plot(log(n.set), log(rmse.true.mat[,5]), type = "l")
#plot(log(n.set), log(rmse.true.mat[,6]), type = "l")
#plot(log(n.set), log(rmse.true.mat[,7]), type = "l")



plot(log(n.set), bias.true.mat[,1], type = "l")
plot(log(n.set), bias.true.mat[,2], type = "l")
plot(log(n.set), bias.true.mat[,3], type = "l")
plot(log(n.set), bias.true.mat[,4], type = "l")
plot(log(n.set), bias.true.mat[,5], type = "l")
#plot(log(n.set), bias.true.mat[,6], type = "l")
#plot(log(n.set), bias.true.mat[,7], type = "l")








# Allow for different slopes by connections to each degree

# Q = P(T = t|C = c)
Tr = 2 # number of traits
K = 2
# just one example of an easy to classify trait set
input.q = c(1,rep(0,K))
Q = matrix(rep(input.q, Tr)[1:(Tr*K)], nrow = Tr, ncol = K, byrow = T)
Q = t(t(Q)/colSums(Q))
K = 2

pi = rep(1/K,K)

set.seed(1)
P = matrix(0.1, nrow = K, ncol = K )
P = (P + t(P))/2
diag(P) = rep(0.4,K)

n.sims = 200
m = 5
n.set = round(exp(seq(log(50),log(1000),length.out = m)))

# Current block for this variance
prob = seq(0.3,0.7,length.out = K)

data.sim.int <- matrix(nrow = n.sims, ncol = length(n.set))
data.sim.Z <- matrix(nrow = n.sims, ncol = length(n.set))
data.sim.a <- matrix(nrow = n.sims, ncol = length(n.set))
data.sim.dt <- matrix(nrow = n.sims, ncol = length(n.set))

beta.Z = seq(0,1, length.out = K)
beta.a = 1

noise = 0.1
# Current Setup for K = 2
beta.V.Z = matrix(c(1,5,5,1), nrow = K, ncol = K)

# beta1 is the intercept
true.param.vec <- c(beta.Z, beta.a, beta.V.Z)

p = length(true.param.vec)
results <- array(NA, c(n.sims,p,m))
results.true <- results


set.seed(1)

for(i in seq(n.set)){
  n = n.set[i]
  print(paste0(i,"/",m))
  #beta.dtau = 50/sqrt(n)
  for(sim in seq(n.sims)){
    cat(paste0("Simulation ", sim,"/",n.sims), end = "\r")

    # TODO: can we make this faster?
    g.sbm = generateSBM(n,P,pi)


    G = g.sbm$G
    Z = g.sbm$groups


    traits = traitsSample(Z,Q)

    X = computeARD(traits, G)
    d.vec <- rowSums(X) #exhaustive ard

    # Assuming correct number of known clusters
    model.kmeans = clusterARDKmeans(X, K)
    Z.hat = model.kmeans$cluster

    #
    Z.hat = labelSwitching(Z,Z.hat)

    a = treatmentAssignment(Z.hat, prob)


    ## Simulate Outcomes

    diag.Z.blocks = array(0, c(n,n,K))
    d.vec = rowSums(X)
    for(k in seq(K)){
      idx = which(Z == k)
      block.mat = array(0,c(n,n))
      diag(block.mat)[idx] = 1
      diag.Z.blocks[,,k] = block.mat
    }

    # TODO: can we make this faster?
    V = array(0,c(n,K))
    for(k in seq(K)){
      V[,k] = as.numeric(G %*% diag.Z.blocks[,,k] %*% a)
    }
    V = V/d.vec

    # new means
    means = beta.Z[Z] + beta.a*a +  rowSums(beta.V.Z[Z,] * V)
    eps = rnorm(n, sd = noise)
    Y = means + eps

    Q.hat = estimateQ(traits, Z.hat)
    Qp.hat = reverseConditionalQ(traits, Z.hat)

    #TODO: fix this in the code
    Qp.hat[is.nan(Qp.hat)] = 0

    V.ave = X %*% t(Qp.hat)
    V.ave = t(t(V.ave) *prob) /d.vec

    data = data.frame("Y" = Y, "A" = a, "Z" = as.factor(Z.hat), "V" = V.ave)

    fmla <- formula("Y ~ Z + A + V:Z")
    model <- lm(fmla,data)

    data.true = data.frame("Y" = Y, "A" = a, "Z" = as.factor(Z), "V" = V)

    model.true <- lm(fmla,data.true)
    #data.sim[sim,] = model$coefficients
    results[sim,,i] = model$coefficients
    results.true[sim,,i] = model.true$coefficients
  }
}



#true.coefficients <- c(0,1,1,1)
#beta.dt.true = 50/sqrt(n.set)

#TODO: Fix this section
rmse.mat = matrix(NA, ncol = p, nrow = m)
bias.mat = matrix(NA, ncol = p, nrow = m)

rmse.true.mat = matrix(NA, ncol = p, nrow = m)
bias.true.mat = matrix(NA, ncol = p, nrow = m)


for(i in seq(m)){
  for(j in seq(p)){
    rmse.mat[i,j] = sqrt(mean((results[,j,i] - true.param.vec[j])^2))
    bias.mat[i,j] = mean(results[,j,i] - true.param.vec[j])
  }
}


for(i in seq(m)){
  for(j in seq(p)){
    rmse.true.mat[i,j] = sqrt(mean((results.true[,j,i] - true.param.vec[j])^2))
    bias.true.mat[i,j] = mean(results.true[,j,i] - true.param.vec[j])
  }
}

#rmse.dt <- sqrt(colMeans((t(t(data.sim.dt) - beta.dt.true))**2))

plot(log(n.set), log(rmse.mat[,1]), type = "l")
plot(log(n.set), log(rmse.mat[,2]), type = "l")
plot(log(n.set), log(rmse.mat[,3]), type = "l")
plot(log(n.set), log(rmse.mat[,4]), type = "l")
plot(log(n.set), log(rmse.mat[,5]), type = "l")
plot(log(n.set), log(rmse.mat[,6]), type = "l")
plot(log(n.set), log(rmse.mat[,7]), type = "l")



plot(log(n.set), bias.mat[,1], type = "l")
plot(log(n.set), bias.mat[,2], type = "l")
plot(log(n.set), bias.mat[,3], type = "l")
plot(log(n.set), bias.mat[,4], type = "l")
plot(log(n.set), bias.mat[,5], type = "l")
plot(log(n.set), bias.mat[,6], type = "l")
plot(log(n.set), bias.mat[,7], type = "l")


plot(log(n.set), log(rmse.true.mat[,1]), type = "l")
plot(log(n.set), log(rmse.true.mat[,2]), type = "l")
plot(log(n.set), log(rmse.true.mat[,3]), type = "l")
plot(log(n.set), log(rmse.true.mat[,4]), type = "l")
plot(log(n.set), log(rmse.true.mat[,5]), type = "l")
plot(log(n.set), log(rmse.true.mat[,6]), type = "l")
plot(log(n.set), log(rmse.true.mat[,7]), type = "l")



plot(log(n.set), bias.true.mat[,1], type = "l")
plot(log(n.set), bias.true.mat[,2], type = "l")
plot(log(n.set), bias.true.mat[,3], type = "l")
plot(log(n.set), bias.true.mat[,4], type = "l")
plot(log(n.set), bias.true.mat[,5], type = "l")
plot(log(n.set), bias.true.mat[,6], type = "l")
plot(log(n.set), bias.true.mat[,7], type = "l")






