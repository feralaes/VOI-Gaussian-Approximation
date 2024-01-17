# This is the OpenBUGs code.  I ended up not using the R code, because the shape parameter
# does not converge as can be shown by running the openbugs code below.
# For future reference, if I need to use the R code, make sure edit the inputs of the 
# JAGS function so that it sends a vector mu, and a matrix Tau.  The vector of x's stay the same.
# Also, retrieve both shape and scale thetas.  
# model{
#   theta0[1:2] ~ dmnorm(mu[], Tau[,])
#   theta[1] <- abs(theta0[1]) #shape
#   theta[2] <- abs(theta0[2]) #scale
#   lambda <- 1/pow(abs(theta[2]), theta[1])
#   for (i in 1:n)
#   {
#     x[i] ~ dweib(theta[1], lambda)
#   }
# }
# 
# 
# list(x = c(.139, .139, .139, .139, .139, .139, .139, .139, .139, .139), 
#      n = 10,
#      mu = c(1, 1), 
#      Tau = structure(.Data = c(1, 0, 0, 1), .Dim = c(2,2)))
# 
# list(theta0 = c(1, 1))
# list(theta0 = c(.5, .5))
# list(theta0 = c(5, 5))

library("R2jags")

max.GOF <- NULL
n0.pred <- NULL
n0list <- NULL

### Independent assumption: find ESS for each of the 4 parameters
# Scale of imatinib
test.model <- function() {
  #theta[1:2] ~ dmnorm(mu[], Tau[,])
  #v <- abs(theta[1]) #shape
  #theta ~ dnorm(1, 0.001)
  theta ~ dgamma(0.001, 0.001)
  v <- shape
  lambda <- 1/pow(theta, v)
  for (i in 1:n)
  {
    x[i] ~ dweib(v, lambda)
  }
}

### ============================================================================================
###   SCALE parameters

shape <- 1 #mean(X[,33]) #mean of shape parameter
nBurnIn = 10000
nIter = 100000 + nBurnIn

n0list <- seq(20,40,1) # list of prior samples to be taken
#n0list2 <- seq(130,150,1) # list of prior samples to be taken
#n0list3 <- seq(75,85,1) # list of prior samples to be taken
#n0listAll = list(n0list1, n0list2, n0list3)


#param <- 32 #scale imatinib
#for (param in 30:32){
#paramIter <- param - 29
#n0list <- n0listAll[[paramIter]]


#prior.theta0 <- X[,param]
prior.theta.mean <- 1 #mean(prior.theta0)
prior.theta.sd <- 0.2 #sd(prior.theta0)


prior.theta <- rnorm(nIter - nBurnIn, prior.theta.mean, prior.theta.sd)

lineinits <- NULL

#z <- 1/rgamma(1000000, shape = 1, scale = 100)
#hist(z, 10000)
#prior.theta <- rgamma(nIter - nBurnIn, shape = alpha, scale = beta) # prior distribution
hist(prior.theta, 100)
prior.theta.sorted <- sort(prior.theta) #sort the prior
#prior.theta.mean <- mean(prior.theta) # prior mean of the parameter needed for R2


len.n0list <- length(n0list)
gof <- NULL


for (i in 1:len.n0list){ # iterate through the sample sizes and find the one that maximizes R2.
  n0 <- n0list[i]
  print(n0)
  #n0 <- 20
  x <- rep(prior.theta.mean * factorial(1/shape), n0)
  #x <- prior.theta.mean
  linedata <- list(x = x, n = n0, shape = shape) 
  test.sim <- try(jags(data = linedata, inits = lineinits, model.file = test.model,
                       parameters.to.save = c("theta"),
                       n.chains = 1, n.iter = nIter, n.burnin = nBurnIn, n.thin = 1,
                       progress.bar = "text"))
  if(inherits(test.sim, "try-error")) # if error occurs, then iterate to the next sample size.
  {
    #error handling code, maybe just skip this iteration using
    gof[i] <- 0
    next
  }
  #rest of iteration for case of no error
  post.theta <- test.sim$BUGSoutput$sims.matrix[,2] #get the posterior distribution sample
  post.theta.sorted <- sort(post.theta) #sort the posterior
  SSerr <- sum((post.theta.sorted - prior.theta.sorted)^2) # sum of squared errors
  SStot <- sum((post.theta.sorted - prior.theta.mean)^2) # sum of total squared
  gof[i] <- 1 - SSerr/SStot # R-squared. 
}
hist(post.theta, 100)
plot(n0list, gof) # plot R2 against the sample sizes
max.GOF <- max(gof[!is.na(gof)]) # maximum R2
n0.pred <- n0list[which.max(gof)] # n0 that mazimizes R2.
#}

plot(post.theta.sorted)
cdfpost <- ecdf(post.theta.sorted)
plot(cdfpost)

# =====================================================================
# Take additional samples to update the prior. conduct data collection
# =====================================================================
test.model <- function() {
  #theta[1:2] ~ dmnorm(mu[], Tau[,])
  #v <- abs(theta[1]) #shape
  tau <- 1/pow(0.2,2)
  theta ~ dnorm(1, tau)
  #theta ~ dgamma(0.001, 0.001)
  v <- shape
  lambda <- 1/pow(theta, v)
  for (i in 1:n)
  {
    x[i] ~ dweib(v, lambda)
  }
}

nBurnIn = 1000
nIter = 10000 + nBurnIn


nlist <- c(10,100,1000)
nSample <- 1000
post.theta.mean <- matrix(0,nSample,3)
post.theta.var <- matrix(0,nSample,3)
for (i in 1:3){ # iterate through n.
  for (j in 1:nSample){
    n <- nlist[i]
    print(n)
    print(j)
    x <- NULL
    #n0 <- 20
    x0 <- rnorm(1, 1, 0.2) #rep(prior.theta.mean * factorial(1/shape), n)
    for (k in 1:n){
      x[k] <- rweibull(1,1,x0)
    }
    #x <- prior.theta.mean
    lineinits <- function(){list("theta" = x0)}
    linedata <- list(x = x, n = n, shape = shape) 
    test.sim <- try(jags(data = linedata, inits = lineinits, model.file = test.model,
                         parameters.to.save = c("theta", ),
                         n.chains = 1, n.iter = nIter, n.burnin = nBurnIn, n.thin = 1,
                         progress.bar = "none")) #"text"
    if(inherits(test.sim, "try-error")) # if error occurs, then iterate to the next sample size.
    {
      #error handling code, maybe just skip this iteration using
      gof[i] <- 0
      next
    }
    #rest of iteration for case of no error
    post.theta <- test.sim$BUGSoutput$sims.array[,,2] #matrix[,2] #get the posterior distribution sample
    #post.theta.sorted[,i] <- sort(post.theta) #sort the posterior
    post.theta.mean[j,i] <- mean(post.theta)
    post.theta.var[j,i] = var(post.theta)
    
  }
}
write.csv(post.theta.mean, "BUNormWeib.csv")
### Beta-binomial proof:

post.theta.var# = var(post.theta)
prior.theta.var <- prior.theta.sd ^ 2
N <- matrix(1,nSample,1) %*% nlist 
predn0 <- N * post.theta.var / (prior.theta.var - post.theta.var)
colMeans(predn0)
hist(predn0[,1], xlim = c(0,100), breaks = 300)
mean(predn0[predn0[,3]>0,3])

# =====================================================================
#       Testing n0 equation
# =====================================================================
test.model <- function() {
  #theta[1:2] ~ dmnorm(mu[], Tau[,])
  #v <- abs(theta[1]) #shape
  tau <- 1/pow(0.2,2)
  theta ~ dnorm(1, tau)
  #theta ~ dgamma(0.001, 0.001)
  v <- shape
  lambda <- 1/pow(theta, v)
  for (i in 1:n)
  {
    x[i] ~ dweib(v, lambda)
  }
}

nBurnIn = 1000
nIter = 10000 + nBurnIn


nlist <- c(10,100,1000)
nSample <- 1000
post.theta.mean <- matrix(0,nSample,3)
post.theta.var <- matrix(0,nSample,3)
for (i in 1:3){ # iterate through n.
  for (j in 1:nSample){
    n <- nlist[i]
    print(n)
    print(j)
    x <- NULL
    #n0 <- 20
    x0 <- rnorm(1, 1, 0.2) #rep(prior.theta.mean * factorial(1/shape), n)
    for (k in 1:n){
      x[k] <- rweibull(1,1,x0)
    }
    #x[1:n] <- mean(rweibull(10000,1,1))
    
    #x <- prior.theta.mean
    linedata <- list(x = x, n = n, shape = shape) 
    test.sim <- try(jags(data = linedata, inits = lineinits, model.file = test.model,
                         parameters.to.save = c("theta"),
                         n.chains = 1, n.iter = nIter, n.burnin = nBurnIn, n.thin = 1,
                         progress.bar = "none")) #"text"
    if(inherits(test.sim, "try-error")) # if error occurs, then iterate to the next sample size.
    {
      #error handling code, maybe just skip this iteration using
      gof[i] <- 0
      next
    }
    #rest of iteration for case of no error
    post.theta <- test.sim$BUGSoutput$sims.array[,,2] #matrix[,2] #get the posterior distribution sample
    #post.theta.sorted[,i] <- sort(post.theta) #sort the posterior
    post.theta.mean[j,i] <- mean(post.theta)
    post.theta.var[j,i] = var(post.theta)
  }
}
write.csv(post.theta.mean, "BUNormWeib.csv")
### Beta-binomial proof:

post.theta.var <- var(post.theta)
prior.theta.var <- prior.theta.sd ^ 2
N <- matrix(1,nSample,1) %*% nlist 
predn0 <- N * post.theta.var / (prior.theta.var - post.theta.var)
colMeans(predn0)
hist(predn0[,1], xlim = c(0,100), breaks = 300)
mean(predn0[predn0[,3]>0,3])

# =====================================================================
#       Testing n0 equation iterate through samples in BUGS instead of R.
# =====================================================================
test.model <- function() {
  #theta[1:2] ~ dmnorm(mu[], Tau[,])
  #v <- abs(theta[1]) #shape
  tau <- 1/pow(0.2,2)
  #tau0 ~ dgamma(0.001, 0.001)
  #mu ~ dnorm(1, tau0)
  v <- shape
  for (i in 1:m){
    theta[i] ~ dnorm(1, tau)
    lambda[i] <- 1/pow(theta[i], v)
    for (j in 1:n)
    {
      x[i,j] ~ dweib(v, lambda[i])
    }
  }
}

n <- 100


nBurnIn = 1000
nIter = 10000 + nBurnIn


nlist <- c(10,100,1000)
nSample <- 1000
post.theta.mean <- matrix(0,nSample,3)
post.theta.var <- matrix(0,nSample,3)
    n <- nlist[i]
    print(n)
    print(j)
    theta <- NULL
    #n0 <- 20
    theta <- rnorm(nSample, 1, 0.2) #rep(prior.theta.mean * factorial(1/shape), n)
    x <- matrix(0, nSample, n)
    for (i in 1:nSample){
      x[i, ] <- rweibull(n, 1, theta[i])
    }
    
    # var(x) 
    
    #x[1:n] <- mean(rweibull(10000,1,1))
    
    #x <- prior.theta.mean
    lineinits <- NULL
    
    linedata <- list(x = x, n = n, m = nSample, shape = shape) 
    test.sim <- try(jags(data = linedata, inits = lineinits, model.file = test.model,
                         parameters.to.save = c("theta"),
                         n.chains = 1, n.iter = nIter, n.burnin = nBurnIn, n.thin = 1,
                         progress.bar = "text")) 

    #rest of iteration for case of no error
    post.theta <- test.sim$BUGSoutput$sims.list$theta
    prepost <- colMeans(post.theta)
    prepost.var <- var(prepost)
    prior.var <- var(theta)
    
    n * (prior.var / prepost.var - 1)
    
    
    post.mu <- test.sim$BUGSoutput$sims.array[,,2] #matrix[,2] #get the posterior distribution sample
    post.mu.prepost <- var(post.mu)
    #post.theta.sorted[,i] <- sort(post.theta) #sort the posterior
    post.theta.mean[j,i] <- mean(post.theta)
    post.theta.var[j,i] = var(post.theta)
    
    n0 <- n * (var.mu.prior / var.mu.prepost - 1)
    
write.csv(post.theta.mean, "BUNormWeib.csv")
### Beta-binomial proof:

post.theta.var <- var(post.theta)
prior.theta.var <- prior.theta.sd ^ 2
N <- matrix(1,nSample,1) %*% nlist 
predn0 <- N * post.theta.var / (prior.theta.var - post.theta.var)
colMeans(predn0)
hist(predn0[,1], xlim = c(0,100), breaks = 300)
mean(predn0[predn0[,3]>0,3])



