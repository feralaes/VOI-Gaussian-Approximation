#### beta-binomial MCMC ####
library(R2jags)
library(matrixStats)
## 1. Obtain experimental data from the data likelihood in R
a <- 12 # alpha parameter for the beta prior distribution p~beta(a,b)
b <- 8 # beta parameter for the beta distribution
n <- 50 # additional data collections
nSample <- 2000
x <- matrix(0, nSample, 1) # initialize the new data matrix
p <- rbeta(nSample, a, b) # prior distribution of p
for (i in 1:nSample) { # For each prior value, conduct an experiment of size n
  x[i] <- rbinom(1, n, p[i]) # likelihood: x|p
}
## 2. Obtain the posterior distribution & compute the posterior means
test.model <- function() {
  for (i in 1:nSample) {
    p[i] ~ dbeta(a, b) #prior distribution of p
    x[i,1] ~ dbin(p[i], n) # likelihood: x|p
  } }
linedata <- list(x = x, n = n, nSample = nSample, a = a, b = b)
test.sim <- try(jags(data = linedata, inits = NULL, model.file = test.model,
                     parameters.to.save = c("p"),
                     n.chains = 1, n.iter = 11000, n.burnin = 1000, n.thin = 1,
                     progress.bar = "text"))
post.p <- test.sim$BUGSoutput$sims.list$p #get the posterior distribution sample
prepost <- colMeans(post.p) #obtain the preposterior distribution from the posterior mean of p
var.prior <- var(p) #prior variance
var.prepost <- var(prepost) #preposterior variance
n0.hat <- n * (var.prior / var.prepost - 1) #predicted prior sample size
n0.hat #21.17 compared to a+b = 12+8 = 20.