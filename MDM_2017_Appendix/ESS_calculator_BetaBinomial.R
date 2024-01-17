library("R2jags")

test.model <- function() {
  p ~ dbeta(1,1)
  x ~ dbin(p, n)
}

curve(dbeta(x, 1,1))

x <- 30
lineinits <- NULL
nBurnIn = 1000
nIter = 10000 + nBurnIn

p0 <- rbeta(nIter - nBurnIn, 30, 70) # prior distribution
p0sorted <- sort(p0) #sort the prior
meanp0 <- mean(p0) # prior mean of the parameter needed for R2
n0list <- seq(1,200) # list of prior samples to be taken
gof <- NULL
for (n0 in n0list){ # iterate through the sample sizes and find the one that maximizes R2.
  print(n0)
linedata <- list(x = x, n = n0) 
test.sim <- try(jags(data = linedata, inits = lineinits, model.file = test.model,
                 parameters.to.save = c("p"),
                 n.chains = 1, n.iter = nIter, n.burnin = nBurnIn, n.thin = 1,
                 progress.bar = "none"))
if(inherits(test.sim, "try-error")) # if error occurs, then iterate to the next sample size.
{
  #error handling code, maybe just skip this iteration using
  next
}
#rest of iteration for case of no error
postp <- test.sim$BUGSoutput$sims.matrix[,2] #get the posterior distribution sample
var.postp <- var(postp)

p1sorted <- sort(postp) #sort the posterior
SSerr = sum((p1sorted - p0sorted)^2) # sum of squared errors
SStot = sum((p1sorted - meanp0)^2) # sum of total squared
gof[n0] = 1 - SSerr/SStot # R-squared. 
}

plot(n0list, gof) # plot R2 against the sample sizes
maxGOF <- max(gof[!is.na(gof)]) # maximum R2
n0pred <- n0list[which.max(gof)] # n0 that mazimizes R2.

## =====================================================================
# TEST N0 Posterior EQUATION
## ===================================================================

library("R2jags")

test.model <- function() {
  p ~ dbeta(a, b)
  x ~ dbin(p, n)
}

curve(dbeta(x, 1,1))

x <- 30
lineinits <- NULL
nBurnIn = 1000
nIter = 10000 + nBurnIn
nSample <- nIter - nBurnIn
p0 <- rbeta(nIter - nBurnIn, 30, 70) # prior distribution
p0sorted <- sort(p0) #sort the prior
meanp0 <- mean(p0) # prior mean of the parameter needed for R2


print(n0)
n <- 60
  
#x <- NULL
a <- 3
b <- 7
p0 <- rbeta(nSample, a, b) # prior distribution
for (i in 1:nSample){
  x[i] <- rbinom(1, n, p0[i])
}
#for (i in 1:nSample){
  linedata <- list(x = round(mean(x)), n = n, a = a, b = b) 
  test.sim <- try(jags(data = linedata, inits = lineinits, model.file = test.model,
                       parameters.to.save = c("p"),
                       n.chains = 1, n.iter = nIter, n.burnin = nBurnIn, n.thin = 1,
                       progress.bar = "none"))
  postp <- test.sim$BUGSoutput$sims.matrix[,2] #get the posterior distribution sample
  
#}
#rest of iteration for case of no error

var.post <- var(postp)
var.prior <- var(p0)
n0pred <- n * (var.prior / var.prepost - 1)
print(n0pred)

n * var.post / (var.prior - var.post)

## =====================================================================
# TEST N0 EQUATION: Preposterior equation.
## ===================================================================

library("R2jags")

test.model <- function() {
  p ~ dbeta(a, b)
  x ~ dbin(p, n)
}

curve(dbeta(x, 1,1))

#x <- 30
lineinits <- NULL
nBurnIn = 1000
nIter = 10000 + nBurnIn
nSample <- nIter - nBurnIn
p0 <- rbeta(nIter - nBurnIn, 30, 70) # prior distribution
p0sorted <- sort(p0) #sort the prior
meanp0 <- mean(p0) # prior mean of the parameter needed for R2


n <- 60

#x <- NULL
a <- 30
b <- 70
p0 <- rbeta(nSample, a, b) # prior distribution
for (i in 1:nSample){
  x[i] <- rbinom(1, n, p0[i])
}
#for (i in 1:nSample){
postp <- matrix(0, nSample, nSample)
for (i in 1:nSample){
  print(i)
  linedata <- list(x = x[i], n = n, a = a, b = b) 
  test.sim <- try(jags(data = linedata, inits = lineinits, model.file = test.model,
                       parameters.to.save = c("p"),
                       n.chains = 1, n.iter = nIter, n.burnin = nBurnIn, n.thin = 1,
                       progress.bar = "none"))
  postp[i,] <- test.sim$BUGSoutput$sims.matrix[,2] #get the posterior distribution sample
}
#}
#rest of iteration for case of no error
prepost <- rowMeans(postp)

var.post <- rowSds(postp)
var.prior <- var(p0)
var.prepost <- var(prepost)
n * (var.prior / var.prepost - 1)

n * var.post / (var.prior - var.post)

## =====================================================================
# TEST N0 EQUATION: Preposterior equation. ALL IN JAGS
## ===================================================================

library("R2jags")
library(matrixStats)
test.model <- function() { 
  for (i in 1:nSample){
    p[i] ~ dbeta(a, b) #prior distribution of p
    x[i,1] ~ dbin(p[i], n) # likelihood: x|p
  }
}

a <- 12 # alpha parameter for the beta distribution p~beta(a,b)
b <- 8
n <- 50 # additional data collections

nSample <- 10000
x <- matrix(0, nSample, 1) # initialize the new data
p <- rbeta(nSample, a, b) # prior distribution of p
for (i in 1:nSample){ # For each prior value, conduct an experiment of size n
  x[i] <- rbinom(1, n, p[i]) # likelihood: x|p
}
linedata <- list(x = x, n = n, nSample = nSample, a = a, b = b) 
test.sim <- try(jags(data = linedata, inits = NULL, model.file = test.model,
                     parameters.to.save = c("p"),
                     n.chains = 1, n.iter = 11000, n.burnin = 1000, n.thin = 1,
                     progress.bar = "text"))
post.p <- test.sim$BUGSoutput$sims.list$p #get the posterior distribution sample
  
prepost <- colMeans(post.p) #calculate the preposterior = posterior mean of p

var.prior <- var(p) #prior variance
var.prepost <- var(prepost) #preposterior variance
n0.hat <- n * (var.prior / var.prepost - 1) #predicted prior sample size
n0.hat #21.17 compared to a+b = 12+8 = 20.
