library("R2jags")
library(matrixStats)

test.model <- function() {
  for (i in 1:nSample){
    pC[i] ~ dbeta(a, b)
    LOR[i] ~ dnorm(mu, tau)
    OR[i] <- exp(LOR[i])
    pT[i] <- OR[i] * pC[i] / (1 + pC[i] * (OR[i] - 1))
    xC[i] ~ dbin(pC[i], n)
    xT[i] ~ dbin(pT[i], n)
  }
}

nSample <- 5000
a <- 15
b <- 85
mu <- -1.5
tau <- 3
sigma <- sqrt(1/tau)
pC <- rbeta(nSample, a, b) # prob critical event. 
LOR <- rnorm(nSample, mu, sigma) 
OR <- exp(LOR)
pT <- OR * pC / (1 + pC * (OR - 1))  
mean(pT) 
# 0.433418 

n <- 200
xC <- NULL; xT <- NULL
for (i in 1:nSample){
  xC[i] <- rbinom(1, n, pC[i])
  xT[i] <- rbinom(1, n, pT[i])
}

### Run JAGS: ===============================================

lineinits <- NULL
nBurnIn = 1000
nIter = 10000 + nBurnIn
nmcmc <- nIter - nBurnIn

linedata <- list(xC = xC, xT = xT, n = n, nSample = nSample, a = a, b = b, mu = mu, tau = tau) 
test.sim <- try(jags(data = linedata, inits = lineinits, model.file = test.model,
                     parameters.to.save = c("pT", "pC", "LOR"),
                     n.chains = 1, n.iter = nIter, n.burnin = nBurnIn, n.thin = 1,
                     progress.bar = "text"))
post.pT <- test.sim$BUGSoutput$sims.list$pT #get the posterior distribution sample
post.LOR <- test.sim$BUGSoutput$sims.list$LOR
hist(post.pT)
hist(post.LOR)
### Compute n0 for pT: ===============================================
prepost.pT.BU <- colMeans(post.pT)
hist(prepost.pT.BU)

var.post <- (colSds(post.pT))^2
var.prior <- var(pT)
var.prepost.pT.BU <- var(prepost.pT.BU)
n0 <- n * (var.prior / var.prepost.pT.BU - 1)

# n * mean(var.post) / (var.prior - mean(var.post))
# mean(n * var.post / (var.prior - var.post))

### Compute preposterior GA = f(n0,OR,pC) or just f(n0,pT)?: ===============================================
sqrtv <- sqrt(n/(n + n0))
prepost.pT.GA <- sqrtv * pT + (1-sqrtv) * mean(pT)

# hist(prepost2)
### compare the CDFs for the preposteriors.

cdf.prepost.pT.BU <- ecdf(prepost.pT.BU)
cdf.prepost.pT.GA <- ecdf(prepost.pT.GA)
plot(cdf.prepost.pT.BU, col = "blue")
lines(cdf.prepost.pT.GA, col = "red")

### Compute n0 for LOR: ===============================================
prepost.LOR.BU <- colMeans(post.LOR)
hist(prepost.LOR.BU)

var.LOR.post <- (colSds(post.LOR))^2
var.LOR.prior <- var(LOR)
var.prepost.LOR.BU <- var(prepost.LOR.BU)
n0 <- n * (var.LOR.prior / var.prepost.LOR.BU - 1)

# n * mean(var.post) / (var.prior - mean(var.post))
# mean(n * var.post / (var.prior - var.post))

### Compute preposterior GA = f(n0,OR)?: ===============================================
sqrtv <- sqrt(n/(n + n0))
prepost.LOR.GA <- sqrtv * LOR + (1-sqrtv) * mean(LOR)

# hist(prepost2)
### compare the CDFs for the preposteriors.

cdf.prepost.LOR.bu <- ecdf(prepost.LOR.BU)
cdf.prepost.LOR.GA <- ecdf(prepost.LOR.GA)
plot(cdf.prepost.LOR.bu, col = "blue")
lines(cdf.prepost.LOR.GA, col = "red")
