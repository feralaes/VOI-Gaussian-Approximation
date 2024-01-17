#### beta-binomial ####
### 1. Obtain samples of theta from the prior distribution in R
a <- 12 # alpha parameter for the beta prior distribution theta~beta(a,b)
b <- 8  # beta parameter for the beta distribution
nSample <- 1000 # Number of samples from the prior
## Draw sample
theta <- rbeta(nSample, a, b)
## Compute variance of theta
var.theta <- var(theta)
### 2. Obtain experimental data from the data likelihood in R
S <- numeric(nSample) # Initialize summary statistic S vector
n <- 50 # Additional data collections
## Generate data and compute summary statistic S
for (i in 1:nSample) {
  k <- rbinom(1, n, theta[i])
  S[i] <- k/n
}
## Compute variance of S
var.S <- var(S)
## Estimate n0
n0.hat <- n*(var.S/var.theta - 1)
n0.hat
