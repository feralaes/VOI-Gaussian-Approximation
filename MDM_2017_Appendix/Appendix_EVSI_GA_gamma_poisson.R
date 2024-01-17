#### Gamma-Poisson ####
### 1. Obtain samples of theta from the prior distribution in R
a <- 50 # alpha parameter for the gamma prior distribution theta~gamma(a,b)
b <- 0.01  # beta parameter for the gamma distribution
nSample <- 10000 # Number of samples from the prior
## Draw sample
theta <- rgamma(nSample, shape = a, scale = b)
## Compute variance of theta
var.theta <- var(theta)
### 2. Obtain experimental data from the data likelihood in R
S <- numeric(nSample) # Initialize summary statistic S vector
n <- 50 # Additional data collections
## Generate data and compute summary statistic S
for (i in 1:nSample) {
  k <- rpois(n, theta[i])
  S[i] <- mean(k)
}
## Compute variance of S
var.S <- var(S)
## Estimate n0
n0.hat <- n*(var.S/var.theta - 1)
n0.hat