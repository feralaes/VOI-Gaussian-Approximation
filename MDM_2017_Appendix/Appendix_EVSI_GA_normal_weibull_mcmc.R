#### normal-weibull MCMC ####
require(R2jags)
shape <- 1 # Weibull shape parameter
n <- 100 # additional data collections
m <- 1000 #number of prior samples
q <- 2000 #number of posterior samples
n.burnin <- 200 # number of burnin posterior values
x <- matrix(0, m, n) # initialize the new data matrix
theta <- rnorm(m, 1, 0.2) # prior distribution of theta
for (i in 1:m) {  # For each prior value, conduct an experiment of size n
  x[i, ] <- rweibull(n, 1, theta[i]) # likelihood: x|p
}

## 2. Obtain the posterior distribution & compute the posterior means
test.model <- function() {
  tau <- 1/pow(0.2,2)
  for (i in 1:m){
    theta[i] ~ dnorm(1, tau) #prior distribution of theta
    lambda[i] <- 1/pow(theta[i], v)
    for (j in 1:n)
    {
      x[i,j] ~ dweib(shape, lambda[i]) # likelihood: x|p
    } 
  }
}
linedata <- list(x = x, n = n, m = m, shape = shape)
test.sim <- try(jags(data = linedata, inits = NULL, model.file = test.model,
                     parameters.to.save = c("theta"),
                     n.chains = 1, n.iter = q + n.burnin, n.burnin = n.burnin, n.thin = 1,
                     progress.bar = "text"))
post.p <- test.sim$BUGSoutput$sims.list$theta #get the posterior distribution sample
prepost <- colMeans(post.p) #obtain the preposterior distribution from the posterior mean of p
var.prior <- var(theta) #prior variance
var.prepost <- var(prepost) #preposterior variance
n0.hat <- n * (var.prior / var.prepost - 1) #predicted prior sample size
n0.hat #24.8