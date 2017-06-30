rm(list = ls())
### Load PSA
psa <- read.csv("Data/psa_demo_GA.csv")
head(psa)

### Create following matrices from your PSA
## Net monetary benefit (NMB) matrix
nmb <- psa[, 5:7]
head(nmb)
## Matrix of model parameter inputs values theta
theta <- psa[, 1:4]
head(theta)
## Number of simulations
n.sim        <- nrow(nmb)
## Number of strategies
n.strategies <- ncol(nmb)

### Load required packages and functions
## For column and row stats
library(matrixStats)  
## To fit spline models
library(mgcv)
## Functions to calculate the conditional loss by computing the 
## preposterior of each of the basis functions of the GAM model
source("Rcode/Appendix_Final/GA_functions.R")

### Find optimal strategy (d*) based on the highest expected NMB
d.star <- which.max(colMeans(nmb))
d.star

### Define the Loss matrix
loss <- nmb - nmb[, d.star]

### EVPI
evpi <- mean(rowMaxs(as.matrix(loss)))
evpi

#========================#
#### Single parameter ####
#========================#
### Generate linear metamodel of one parameter for each opportunity loss
## Selected parameter for EVPPI & EVSI
sel.param <- 3 
lmm1 <- gam(as.formula(paste("loss[, 1] ~ s(", colnames(theta)[sel.param], ")")),
            data = theta)
lmm2 <- gam(as.formula(paste("loss[, 2] ~ s(", colnames(theta)[sel.param], ")")),
            data = theta)
lmm3 <- gam(as.formula(paste("loss[, 3] ~ s(", colnames(theta)[sel.param], ")")),
            data = theta)


#### Compute EVPPI on one parameter ####
## Compute estimated losses
loss.hat <- cbind(lmm1$fitted, lmm2$fitted, lmm3$fitted)

### Apply EVPPI equation
evppi <- mean(rowMaxs(loss.hat))
evppi

#### Compute EVSI on one parameter ####
## Initial sample size
n0 <- 10 
## Additional sample size
n <- 100

### Compute expected conditional loss for each strategy
Ltilde1 <- predict.ga(lmm1, n = n, n0 = n0)
Ltilde2 <- predict.ga(lmm2, n = n, n0 = n0)
Ltilde3 <- predict.ga(lmm3, n = n, n0 = n0)
## Combine losses into one matrix
loss.tilde <- cbind(Ltilde1, Ltilde2, Ltilde3)

### Apply EVSI equation
evsi <- mean(rowMaxs(loss.tilde))
evsi

#======================#
#### Two parameters ####
#======================#
### Generate linear metamodel of two parameters for each opportunity loss
## Select parameters for EVPPI & EVSI
sel.param1 <- 1
sel.param2 <- 3
sel.params <- c(sel.param1, sel.param2)
### Estimate linear metamodel of two parameters
lmm1 <- gam(as.formula(paste("loss[, 1] ~ s(",
                             paste(colnames(theta[, sel.params]), collapse= ") + s("), 
                             ") + ti(",
                             paste(colnames(theta[, sel.params]), collapse= ", "),
                             ")")), data = theta)
lmm2 <- gam(as.formula(paste("loss[, 2] ~ s(",
                             paste(colnames(theta[, sel.params]), collapse= ") + s("), 
                             ") + ti(",
                             paste(colnames(theta[, sel.params]), collapse= ", "),
                             ")")), data = theta)
lmm3 <- gam(as.formula(paste("loss[, 3] ~ s(",
                             paste(colnames(theta[, sel.params]), collapse= ") + s("), 
                             ") + ti(",
                             paste(colnames(theta[, sel.params]), collapse= ", "),
                             ")")), data = theta)

#### Compute EVPPI on two parameters ####
## Compute estimated losses
loss.hat <- cbind(lmm1$fitted, lmm2$fitted, lmm3$fitted)

### Apply EVPPI equation
evppi <- mean(rowMaxs(loss.hat))
evppi

#### Compute EVSI on two parameters Case 1: Same n0 and n ####
## Initial sample size of each parameter
n0 <- c(10, 10)
## Additional sample size
n  <- c(100, 100)

### Compute expected conditional loss for each strategy
Ltilde1 <- predict.ga(lmm1, n = n, n0 = n0)
Ltilde2 <- predict.ga(lmm2, n = n, n0 = n0)
Ltilde3 <- predict.ga(lmm3, n = n, n0 = n0)
## Combine losses into one matrix
loss.tilde <- cbind(Ltilde1, Ltilde2, Ltilde3)

### Apply EVSI equation
evsi <- mean(rowMaxs(loss.tilde))
evsi

#### Compute EVSI on two parameters Case 2: Different n0, same n ####
## Initial sample size of each parameter
n0 <- c(10, 50)
## Additional sample size
n  <- c(100, 100)

### Compute expected conditional loss for each strategy
Ltilde1 <- predict.ga(lmm1, n = n, n0 = n0)
Ltilde2 <- predict.ga(lmm2, n = n, n0 = n0)
Ltilde3 <- predict.ga(lmm3, n = n, n0 = n0)
## Combine losses into one matrix
loss.tilde <- cbind(Ltilde1, Ltilde2, Ltilde3)

### Apply EVSI equation
evsi <- mean(rowMaxs(loss.tilde))
evsi

#### Compute EVSI on two parameters Case 3: Different n0 and n ####
## Initial sample size of each parameter
n0 <- c(10, 50)
## Additional sample size
n  <- c(100, 1000)

### Compute expected conditional loss for each strategy
Ltilde1 <- predict.ga(lmm1, n = n, n0 = n0)
Ltilde2 <- predict.ga(lmm2, n = n, n0 = n0)
Ltilde3 <- predict.ga(lmm3, n = n, n0 = n0)
## Combine losses into one matrix
loss.tilde <- cbind(Ltilde1, Ltilde2, Ltilde3)

### Apply EVSI equation
evsi <- mean(rowMaxs(loss.tilde))
evsi