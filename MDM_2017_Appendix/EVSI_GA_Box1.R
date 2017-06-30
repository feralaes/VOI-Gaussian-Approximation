rm(list = ls())
#### Step 1 ####
### Load required packages and functions
library(mgcv); library(matrixStats)
source("Rcode/Appendix_Final/GA_functions.R")

### Load PSA
psa <- read.csv("Data/psa_demo_GA.csv")
n.sim <- nrow(psa)
head(psa)

### Create following matrices from your PSA
## Net monetary benefit (NMB) matrix
nmb <- psa[, 5:7]
## Parameter of interest is column 1
theta_I <- psa[, 1]

#### Step 2 ####
### Find optimal strategy (d*) based on the highest expected NMB
d.star <- which.max(colMeans(nmb))

#### Step 3 ####
### Define the Loss matrix
loss <- nmb - nmb[, d.star]

#### Step 4 ####
lmm1 <- gam(loss[, 1] ~ s(theta_I))
lmm2 <- gam(loss[, 2] ~ s(theta_I))
lmm3 <- gam(loss[, 3] ~ s(theta_I))

#### Step 5 ####
## Compute estimated losses
Lhat <- cbind(lmm1$fitted, lmm2$fitted, lmm3$fitted)
### Compute EVPPI on one parameter
evppi <- mean(rowMaxs(Lhat))
evppi

#### Step 6 ####
## Initial sample size of each parameter
n0 <- 10
## Additional sample size
n <- 100

#### Step 7 ####
## Compute predicted loss from the preposterior mean phi
Ltilde1 <- predict.ga(lmm1, n = n, n0 = n0)
Ltilde2 <- predict.ga(lmm2, n = n, n0 = n0)
Ltilde3 <- predict.ga(lmm3, n = n, n0 = n0)
loss.predicted <- cbind(Ltilde1, Ltilde2, Ltilde3)

#### Step 8 ####
## Compute EVSI
evsi <- mean(rowMaxs(loss.predicted))
evsi
