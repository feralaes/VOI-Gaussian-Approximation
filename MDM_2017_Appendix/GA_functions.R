predict.ga <- function(object, n, n0, verbose = T){
  #### Function to compute the preposterior for each of the 
  #### basis functions of the GAM model.
  #### Inputs:
  #### - object: gam object
  #### - n: scalar or vector of new sample size to compute evsi on
  #### - n0: scalar or vector of effective prior sample size
  #### - verbose: Prints the variance reduction factor for each parameter
  
  ### Name of parameters
  names.data <- colnames(object$model)
  ### Create dataframe with parameter values
  data <- data.frame(object$model[,-1])
  ## Name columns of dataframe 
  colnames(data) <- names.data[-1]
  
  ### Number of parameters
  n.params <- ncol(data)
  
  ### Sanity checks
  if(!(length(n)==1 | length(n)==n.params)){
    stop("Variable 'n' should be either a scalar or a vector 
         the same size as the number of parameters")
  }
  if(!(length(n0)==1 | length(n0)==n.params)){
    stop("Variable 'n0' should be either a scalar or a vector 
         the same size as the number of parameters")
  }
  
  ### Make n & n0 consistent with the number of parameters
  if(length(n) == 1){
    n <- rep(n, n.params)
  }
  if(length(n0) == 1){
    n0 <- rep(n0, n.params)
  }
   
  ### Compute variance reduction factor
  v.ga <- sqrt(n/(n+n0))
  if (verbose){
    print(paste("Variance reduction factor =", round(v.ga, 3)))
  }
  
  ### Number of smoothers
  n.smooth <- length(object$smooth)
  ### Number of total basis functions
  n.colX <- length(object$coefficients)
  ### Number of observations 
  n.rowX <- nrow(object$model)
  
  ### Initialize matrix for preposterior of total basis functions
  X <- matrix(NA, n.rowX, n.colX)
  X[, 1] <- 1
  
  for (k in 1:n.smooth) { # k <- 1
    klab <- substr(object$smooth[[k]]$label, 1, 1)
    if (klab == "s"){
      Xfrag <- Predict.smooth.ga(object$smooth[[k]], data, v.ga[k])
    } else {
      Xfrag <- Predict.matrix.tensor.smooth.ga(object$smooth[[k]], data, v.ga)
    }
    X[, object$smooth[[k]]$first.para:object$smooth[[k]]$last.para] <- Xfrag
  }
  
  ### Coefficients of GAM model
  Beta <- coef(object)
  
  ### Compute conditional Loss
  Ltilde <- X %*% Beta
  
  return(Ltilde)
}

Predict.smooth.ga <- function (object, data, v.ga = 1) {
  #### Function to compute the preposterior for each of the 
  #### basis functions of a smooth for one parameter
  
  ### Produce basis functions for one parameter
  X <- PredictMat(object, data)
  ## Number of observations
  n.obs <- nrow(X)
  
  ### Apply variance reduction to compute the preposterior 
  ### for each of the basis functions
  ## Vector of ones
  ones <- matrix(1, n.obs, 1)
  ## Compute phi on each of the basis function
  X.ga <- v.ga*X + (1-v.ga)*(ones %*% colMeans(X))
  
  return(X.ga)
}

Predict.matrix.tensor.smooth.ga <- function (object, 
                                             data, 
                                             v.ga = rep(1, ncol(data))){
  #### Function to compute the preposterior for each of the 
  #### basis functions for one or more parameters and calculates
  #### the tensor product if more than one parameter is selected
  #### (Heavily based on function Predict.matrix.tensor.smooth from
  #### mgcv package)
  
  m <- length(object$margin)
  X <- list()
  for (i in 1:m) { # i <- 1
    term <- object$margin[[i]]$term
    dat <- list()
    for (j in 1:length(term)) { # j <- 1
      dat[[term[j]]] <- data[[term[j]]]
    }
    X[[i]] <- if (!is.null(object$mc[i])) # before: object$mc[i]
      PredictMat(object$margin[[i]], dat, n = length(dat[[1]]))
    else Predict.matrix(object$margin[[i]], dat)
    n.obs <- nrow(X[[i]])
  } # end for 'i'
  mxp <- length(object$XP)
  if (mxp > 0) 
    for (i in 1:mxp) if (!is.null(object$XP[[i]])) 
      X[[i]] <- X[[i]] %*% object$XP[[i]]
  
  ### Apply variance reduction to compute the preposterior 
  ### for each of the basis functions
  ## Vector of ones
  ones <- matrix(1, n.obs, 1)
  ## Initialize and fill list with preposterior of basis functions 
  ## for each parameter
  X.ga <- list()
  for (i in 1:m) { # i <- 1
    X.ga[[i]] <- v.ga[i]*X[[i]] + (1-v.ga[i])*(ones %*% colMeans(X[[i]]))
  }
  
  ### Compute tensor product
  T.ga <- tensor.prod.model.matrix(X.ga)
  
  return(T.ga)
}