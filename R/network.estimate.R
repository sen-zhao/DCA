# This function estimates high-dimensional networks, as well as their intersections.
# Author:       Sen Zhao
# Email:        sen-zhao@sen-zhao.com
# ----------------------------------------------------------------------------
# Arguments:
# X:            a list of standardized design matrices with the same number of columns.
# method:       method for network estimation: neighborhood selection ("MB").
# rule:         "AND" or "OR" rule for the estimation.
# lambda:       a vector of tuning parameters.
# ----------------------------------------------------------------------------
# Outputs:
# est:          a list of estimated adjacency matrices.
# common:       adjacency matrix of the common network.
# ----------------------------------------------------------------------------


network.estimate <- function(X, lambda = NULL, method = "MB", rule = "OR"){
  m <- length(X)
  est <- list()
  for(i in 1:m){
    data <- X[[i]]
    n <- nrow(data)
    p <- ncol(data)
    adjacency <- matrix(0, nrow = p, ncol = p)
    for(j in 1:p){
      lammin <- cv.glmnet(data[, -j], data[, j], lambda = lambda)$lambda.min
      b <- glmnet(data[, -j], data[, j], lambda = lammin)$beta
      index <- which(b != 0)
      if(length(index) != 0){
        index[index >= j] <- index[index >= j] + 1
      }
      adjacency[index, j] <- 1
    }
    adjacency <- adjacency + t(adjacency)
    if(rule == "AND"){
      adjacency[adjacency < 2] <- 0
      adjacency[adjacency == 2] <- 1
    }else if(rule == "OR"){
      adjacency[adjacency >= 1] <- 1
    }else{
      stop("Wrong rule. Must be either AND or OR")
    }
    est[[i]] <- adjacency
  }
  common <- est[[1]]
  for(i in 1:m){
    common <- common * est[[i]]
  }
  return(list(est = est, common = common))
}
