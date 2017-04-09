# This function tests differential connectivity in high-dimensional networks.
# Author:             Sen Zhao
# Email:              sen-zhao@sen-zhao.com
# ----------------------------------------------------------------------------
# Arguments:
# X:                  a list of standardized design matrices with the same number of columns.
# est.method:         method for network estimation: neighborhood selection ("MB").
# test.method:        method for hypothesis testing: GraceI ("GraceI"), LDPE ("LDPE"), 
#                     ridge ("ridge") or SKAT ("SKAT").
# sample.split        whether samples need to be randomly splitted for estimation and testing.
# alpha:              alpha level of type-I error rate.
# test.level:         node-wise ("node") or edge-wise ("edge").
# rule:               "AND" or "OR" rule for the estimation.
# ----------------------------------------------------------------------------
# Outputs:
# diffnet:            a matrix or vector of differential connections.

network.test <- function(X, est.method = "MB", test.method = "GraceI", sample.split = FALSE, alpha = 0.05, test.level = "node", rule = "OR"){
  m <- length(X)
  if(!sample.split){
    est <- list()
    for(i in 1:m){
      data <- X[[i]]
      n <- nrow(data)
      p <- ncol(data)
      adjacency <- matrix(0, nrow = p, ncol = p)
      for(j in 1:p){
        lammin <- cv.glmnet(data[, -j], data[, j])$lambda.min
        b <- glmnet(data[, -j], data[, j], lambda = lammin)$beta
        index <- which(b != 0)
        if(length(index) != 0){
          index[index >= j] <- index[index >= j] + 1
        }
        adjacency[index, j] <- 1
      }
      est[[i]] <- adjacency
    }
    common <- est[[1]]
    for(i in 1:m){
      common <- common * est[[i]]
    }
    
    p.thresh <- 1 - (1 - alpha)^(1 / m)
    
    if(test.level == "edge"){
      diffnet <- matrix(0, ncol = p, nrow = p)
      if(test.method == "GraceI"){
        for(i in 1:m){
          data <- X[[i]]
          n <- nrow(data)
          p <- ncol(data)
          for(j in 1:p){
            MJC <- which(common[j, -j] == 0)
            if(length(MJC) != 0){
              Graceres <- graceI.test(data[, j], data[, -j], lambda.2 = exp(seq(from = -5, to = 10, length.out = 30)))
              p.Grace <- Graceres$pvalue[MJC]
              MJC[MJC >= j] <- MJC[MJC >= j] + 1
              diffnet[j, MJC] <- diffnet[j, MJC] + (p.Grace < p.thresh)
            }
          }
        }
      }
      else if(test.method == "LDPE"){
        for(i in 1:m){
          data <- X[[i]]
          n <- nrow(data)
          p <- ncol(data)
          for(j in 1:p){
            MJC <- which(common[j, -j] == 0)
            if(length(MJC) != 0){
              LDPEres <- lasso.proj(data[, -j], data[, j])
              p.LDPE <- LDPEres$pval[MJC]
              MJC[MJC >= j] <- MJC[MJC >= j] + 1
              diffnet[j, MJC] <- diffnet[j, MJC] + (p.LDPE < p.thresh)
            }
          }
        }
      }
      else if(test.method == "ridge"){
        for(i in 1:m){
          data <- X[[i]]
          n <- nrow(data)
          p <- ncol(data)
          for(j in 1:p){
            MJC <- which(common[j, -j] == 0)
            if(length(MJC) != 0){
              ridgeres <- ridge.proj(data[, -j], data[, j])
              p.ridge <- ridgeres$pval[MJC]
              MJC[MJC >= j] <- MJC[MJC >= j] + 1
              diffnet[j, MJC] <- diffnet[j, MJC] + (p.ridge < p.thresh)
            }
          }
        }
      }
      else if(test.method == "SKAT"){
        stop("test.level should be set to node with SKAT")
      }
      diffnet <- diffnet + t(diffnet)
      if(rule == "AND"){
        diffnet[diffnet < 2] <- 0
        diffnet[diffnet == 2] <- 1
      }else if(rule == "OR"){
        diffnet[diffnet >= 1] <- 1
      }else{
        stop("Wrong rule. Must be either AND or OR")
      }
    }
    else if(test.level == "node"){
      diffnet <- rep(0, p)
      if(test.method == "GraceI"){
        for(i in 1:m){
          data <- X[[i]]
          n <- nrow(data)
          p <- ncol(data)
          for(j in 1:p){
            MJC <- which(common[j, -j] == 0)
            if(length(MJC) != 0){
              Graceres <- graceI.test(data[, j], data[, -j], lambda.2 = exp(seq(from = -5, to = 10, length.out = 30)))
              p.Grace <- Graceres$pvalue[MJC]
              p.Grace.adj <- p.adjust(p.Grace, method = "holm")
              diffnet[j] <- diffnet[j] + (sum(p.Grace.adj < p.thresh) > 0)
            }
          }
        }
      }
      else if(test.method == "LDPE"){
        for(i in 1:m){
          data <- X[[i]]
          n <- nrow(data)
          p <- ncol(data)
          for(j in 1:p){
            MJC <- which(common[j, -j] == 0)
            if(length(MJC) != 0){
              LDPEres <- lasso.proj(data[, -j], data[, j])
              p.LDPE <- LDPEres$pval[MJC]
              p.LDPE.adj <- p.adjust(p.LDPE, method = "holm")
              diffnet[j] <- diffnet[j] + (sum(p.LDPE.adj < p.thresh) > 0)
            }
          }
        }
      }
      else if(test.method == "ridge"){
        for(i in 1:m){
          data <- X[[i]]
          n <- nrow(data)
          p <- ncol(data)
          for(j in 1:p){
            MJC <- which(common[j, -j] == 0)
            if(length(MJC) != 0){
              ridgeres <- ridge.proj(data[, -j], data[, j])
              p.ridge <- ridgeres$pval[MJC]
              p.ridge.adj <- p.adjust(p.ridge, method = "holm")
              diffnet[j] <- diffnet[j] + (sum(p.ridge.adj < p.thresh) > 0)
            }
          }
        }
      }
      else if(test.method == "SKAT"){
        for(i in 1:m){
          data <- X[[i]]
          n <- nrow(data)
          p <- ncol(data)
          for(j in 1:p){
            MJC <- which(common[j, -j] == 0)
            if(length(MJC) == p - 1){
              skatmodel <- SKAT_Null_Model(data[, j] ~ 1)
            }else{
              skatmodel <- SKAT_Null_Model(data[, j] ~ data[, -j][, -MJC])
            }
            p.SKAT <- SKAT(data[, -j][, MJC], skatmodel)$p.value
            diffnet[j] <- diffnet[j] + (sum(p.SKAT < p.thresh) > 0)
          }
        }
      }
    }
  }else{
    for(i in 1:m){
      n <- nrow(X[[i]])
      randidx <- sample.int(n, n)
      X[[i]] <- X[[i]][randidx, ]
    }
    est <- list()
    for(i in 1:m){
      data <- X[[i]]
      n <- nrow(data)
      idx <- 1:floor(n / 2)
      p <- ncol(data)
      adjacency <- matrix(0, nrow = p, ncol = p)
      for(j in 1:p){
        lammin <- cv.glmnet(data[idx, -j], data[idx, j])$lambda.min
        b <- glmnet(data[idx, -j], data[idx, j], lambda = lammin)$beta
        index <- which(b != 0)
        if(length(index) != 0){
          index[index >= j] <- index[index >= j] + 1
        }
        adjacency[index, j] <- 1
      }
      est[[i]] <- adjacency
    }
    common <- est[[1]]
    for(i in 1:m){
      common <- common * est[[i]]
    }
    
    p.thresh <- 1 - (1 - alpha)^(1 / m)
    
    if(test.level == "edge"){
      diffnet <- matrix(0, ncol = p, nrow = p)
      if(test.method == "GraceI"){
        for(i in 1:m){
          data <- X[[i]]
          n <- nrow(data)
          idx <- floor(n / 2 + 1):n
          p <- ncol(data)
          for(j in 1:p){
            MJC <- which(common[j, -j] == 0)
            if(length(MJC) != 0){
              Graceres <- graceI.test(data[idx, j], data[idx, -j], lambda.2 = exp(seq(from = -5, to = 10, length.out = 30)))
              p.Grace <- Graceres$pvalue[MJC]
              MJC[MJC >= j] <- MJC[MJC >= j] + 1
              diffnet[j, MJC] <- diffnet[j, MJC] + (p.Grace < p.thresh)
            }
          }
        }
      }
      else if(test.method == "LDPE"){
        for(i in 1:m){
          data <- X[[i]]
          n <- nrow(data)
          idx <- floor(n / 2 + 1):n
          p <- ncol(data)
          for(j in 1:p){
            MJC <- which(common[j, -j] == 0)
            if(length(MJC) != 0){
              LDPEres <- lasso.proj(data[idx, -j], data[idx, j])
              p.LDPE <- LDPEres$pval[MJC]
              MJC[MJC >= j] <- MJC[MJC >= j] + 1
              diffnet[j, MJC] <- diffnet[j, MJC] + (p.LDPE < p.thresh)
            }
          }
        }
      }
      else if(test.method == "ridge"){
        for(i in 1:m){
          data <- X[[i]]
          n <- nrow(data)
          idx <- floor(n / 2 + 1):n
          p <- ncol(data)
          for(j in 1:p){
            MJC <- which(common[j, -j] == 0)
            if(length(MJC) != 0){
              ridgeres <- ridge.proj(data[idx, -j], data[idx, j])
              p.ridge <- ridgeres$pval[MJC]
              MJC[MJC >= j] <- MJC[MJC >= j] + 1
              diffnet[j, MJC] <- diffnet[j, MJC] + (p.ridge < p.thresh)
            }
          }
        }
      }
      else if(test.method == "SKAT"){
        stop("test.level should be set to node with SKAT")
      }
      diffnet <- diffnet + t(diffnet)
      if(rule == "AND"){
        diffnet[diffnet < 2] <- 0
        diffnet[diffnet == 2] <- 1
      }else if(rule == "OR"){
        diffnet[diffnet >= 1] <- 1
      }else{
        stop("Wrong rule. Must be either AND or OR")
      }
    }
    else if(test.level == "node"){
      diffnet <- rep(0, p)
      if(test.method == "GraceI"){
        for(i in 1:m){
          data <- X[[i]]
          n <- nrow(data)
          idx <- floor(n / 2 + 1):n
          p <- ncol(data)
          for(j in 1:p){
            MJC <- which(common[j, -j] == 0)
            if(length(MJC) != 0){
              Graceres <- graceI.test(data[idx, j], data[idx, -j], lambda.2 = exp(seq(from = -5, to = 10, length.out = 30)))
              p.Grace <- Graceres$pvalue[MJC]
              p.Grace.adj <- p.adjust(p.Grace, method = "holm")
              diffnet[j] <- diffnet[j] + (sum(p.Grace.adj < p.thresh) > 0)
            }
          }
        }
      }
      else if(test.method == "LDPE"){
        for(i in 1:m){
          data <- X[[i]]
          n <- nrow(data)
          idx <- floor(n / 2 + 1):n
          p <- ncol(data)
          for(j in 1:p){
            MJC <- which(common[j, -j] == 0)
            if(length(MJC) != 0){
              LDPEres <- lasso.proj(data[idx, -j], data[idx, j])
              p.LDPE <- LDPEres$pval[MJC]
              p.LDPE.adj <- p.adjust(p.LDPE, method = "holm")
              diffnet[j] <- diffnet[j] + (sum(p.LDPE.adj < p.thresh) > 0)
            }
          }
        }
      }
      else if(test.method == "ridge"){
        for(i in 1:m){
          data <- X[[i]]
          n <- nrow(data)
          idx <- floor(n / 2 + 1):n
          p <- ncol(data)
          for(j in 1:p){
            MJC <- which(common[j, -j] == 0)
            if(length(MJC) != 0){
              ridgeres <- ridge.proj(data[idx, -j], data[idx, j])
              p.ridge <- ridgeres$pval[MJC]
              p.ridge.adj <- p.adjust(p.ridge, method = "holm")
              diffnet[j] <- diffnet[j] + (sum(p.ridge.adj < p.thresh) > 0)
            }
          }
        }
      }
      else if(test.method == "SKAT"){
        for(i in 1:m){
          data <- X[[i]]
          n <- nrow(data)
          idx <- floor(n / 2 + 1):n
          p <- ncol(data)
          for(j in 1:p){
            MJC <- which(common[j, -j] == 0)
            if(length(MJC) == p - 1){
              skatmodel <- SKAT_Null_Model(data[idx, j] ~ 1)
            }else{
              skatmodel <- SKAT_Null_Model(data[idx, j] ~ data[idx, -j][, -MJC])
            }
            p.SKAT <- SKAT(data[idx, -j][, MJC], skatmodel)$p.value
            diffnet[j] <- diffnet[j] + (sum(p.SKAT < p.thresh) > 0)
          }
        }
      }
    }
  }
  return(list(diffnet = diffnet))
}
