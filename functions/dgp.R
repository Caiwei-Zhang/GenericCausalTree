# without confounding 
# randomized trial
no.conf <- function(N, p = 5, dgp = c("iid", "cov", "rand", "mix"), N.test = 1000, p.cont) {
    
    Sigma <- matrix(0, p, p)
    Sigma <- ifelse(col(Sigma)==row(Sigma), 1, 0.3)
    X <- mvrnorm(N, rep(0, p), Sigma)
    colnames(X) <- paste0("X", 1:p, sep="")
    for (i in 1:p) {assign(paste0("X", i, sep=""), X[,i])}
    
    e <- 0.5
    A <- rbinom(N, 1, e)
    Y <- 2 + 2*A + X1^2 + exp(X2) + 2*X3 + 3*A*ifelse(X4>0,1,0) + X5^3 + rnorm(N)
    
    X.test <- matrix(rnorm(N.test*p), N.test, p)
    colnames(X.test) <- paste0("X", 1:p, sep="")
    
    data <- data.frame(X, Y, A)
    data.test <- data.frame(X.test)
    
    true.trt.eff <- ifelse(data.test[, 4] > 0, 5, 2)
    
    out <- list(data = data, data.test = data.test,
                true.trt.eff = true.trt.eff, trueGroupNum = 2, eff.modifier = "X4", 
                type.eff.mod = "continuous", splt.point = 0, 
                noise.splt.var = paste("X", c(1:3,5), sep = ""))
    
}


## with confounding 
# heterogeneous treatment effect
hetero.DataGen <- function(N, p = 5, dgp = c("iid", "cov", "mix"), N.test = 1000, p.cont) {
  
  if (dgp == "iid") { # covariates iid generate from N(0, 1)
    
    X <- matrix(rnorm(N*p), N, p)
    colnames(X) <- paste0("X", 1:p, sep="")
    for (i in 1:p) {assign(paste0("X", i, sep=""), X[,i])}
    
    pscore.coeff <- matrix(c(0.5, -0.5, 0.5, -0.5, 0.5), 5, 1)
    e <- exp(X[, 1:5] %*% pscore.coeff)/(1 + exp(X[, 1:5] %*% pscore.coeff))
    A <- rbinom(N, 1, e)
    
    Y <- 2 + 2*A + X1^2 + exp(X2) + 2*X3 + 3*A*ifelse(X4>0,1,0) + X5^3 + rnorm(N)
    # Y <- 2 + 2*A + X1^2 + X1 * X2 + 2*X3 + 3*A*ifelse(X4>0,1,0) + X5^3 + rnorm(N)
    
    X.test <- matrix(rnorm(N.test*p), N.test, p)
    colnames(X.test) <- paste0("X", 1:p, sep="")
    
    true.trt.eff <- ifelse(X.test[, 4] > 0, 5, 2)
    type.eff.mod <- "continuous"
    
  } 
  
  if (dgp == "cov") { # covariates generate from mvnormal distribution with covariance 0.3
    
    Sigma <- matrix(0, p, p)
    Sigma <- ifelse(col(Sigma)==row(Sigma), 1, 0.3)
    X <- mvrnorm(N, rep(0, p), Sigma)
    colnames(X) <- paste0("X", 1:p, sep="")
    for (i in 1:p) {assign(paste0("X", i, sep=""), X[,i])}
    
    pscore.coeff <- matrix(c(0.5, -0.5, 0.5, -0.5, 0.5), 5, 1)
    e <- exp(X[, 1:5] %*% pscore.coeff)/(1 + exp(X[,1:5] %*% pscore.coeff))
    A <- rbinom(N, 1, e)
    Y <- 2 + 2*A + X1^2 + exp(X2) + 2*X3 + 3*A*ifelse(X4>0,1,0) + X5^3 + rnorm(N)

    X.test <- matrix(rnorm(N.test*p), N.test, p)
    colnames(X.test) <- paste0("X", 1:p, sep="")
    
    true.trt.eff <- ifelse(X.test[, 4] > 0, 5, 2)
    type.eff.mod <- "continuous"
    
  }
  
  
  if (dgp == "mix") { # three covariates iid generate from N(0, 1), one covariates is binary, one is category
    
    p.cont <- 3
    
    X.cont <- matrix(rnorm(N*p.cont), N, p.cont)
    X.bin  <- rbinom(N, 1, 0.5)
    X.cate <- rbinom(N, 2, 0.5)
    X <- cbind(X.cont, X.bin, X.cate)
    
    colnames(X) <- paste0("X", 1:p, sep="")
    for (i in 1:p) {assign(paste0("X", i, sep=""), X[,i])}
    
    pscore.coeff <- matrix(c(0.5, -0.5, 0.5, -0.5, 0.5), 5, 1)
    e <- exp(X[, 1:5] %*% pscore.coeff)/(1 + exp(X[, 1:5] %*% pscore.coeff))
    A <- rbinom(N, 1, e)
    
    Y <- 2 + 2*A + X1^2 + exp(X2) + 2*X3 + 3*A*ifelse(X4==0,1,0) + ifelse(X5==1, 1, 0) + rnorm(N)
    
    X.test <- cbind(matrix(rnorm(N.test*p.cont), N.test, p.cont), 
                    rbinom(N.test, 1, 0.5), 
                    rbinom(N.test, 2, 0.5))
    colnames(X.test) <- paste0("X", 1:p, sep="")
    
    true.trt.eff <- ifelse(X.test[, 4] == 0, 5, 2)
    type.eff.mod <- "catergory"
  }
  
  data <- data.frame(X, Y, A)
  data.test <- data.frame(X.test)
  
    
  out <- list(data = data, data.test = data.test, true.trt.eff = true.trt.eff,
              trueGroupNum = 2, eff.modifier = "X4", type.eff.mod = type.eff.mod,
              splt.point = 0, noise.splt.var = paste("X", c(1:3,5), sep = ""))
}
  



## with confounding 
# homogeneous treatment effect
homo.DataGen <- function(N, p = 5, dgp = c("iid", "cov", "mix"), N.test = 1000, p.cont) {
  
  if (dgp == "iid") { # covariates iid generate from N(0, 1)
    
    X <- matrix(rnorm(N*p), N, p)
    colnames(X) <- paste0("X", 1:p, sep="")
    for (i in 1:p) {assign(paste0("X", i, sep=""), X[,i])}
    
    pscore.coeff <- matrix(c(0.5, -0.5, 0.5, -0.5, 0.5), 5, 1)
    e <- exp(X[, 1:5] %*% pscore.coeff)/(1 + exp(X[, 1:5] %*% pscore.coeff))
    A <- rbinom(N, 1, e)
    
    Y <- 2 + 2*A + X1^2 + exp(X2) + 2*X3 + 3*ifelse(X4>0,1,0) + X5^3 + rnorm(N)
    
    X.test <- matrix(rnorm(N.test*p), N.test, p)
    colnames(X.test) <- paste0("X", 1:p, sep="")
  
    
  } 
  
  if (dgp == "cov") { # covariates generate from mvnormal distribution with covariance 0.3
    
    Sigma <- matrix(0, p, p)
    Sigma <- ifelse(col(Sigma)==row(Sigma), 1, 0.3)
    X <- mvrnorm(N, rep(0, p), Sigma)
    colnames(X) <- paste0("X", 1:p, sep="")
    for (i in 1:p) {assign(paste0("X", i, sep=""), X[,i])}
    
    pscore.coeff <- matrix(c(0.5, -0.5, 0.5, -0.5, 0.5), 5, 1)
    e <- exp(X[, 1:5] %*% pscore.coeff)/(1 + exp(X[, 1:5] %*% pscore.coeff))
    A <- rbinom(N, 1, e)
    Y <- 2 + 2*A + X1^2 + exp(X2) + 2*X3 + 3*ifelse(X4>0,1,0) + X5^3 + rnorm(N)
    
    X.test <- matrix(rnorm(N.test*p), N.test, p)
    colnames(X.test) <- paste0("X", 1:p, sep="")
    
  }
  
  if (dgp == "mix") { # three covariates iid generate from N(0, 1), one covariates is binary, one is category
    
    p.cont <- 3
    
    X.cont <- matrix(rnorm(N*p.cont), N, p.cont)
    X.bin  <- rbinom(N, 1, 0.5)
    X.cate <- rbinom(N, 2, 0.5)
    X <- cbind(X.cont, X.bin, X.cate)
    
    colnames(X) <- paste0("X", 1:p, sep="")
    for (i in 1:p) {assign(paste0("X", i, sep=""), X[,i])}
    
    pscore.coeff <- matrix(c(0.5, -0.5, 0.5, -0.5, 0.5), 5, 1)
    e <- exp(X[, 1:5] %*% pscore.coeff)/(1 + exp(X[, 1:5] %*% pscore.coeff))
    A <- rbinom(N, 1, e)
    
    Y <- 2 + 2*A + X1^2 + exp(X2) + 2*X3 + 3*ifelse(X4==0,1,0) + ifelse(X5==1, 1, 0) + rnorm(N)
    
    X.test <- cbind(matrix(rnorm(N.test*p.cont), N.test, p.cont), 
                    rbinom(N.test, 1, 0.5), 
                    rbinom(N.test, 2, 0.5))
    colnames(X.test) <- paste0("X", 1:p, sep="")
    
  }
  
  data <- data.frame(X, Y, A)
  data.test <- data.frame(X.test)
  true.trt.eff <- rep(2, N.test)
  
  out <- list(data = data, data.test = data.test, true.trt.eff = true.trt.eff,
              trueGroupNum = 1, eff.modifier = NULL, type.eff.mod = NULL,
              splt.point = NULL, noise.splt.var = paste("X", 1:p, sep = ""))
}


## noise setting: increase dimensions of noise covariates.