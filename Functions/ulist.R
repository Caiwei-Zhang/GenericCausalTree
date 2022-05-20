level.trt.eff <- function(data.node, x, ux, pscore, cond.outcome) {
  
  # calculate treatment effect of each level of x
  # data.node: a data.frame without column "A"
  # x: current split covariate 
  # ux: the sorted value of x
  
  trt.eff <- numeric(length(ux))
  
  for (i in 1:length(ux)) {
    sub.data.node <- data.node[x == ux[i], ]
    sub.v.hat     <- sub.data.node$A[x == ux[i]] - pscore[x == ux[i]]
    sub.y.star    <- sub.data.node$Y[x == ux[i]] - cond.outcome[x == ux[i]]
    trt.eff[i]    <- sum(sub.v.hat * sub.y.star) / sum(sub.v.hat ^ 2)
  }
  return(trt.eff)
}



# Splitting function
stemp <- function(y, wt, x, parms, continuous) {
  
  # usercode: user-defined splitting function
  # y: vector of response values
  # wt: vector of sample weights
  # x: vector, values of current-splitting covariate 
  # parms: vector of defined parameters, passed forward
  # continuous if TRUE, x should be treated as continuous
  
  n <- length(y)
  
  # extract parameters
  obs.ind      <- y
  X            <- parms$covariates[obs.ind, ]
  A            <- parms$treatment[obs.ind]
  Y            <- parms$outcome[obs.ind] # * wt[obs.ind]
  data.node    <- data.frame(X, A, Y)
  
  pscore        <- parms$pscore[obs.ind]
  pscore.mthd   <- parms$pscore.mthd
  cond.outcome  <- parms$cond.outcome[obs.ind]
  outcome.mthd  <- parms$outcome.mthd
  form.outcome  <- parms$form.outcome
  form.pscore   <- parms$form.pscore
  type.outcome  <- parms$type.outcome
  minsplit      <- parms$minsplit
  minsize       <- parms$minsize
  
  
  if (is.null(pscore) & is.null(cond.outcome)) { # !modelfitBefore
    
    # remove the covariates columns if they have only one level
    data.node.update <- data.node[, c(sapply(data.node[, !colnames(data.node) %in% c("A", "Y")],
                                             function(col) {length(unique(col))}) > 1, TRUE, TRUE)]

    # update formula
    form.outc.r       <- unlist(strsplit(form.outcome, " ~ "))[2]
    form.outc.r.terms <- unlist(strsplit(form.outc.r, " + ", fixed = TRUE))
    form.outc.update  <- paste("Y ~ ", paste0(form.outc.r.terms[form.outc.r.terms %in% colnames(data.node.update)], collapse = " + "), sep = "")
    
    form.psc.r        <- unlist(strsplit(form.pscore, " ~ "))[2]
    form.psc.r.terms  <- unlist(strsplit(form.psc.r, " + ", fixed = TRUE))
    form.psc.update   <- paste("A ~ ", paste0(form.psc.r.terms[form.psc.r.terms %in% colnames(data.node.update)], collapse = " + "), sep = "")
    
    if(nrow(data.node.update) >= 2*minsize*5){ # set lowest crossfit sample size 
      
      pscore       <- est.pscore(form.psc.update, data.node.update, method  = pscore.mthd, folds = 5)
      cond.outcome <- est.outcome(form.outc.update, data.node.update, method = outcome.mthd, type.outcome = "continuous", folds = 5)
      
    } else { # sample size too small to do crossfit
      
      pscore       <- est.pscore(form.psc.update, data.node.update, method  = pscore.mthd, crossfit = FALSE)
      cond.outcome <- est.outcome(form.outc.update, data.node.update, method = outcome.mthd, type.outcome = "continuous", crossfit = FALSE)
      
    }
    
  } 
  
  v.hat  <- A - pscore
  Y.star <- Y - cond.outcome
  
  if (continuous){
    
    # Skip the first and last minsplit splits
    goodness  <- rep(0, n - 1)
    direction <- rep(1, n - 1)
    
    # Ensure the model is at least fitted on minsplit obs
    for (i in (minsplit:(n-minsplit))) { 
      
      # the number of control and treated units must larger than minsize in each node.
      if (min(sum(A[1:i]), sum(A[1:i] == 0)) < minsize | min(sum(A[(i+1):n]), sum(A[(i+1):n] == 0)) < minsize) {
        next
      }
      
      # Calculating the estimator
      trt.eff.l <- sum(v.hat[1:i] * Y.star[1:i]) / sum(v.hat[1:i] ^ 2)
      trt.eff.r <- sum(v.hat[(i+1):n] * Y.star[(i+1):n]) / sum(v.hat[(i+1):n] ^ 2)
    
      delta.trt.eff <- trt.eff.l - trt.eff.r
      
      # calculate the variance estimator i, (n-i)
      var.left  <- sum((Y.star[1:i] - trt.eff.l * v.hat[1:i])^2 * v.hat[1:i]^2) / (sum(v.hat[1:i] ^ 2)) ^ 2
      var.right <- sum((Y.star[(i+1):n] - trt.eff.r * v.hat[(i+1):n])^2 * v.hat[(i+1):n]^2) / (sum(v.hat[(i+1):n] ^ 2)) ^ 2
      
      # Skip this loop if there is negative estimated variance
      if (var.left < 0 | var.right < 0) {
        next
      }
      
      goodness[i] <- delta.trt.eff^2  / (var.left + var.right)
      direction[i] <- c(sign(trt.eff.l - trt.eff.r))      
    }
      
  } else { # !continuous
    
    ux <- sort(unique(x))
    
    # Order the levels of the categorical covariates by their treatment effect
    # Currently code only works for outside node and innode model fitting, not for in split
    trt.eff <- level.trt.eff(data.node, x, ux,
                             pscore         = pscore, # pscore.mthd    = NULL,
                             cond.outcome   = cond.outcome) # outcome.mthd   = NULL)
    
    
    ord.ux <- order(trt.eff)
    goodness  <- rep(0, length(ux) - 1)
    
    # set each level of x as a split point
    for (i in 1:(length(ux) - 1)){
      
      ind.l <- x %in% ux[ord.ux[1:i]]
      ind.r <- x %in% ux[ord.ux[(i+1):length(ux)]] 
      
      if (sum(ind.l) == 1 | sum(ind.r) == 1) {
        next
      }
      
      # the number of control and treated units must larger than minsize in each node.
      if (min(sum(A[ind.l]), sum(A[ind.l] == 0)) < minsize | min(sum(A[ind.r]), sum(A[ind.r] == 0)) < minsize) {
        next
      }
      
      # Calculating the estimator
      trt.eff.l <- sum(v.hat[ind.l] * Y.star[ind.l]) / sum(v.hat[ind.l] ^ 2)
      trt.eff.r <- sum(v.hat[ind.r] * Y.star[ind.r]) / sum(v.hat[ind.r] ^ 2)
      
      delta.trt.eff <- trt.eff.l - trt.eff.r
      
      # calculate the variance estimator length(ind.l) length(ind.r)
      var.left  <- sum((Y.star[ind.l] - trt.eff.l * v.hat[ind.l])^2 * v.hat[ind.l]^2) / (sum(v.hat[ind.l] ^ 2))^2
      var.right <- sum((Y.star[ind.r] - trt.eff.r * v.hat[ind.r])^2 * v.hat[ind.r]^2) / (sum(v.hat[ind.r] ^ 2))^2
      
      # Skip this loop if there is negative estimated variance
      if (var.left < 0 | var.right < 0) {
        next
      }
      
      goodness[i] <- delta.trt.eff^2  / (var.left + var.right)
      
    }
    
    direction <- ux[ord.ux] 
    
  }
  
  list(goodness  = goodness,
       direction = direction)
}




# Evaluation function
etemp <- function(y, wt, parms) {
  
  n <- length(y)
  
  # extract parameters
  obs.ind      <- y
  X            <- parms$covariates[obs.ind, ]
  A            <- parms$treatment[obs.ind]
  Y            <- parms$outcome[obs.ind] # * wt[obs.ind]
  data.node    <- data.frame(X, A, Y)
  
  pscore        <- parms$pscore[obs.ind]
  pscore.mthd   <- parms$pscore.mthd
  cond.outcome  <- parms$cond.outcome[obs.ind]
  outcome.mthd  <- parms$outcome.mthd
  form.outcome  <- parms$form.outcome
  form.pscore   <- parms$form.pscore
  type.outcome  <- parms$type.outcome
  minsize       <- parms$minsize
  
  # get rss (not used)
  wmean <- sum(Y * wt) / sum(wt)
  rss <- sum(wt *(Y-wmean)^2)
  
  # get trt.eff
  if (is.null(pscore) & is.null(cond.outcome)) {

    if (n ==1) {
      return(list(label = NA, deviance = rss))
    }

    data.node.update <- data.node[, c(sapply(data.node[, !colnames(data.node) %in% c("A", "Y")],
                                             function(col) {length(unique(col))}) > 1, TRUE, TRUE)]

    form.outc.r       <- unlist(strsplit(form.outcome, " ~ "))[2]
    form.outc.r.terms <- unlist(strsplit(form.outc.r, " + ", fixed = TRUE))
    form.outc.update  <- paste("Y ~ ", paste0(form.outc.r.terms[form.outc.r.terms %in% colnames(data.node.update)], collapse = " + "), sep = "")
    
    form.psc.r        <- unlist(strsplit(form.pscore, " ~ "))[2]
    form.psc.r.terms  <- unlist(strsplit(form.psc.r, " + ", fixed = TRUE))
    form.psc.update   <- paste("A ~ ", paste0(form.psc.r.terms[form.psc.r.terms %in% colnames(data.node.update)], collapse = " + "), sep = "")
    
    
    if(nrow(data.node.update) >= 2*minsize*5){ # set lowest crossfit sample size 
      
      pscore       <- est.pscore(form.psc.update, data.node.update, method  = pscore.mthd, folds = 5)
      cond.outcome <- est.outcome(form.outc.update, data.node.update, method = outcome.mthd, type.outcome = "continuous", folds = 5)
      
    } else { # sample size too small to do crossfit
      
      pscore       <- est.pscore(form.psc.update, data.node.update, method  = pscore.mthd, crossfit = FALSE)
      cond.outcome <- est.outcome(form.outc.update, data.node.update, method = outcome.mthd, type.outcome = "continuous", crossfit = FALSE)
      
    }
    
  }
  
  v.hat  <- A - pscore
  Y.star <- Y - cond.outcome
  
  trt.eff <- sum(v.hat * Y.star) / sum(v.hat ^ 2)
  #var  <-  sum((Y.star - trt.eff * v.hat)^2 * v.hat^2) / (sum(v.hat ^ 2))^2

  return(list(label = trt.eff, deviance = rss))
  
}




# initialization function
itemp <- function(y, offset, parms, wt) {
  
  # The initialization function needs to perform any data checks, 
  # deal with the offset vector if present, and return a list
  # y: the response value as found in the formula.
  # offset: the offset term, if any, found on the right hand side of the formula
  # parms: the vector or list (if any) supplied by the user as a parms argument to the call.
  # wt: the weight vector from the call, if any
  
  if (is.matrix(y) && ncol(y) > 1)
    stop("Matrix response not allowed")
  # if (!missing(parms) && length(parms) > 0)
  #   warning("parameter argument ignored")
  if (length(offset)) y <- y - offset
  
  sfun <- function(yval, dev, wt, ylevel, digits) {
    paste(" mean=", format(signif(yval, digits)),
          ", MSE=" , format(signif(dev/wt, digits)),
          sep = '')
  }
  
  environment(sfun) <- .GlobalEnv
  list(y = c(y), parms = parms, numresp = 1, numy = 1, summary = sfun)
}


# pscore.mthd  <- "GLM"
# outcome.mthd <- "GBM"
# minsplit <- 30
# min.node <- 40
# 
# data.used.covariates  <- data[, !colnames(data) %in% c("Y", "A")]
# form.used             <- as.formula(paste("rownumb ~ ", paste(colnames(data.used.covariates), collapse= "+")))
# data.used.wiz.rownumb <- cbind(rownumb = 1:dim(data)[1], data)
# data.used.xy <- data[, !colnames(data) %in% c("A")]
# data.used.xa <- data[, !colnames(data) %in% c("Y")]
# 
# ulist.used <- list(eval = etemp, split = stemp, init = itemp)
# 
# # predict pscore and cond.outcome in parent node
# parms.used <- list(covariates = data.used.covariates,
#                    treatment = data$A,
#                    outcome = data$Y,
#                    pscore = NULL, # pscore,
#                    cond.outcome = NULL, # cond.outcome,
#                    # pscore.mthd = pscore.mthd, # NULL,
#                    # outcome.mthd = outcome.mthd, #NULL,
#                    mthd = "GBM",
#                    minsplit = minsplit)
# 
# # pscore, cond.outcomez are predicted before building tree
# # cond.outcome <- est.outcome(data.used.xy, method = "RF")
# # pscore  <- est.pscore(data.used.xa, method = "GLM")
# tmp <- model.fit(data.node = data,
#                  method = "GBM", folds = 5,
#                  type.outcome = "continuous")
# 
# cond.outcome <- tmp$cond.outcome
# pscore       <- tmp$pscore
# 
# parms.used <- list(covariates = data.used.covariates,
#                    treatment = data$A,
#                    outcome = data$Y,
#                    pscore = pscore,
#                    cond.outcome = cond.outcome,
#                    pscore.mthd = NULL,
#                    outcome.mthd = NULL,
#                    minsplit = minsplit)
# 
# # Test
# t0 <- Sys.time()
# a <- rpart(form.used, data = data.used.wiz.rownumb,
#            method   = ulist.used,
#            parms    = parms.used,
#            control  = rpart.control(cp = 0, minsplit = minsplit * 2, minbucket = min.node, 
#                                     maxsurrogate = 0, maxcompete = 0))
# t1 <- Sys.time()
# t1-t0
# 
# tree.standard <- causalTree(Y ~ ., data = data, treatment = data$A, method = "anova", 
#                             cp = 0, minbucket = 20, cv.option = "matching", split.option = "CT", xval = 0)
