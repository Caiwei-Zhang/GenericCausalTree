eval.honest.tree <- function(tree, data.test, eff.modifier, true.trt.eff, trueGroupNum, 
                             splt.point = 0, noise.splt.var, type.eff.mod, split.cate = NULL) {
  
  # split.cate used for categorical splits. 
  # Specify which categories should be in the same node for categorical variables
  splits.var <- tree$frame$var[tree$frame.var != "<leaf>"]
  is.leaf    <- (tree$frame$var == "<leaf>")
  num.splt   <- sum(tree$frame$var != "<leaf>")
  num.leaf <- sum(tree$frame$var == "<leaf>")
  
  # Measure 1: first split accuracy
  fir.var <- as.character(tree$frame$var[1])
  if (is.null(eff.modifier)) { # homo: fir.splt.ind = NULL or 0 or equals to corr.tree
    fir.splt.ind <- 0 # sum(fir.var == "<leaf>")
  } else {
    fir.splt.ind <- ifelse(fir.var %in% eff.modifier, 1, 0)
  }
  
  # Measure 2.1: split point absolute error
  if (is.null(type.eff.mod)) {
    abs.err.sp <- 0
  } else if (type.eff.mod == "continuous") {
    if (tree$frame$var[1] == eff.modifier) {
      abs.err.sp <- abs(tree$splits[1, 4] - splt.point)
    } else {
      abs.err.sp <- NA 
    }
  } else { # if eff.modifier not continuous
    abs.err.sp <- 0
  }
  
  # Measure 4: number of noise split
  num.noise.splt  <- sum(tree$frame$var %in% noise.splt.var)
  
  # Measure 3: correct trees (ONLY 1 SPLIT)
  if (is.null(eff.modifier)) { # homo
    corr.tree <- sum(num.splt == 0)
  } else if (num.leaf != trueGroupNum) { # hetero # 1.judge if the number of subgroups is correct
    corr.tree <- 0
  } else {
    
    tmp <- TRUE
    
    if (sum(splits.var %in% eff.modifier) == length(splits.var)) { # 2. jugde if all nodes split at effct modifiers 
      tmp <- tmp & TRUE
    } else {
      tmp <- tmp & FALSE
    }  
    
    if (length(splits.var) == length(unique(splits.var))) { # judge if each effect modifier only split one time 
      tmp <- tmp & TRUE
    } else {
      tmp <- tmp & FALSE
    }
    
    if (!is.null(tree$csplit)) {
      
      # which is the current split
      if (CT) { # in CausalTree, the split.stat doesn't exist; need to use direct way to specify the indices
        split.used <- tree$splits[1, 4]
      } else {
        row.ind <- which(tree$splits[, 3] == tree$frame$split.stat[1])
        split.used <- tree$splits[row.ind, 4]
      }
      
      if (split.used %% 1 == 0) {
        
        tmp.csplit <- tree$csplit[split.used,]
        tmp <- tmp & (length(unique(tmp.csplit[split.cate[[1]]])) == 1)
        
        rest.csplit <- tmp.csplit[-split.cate[[1]]]
        rest.csplit <- rest.csplit[rest.csplit != 2]
        tmp <- tmp & (length(unique(rest.csplit)) == 1)
        
      }
    }
    
    if (tmp) {
      corr.tree <- 1
    } else {
      corr.tree <- 0
    }
    
  }
  
  ####### put test data into the tree
  est.trt.eff <- predict(tree, newdata = data.test)
  
  # Measure 6: pairwise partition similarity
  ind <- 0
  for (i in 1:(nrow(data.test)-1)) {
    for (j in (i+1):nrow(data.test)) {
      Ie <- ifelse(est.trt.eff[i] == est.trt.eff[j], 1, 0)      # subgroup by estimation 
      It <- ifelse(true.trt.eff[i] == true.trt.eff[j], 1, 0)    # true subgroup
      ind <- ind + abs(Ie-It)
    }
  }
  pps <- 1-sum(ind)/choose(nrow(data.test), 2)
  
  
  # Measure 7: mean square error of estimator
  mse <- sum((est.trt.eff - true.trt.eff) ^ 2) / nrow(data.test)
  # if (nrow(tree$frame) > 3) {
  #   pred.tree = predict(tree, newdata = data.test)
  # } else if (nrow(tree$frame) == 3){
  #   pred.tree = rep(NA, nrow(data.test))
  #   split.used = tree$splits[, 4]
  #   var.used = tree$frame$var[1]
  #   col.ind <- which(colnames(data.test) == var.used)
  #   
  #   if (length(split.used) > 1) {
  #     split.used <- split.used[1]
  #     direction <- tree$splits[1, 2]
  #   } else {
  #     direction <- tree$splits[2]
  #   }
  #   
  #   # Need to figure out observations going to the left/right node
  #   if ((split.used %% 1) == 0){   
  #     
  #     # Categorical covariate split
  #     lvls <- levels(data.test[, col.ind])
  #     pred.tree[data.test[, col.ind] %in% lvls[tree$csplit[split.used,] == 1]] <- tree$frame$yval[2]
  #     pred.tree[data.test[, col.ind] %in% lvls[tree$csplit[split.used,] == 3]] <- tree$frame$yval[3]
  #     
  #   } else{
  #     # Continuous covariate split
  #     # Need to take care of left or right
  #     if (direction > 0) {
  #       pred.tree[data.test[,  col.ind] >= split.used] <- tree$frame$yval[2]
  #       pred.tree[data.test[,  col.ind] < split.used]  <- tree$frame$yval[3]
  #     } else {
  #       pred.tree[data.test[,  col.ind] < split.used]  <- tree$frame$yval[2]
  #       pred.tree[data.test[,  col.ind] >= split.used] <- tree$frame$yval[3]
  #     }
  #     
  #   }
  # } else {
  #   pred.tree = rep(tree$frame$yval, nrow(data.test))
  # } 
  # 
  # # mse will be NA if there is only 1 observation in the terminal node
  # = mean((pred.tree - true.trt.eff)^2, na.rm = T)
  
  
  # Measure 9: coverage of confidence interval
  if (is.null(eff.modifier)) { # homo
    
    true.trt.eff.leaf <- 2
    est.trt.eff.leaf <- tree$frame$yval[1]
    
    if (corr.tree == 1) {
      #bl <- (abs(est.trt.eff.leaf - true.trt.eff.leaf) <= qnorm(0.975) * sqrt(tree$frame$variance[1]))
      bl <- (tree$frame$ci.lower[1] <= true.trt.eff.leaf) & (tree$frame$ci.upper[1] >= true.trt.eff.leaf)
      cvr <- ifelse(bl == 1, 1, 0)
      
      # Measure 8: width of confidence interval, remove in Output() if corr.tree == 0
      widthCI <- tree$frame$widthCI[1]
      
    } else {
      
      cvr <- 0
      widthCI <- NA
      
    }
    
    
    # Measure 10: average coverage
    test.frame <- TrueATE(data.test = data.test, tree = tree, true.trt.eff = true.trt.eff)
    avg_cov <- sum((test.frame$ci.lower < test.frame$true.avg.trt.eff) & (test.frame$true.avg.trt.eff < test.frame$ci.upper))/nrow(tree$frame)
    

    results <- data.frame(fir.splt.ind   = fir.splt.ind, 
                          abs.err.sp     = abs.err.sp,
                          corr.tree = corr.tree, 
                          num.leaf       = num.leaf,
                          num.noise.splt = num.noise.splt, 
                          pps  = pps,  mse  = mse, 
                          cvr = cvr, 
                          cov_ratio = avg_cov,
                          widthCI = widthCI) 
  } else {  # hetero 
    
    if (corr.tree == 1) {
      
      est.trt.eff.leaf <- tree$frame$yval[tree$frame$var == "<leaf>"]
      
      if (abs(est.trt.eff.leaf[1] - 2) < abs(est.trt.eff.leaf[1] - 5)){ # left is 2, right is 5
        
        true.trt.eff.leaf <- c(2, 5)
        
      } else {
        
        true.trt.eff.leaf <- c(5, 2)
        
      }
       
      # bl <- abs(est.trt.eff.leaf - true.trt.eff.leaf) <= qnorm(0.975) * sqrt(tree$frame$variance[tree$frame$var == "<leaf>"])
      bl <- (tree$frame$ci.lower[tree$frame$var == "<leaf>"] <= true.trt.eff.leaf) & (tree$frame$ci.upper[tree$frame$var == "<leaf>"] >= true.trt.eff.leaf)
      
      # simultaneous confidence inference: a larger quantile 
      sim_bl <- abs(est.trt.eff.leaf - true.trt.eff.leaf) <= qnorm(1-(1-sqrt(0.95))/2)*sqrt(tree$frame$variance[tree$frame$var == "<leaf>"]) 
      sim_cvr <- ifelse(sum(sim_bl) == 2, 1, 0)
      
      if(true.trt.eff.leaf[1] == 2) {
        
        cvr1 <- ifelse(bl[1] == 1, 1, 0)
        cvr2 <- ifelse(bl[2] == 1, 1, 0)
        
        widthCI1 <- tree$frame$widthCI[tree$frame$var == "<leaf>"][1]
        widthCI2 <- tree$frame$widthCI[tree$frame$var == "<leaf>"][2]
        
      } else {
        
        cvr1 <- ifelse(bl[2] == 1, 1, 0)
        cvr2 <- ifelse(bl[1] == 1, 1, 0)
        
        widthCI1 <- tree$frame$widthCI[tree$frame$var == "<leaf>"][2]
        widthCI2 <- tree$frame$widthCI[tree$frame$var == "<leaf>"][1]
        
      }
      
    } else {
      
      cvr1 <- cvr2 <- 0
      sim_cvr <- 0
      
      widthCI1 <- widthCI2 <- mean(tree$frame$widthCI[tree$frame$var == "<leaf>"])
      
    }
    
    
    # Measure 10: avgerage coverage
    test.frame <- TrueATE(data.test = data.test, tree = tree, true.trt.eff = true.trt.eff)
    avg_cov <- sum((test.frame$ci.lower < test.frame$true.avg.trt.eff) & (test.frame$true.avg.trt.eff < test.frame$ci.upper))/nrow(tree$frame)
    
    results <- data.frame(fir.splt.ind   = fir.splt.ind, 
                      abs.err.sp = abs.err.sp,
                      corr.tree = corr.tree, 
                      num.leaf  = num.leaf,
                      num.noise.splt = num.noise.splt,
                      pps  = pps,  mse  = mse, 
                      cvr1 = cvr1, cvr2 = cvr2, 
                      sim_cvr = sim_cvr, 
                      cov_ratio = avg_cov,
                      widthCI1 = widthCI1, widthCI2 = widthCI2)
  }
}
  
