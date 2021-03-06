eval.honest.tree <- function(tree, data.test, eff.modifier, true.trt.eff, trueGroupNum, 
                             splt.point = 0, noise.splt.var, type.eff.mod, split.cate = NULL) {
  
  # split.cate used for categorical splits. 
  # Specify which categories should be in the same node for categorical variables
  splits.var <- tree$frame$var[tree$frame.var != "<leaf>"]
  is.leaf    <- (tree$frame$var ==  "<leaf>")
  num.splt   <- sum(tree$frame$var !=  "<leaf>")
  num.leaf <- sum(tree$frame$var == "<leaf>")
  
  # Measure 1: first split accuracy
  fir.var <- as.character(tree$frame$var[1])
  if (is.null(eff.modifier)) { # homo: fir.splt.ind = NULL or 0 or equals to corr.tree
    fir.splt.ind <- 0 # sum(fir.var == "<leaf>")
  } else {
    fir.splt.ind <- ifelse(fir.var %in% eff.modifier, 1, 0)
    # fir.splt.acc
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
  } else if (num.leaf != trueGroupNum | num.noise.splt != 0) { # hetero # 1.judge if the number of subgroups is correct
    corr.tree <- 0
  } else {
    
    flag <- TRUE
    
    # if (sum(splits.var %in% eff.modifier) == length(splits.var)) { # 2. jugde if all nodes split at effct modifiers 
    #   flag <- flag & TRUE
    # } else {
    #   flag <- flag & FALSE
    # }  
    
    ### 不需要考虑是否所有的分割变量都是effect modifier，只需要看是否存在noise分割；
    ### 如果存在就说明不是corr.tree
    
    if (length(splits.var) == length(unique(splits.var))) { # judge if each effect modifier only split one time 
      flag <- flag & TRUE
    } else {
      flag <- flag & FALSE
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
        flag <- flag & (length(unique(tmp.csplit[split.cate[[1]]])) == 1)
        
        rest.csplit <- tmp.csplit[-split.cate[[1]]]
        rest.csplit <- rest.csplit[rest.csplit != 2]
        flag <- flag & (length(unique(rest.csplit)) == 1)
        
      }
    }
    
    if (flag) {
      corr.tree <- 1
    } else {
      corr.tree <- 0
    }
    
  }
  
  # Measure 5: number of subgroups accuracy
  num.sbg.acc <- ifelse(num.leaf == trueGroupNum, 1, 0)
  
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
  # mse <- sum((est.trt.eff - true.trt.eff) ^ 2) / nrow(data.test)
  if (nrow(tree$frame) > 3) {
    pred.tree = predict(tree, newdata = data.test)
  } else if (nrow(tree$frame) == 3){
    # If there is one or zero splits there is a weird memory error so need to do manually
    pred.tree = rep(NA, nrow(data.test))
    split.used = tree$splits[, 4]
    var.used = tree$frame$var[1]
    col.ind <- which(colnames(data.test) == var.used)
    
    if (length(split.used) > 1) {
      split.used <- split.used[1]
      direction <- tree$splits[1, 2]
    } else {
      direction <- tree$splits[2]
    }
    
    # Need to figure out observations going to the left/right node
    if ((split.used %% 1) == 0){   
      
      # Categorical covariate split
      lvls <- levels(data.test[, col.ind])
      pred.tree[data.test[, col.ind] %in% lvls[tree$csplit[split.used,] == 1]] <- tree$frame$yval[2]
      pred.tree[data.test[, col.ind] %in% lvls[tree$csplit[split.used,] == 3]] <- tree$frame$yval[3]
      
    } else{
      # Continuous covariate split
      # Need to take care of left or right
      if (direction > 0) {
        pred.tree[data.test[,  col.ind] >= split.used] <- tree$frame$yval[2]
        pred.tree[data.test[,  col.ind] < split.used]  <- tree$frame$yval[3]
      } else {
        pred.tree[data.test[,  col.ind] < split.used]  <- tree$frame$yval[2]
        pred.tree[data.test[,  col.ind] >= split.used] <- tree$frame$yval[3]
      }
      
    }
  } else {
    pred.tree = rep(tree$frame$yval, nrow(data.test))
  } 
  
  # mse will be NA if there is only 1 observation in the terminal node
  mse = mean((pred.tree - true.trt.eff)^2, na.rm = T)
  
  # Measure 8: widthCI
  widthCI <- mean(tree$frame$widthCI[tree$frame$var == "<leaf>"])
  
  # Measure 9: coverage of confidence interval
  if (is.null(eff.modifier)) { # homo
    
    true.trt.eff.leaf <- 2
    est.trt.eff.leaf <- tree$frame$yval[1]
    
    bl <- (abs(est.trt.eff.leaf - true.trt.eff.leaf) <= qnorm(0.975) * sqrt(tree$frame$variance[1]))
    cvr <- ifelse(bl == 1, 1, 0)
    
    results <- data.frame(fir.splt.ind   = fir.splt.ind, 
                          abs.err.sp     = abs.err.sp,
                          corr.tree.prop = corr.tree, 
                          num.leaf       = num.leaf,
                          # num.splt       = num.splt,
                          num.noise.splt = num.noise.splt,
                          num.sbg.acc     = num.sbg.acc,
                          pps = pps,  mse  = mse, 
                          cvr = cvr, widthCI = widthCI)
    
    
  } else {  # hetero
    
    if (corr.tree == 1) {
      est.trt.eff.leaf <- tree$frame$yval[tree$frame$var == "<leaf>"]
      
      if (abs(est.trt.eff.leaf[1] - 2) < abs(est.trt.eff.leaf[1] - 5)){
        true.trt.eff.leaf <- c(2, 5)
      } else {
        true.trt.eff.leaf <- c(5, 2)
      }

      bl <- (abs(est.trt.eff.leaf - true.trt.eff.leaf) <= qnorm(0.975) * sqrt(tree$frame$variance[tree$frame$var == "<leaf>"]))
      if(true.trt.eff.leaf[1] == 2) {
        cvr1 <- ifelse(bl[1] == 1, 1, 0)
        cvr2 <- ifelse(bl[2] == 1, 1, 0)
      } else {
        cvr1 <- ifelse(bl[2] == 1, 1, 0)
        cvr2 <- ifelse(bl[1] == 1, 1, 0)
      }
      
    } else {
      cvr1 <- cvr2 <- 0
    }
    
    results <- data.frame(fir.splt.ind   = fir.splt.ind, 
                      abs.err.sp     = abs.err.sp,
                      corr.tree.prop = corr.tree, 
                      num.leaf       = num.leaf,
                      # num.splt       = num.splt,
                      num.noise.splt = num.noise.splt,
                      num.sbg.acc     = num.sbg.acc,
                      pps  = pps,  mse  = mse, 
                      cvr1 = cvr1, cvr2 = cvr2,
                      widthCI = widthCI)
    
    
  }

}
  

# I was wondering how they get the coverage of confidence intervals, too. 
# I guess they calculate the CI in each leaf,  and consider whether the CI covers true individual treatment effect of 
# the observations that fall into this leaf. 
# Simulate for T times, we can get CI's coverage of each observation. Then average? 
### simultaneous ci
# if (corr.tree == 1) {
#   
#   est.trt.eff.leaf <- tree$frame$yval[tree$frame$var == "<leaf>"]
#   if (is.null(eff.modifier)) { #homo
#     true.trt.eff.leaf <- 2
#   } else { # hetero
#     if (abs(est.trt.eff.leaf[1] - 2) < abs(est.trt.eff.leaf[1] - 5)){
#       true.trt.eff.leaf <- c(2, 5)
#     } else {
#       true.trt.eff.leaf <- c(5, 2)
#     }
#     
#     
#   }
#   
#   # bl <- (abs(est.trt.eff.leaf - true.trt.eff.leaf) <= qnorm(1-(1-sqrt(0.9))/2) * sqrt(tree$frame$variance))
#   bl <- (abs(est.trt.eff.leaf - true.trt.eff.leaf) <= qnorm(0.975) * sqrt(tree$frame$variance[tree$frame$var == "<leaf>"]))
#   if(true.trt.eff.leaf[1] == 2) {
#     cvr1 <- ifelse(bl[1] == 1, 1, 0)
#     cvr2 <- ifelse(bl[2] == 1, 1, 0)
#   } else {
#     cvr1 <- ifelse(bl[2] == 1, 1, 0)
#     cvr2 <- ifelse(bl[1] == 1, 1, 0)
#   }
#   
# } else {
#   cvr1 <- cvr2 <- 0
# }
## individual ci
# bl <- (abs(pred.tree - true.trt.eff) <= qnorm(0.95) * sqrt(tree$frame$variance))
# cvr <- ifelse(sum(bl) == trueGroupNum, 1, 0)