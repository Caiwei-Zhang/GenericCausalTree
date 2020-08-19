eval.tree <- function(tree, data.test, eff.modifier, true.trt.eff, trueGroupNum, 
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
  
  # Measure 4: number of noise split
  num.noise.splt  <- sum(tree$frame$var %in% noise.splt.var)
  
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
  
  
  
  results <- data.frame(fir.splt.ind = fir.splt.ind, 
                        abs.err.sp    = abs.err.sp,
                        corr.tree.prop= corr.tree, 
                        num.leaf       = num.leaf,
                        num.splt       = num.splt,
                        num.noise.splt = num.noise.splt,
                        num.sbg.acc     = num.sbg.acc,
                        pps = pps, mse  = mse)
}





# Judge correct tree
# tmp <- TRUE
# 
# for (svar in splits.var) {
#   
#   if(svar %in% eff.modifier) {
#     tmp <- tmp & TRUE
#   } else{
#     tmp <- tmp & FALSE
#     break
#   }
#   
#   if (sum(eff.modifier %in% splits.var) == 1) {
#     tmp <- tmp & TRUE
#   }
#   
# }



# if (num.leaf != trueGroupNum){ 
#   corr.tree <- 0
#   
#   # Measure 5: number of subgroups accuracy
#   num.group.acc <- 0
#   
# } else { # num.leaf == trueGroupNum
#   
#   if(!is.null(tree$csplit)) { # binary or categorical
#     
#     tree$cspli
#     
#   } else { # continuous
#     
#     if (sum(!tree$frame$var %in% eff.modifier) != 1) {
#       corr.tree <- 0
#     } 
#     
#     # judge if the split times of eff.modifier is 1
#     for (em in eff.modifier) {
#       if (sum(tree$frame$var == em) != 1) {
#         corr.tree <- 0
#       } 
#     }
#     
#   }
#   
# }
