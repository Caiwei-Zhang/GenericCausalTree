honest.Treebuild <- function(data, form.outcome, form.pscore, weights,  # subset, x, y,
                      type.outcome, train.frac = 0.8, # na.action = na.rpart,
                      modelFitBefore = FALSE, est.modelFitBefore = FALSE,
                      outcome.mthd, pscore.mthd, cp,
                      minsplit = 30, minsize = 10) {
  
  # data.train: training data used to build the tree
  # form.outcome: the form of outcome model, contains outcome and covariates but no treatment.
  # form.pscore: the form of propensity score model, contains treatment and covariates but no outcome.
  # treatment: specify the treatment variable
  # weights: sample weights in tree building. If NULL, each observation receives the same weight.
  # type.outcome: the class of outcome, transform by "method" in rpart.
  # train.frac: the fraction of observations used to build a tree and do honest estimation.
  # modeFitBefore: binary, predicted conditional mean of outcome and propensity score before tree building (TRUE) or not (FALSE).
  # est.modelFitBefore: binary, fit models before put data.est into tree, or fit models in nodes after put data.est into tree. 
  # outcome.mthd: machine learning method used to estimate conditional mean of outcome.
  # pscore.mthd: machine learning method used to estimate propensity score.
  # minsplit: minimum sample size when a split is attempted.
  # minsize: minimum sample size in terminal <leaf> node.
  
  data.tr   <- data[1:(nrow(data)/2), ]
  train.idx <- 1 : (train.frac * nrow(data.tr))
  data.train <- data.tr[train.idx, ]
  data.valid <- data.tr[-train.idx, ]
  
  data.est   <- data[(nrow(data)/2 + 1) : nrow(data), ]
  
  ## get formula
  if (missing(form.outcome) | missing(form.pscore)) { 
    
    form.outcome <- "Y ~ ."
    form.pscore  <- "A ~ ."
    
  }
  
  if(form.outcome == "Y ~ ." | form.pscore == "A ~ .") {
    
    form.r.terms <- colnames(data.train)[!colnames(data.train) %in% c("Y", "A")]
    
    form.outcome <- paste("Y ~ ", paste0(form.r.terms, collapse = " + "), sep = "")
    form.pscore  <- paste("A ~ ", paste0(form.r.terms, collapse = " + "), sep = "")
    
  }
  
  data.covariates  <- data.train[, !colnames(data.train) %in% c("Y", "A")]
  form.rpart       <- as.formula(paste("rownumb ~ ", paste(colnames(data.covariates), collapse = "+"))) 
  data.rownumb     <- cbind(rownumb = 1:nrow(data.train), data.train)
  
  ## get type.outcome
  if (length(levels(as.factor(data.train$Y) == 2))) {
    type.outcome <- "binary"
  } else if(is.matrix(data.train$Y)) {
    stop("outcome cannot be matrix.")
  } else if (inherits(data.train$Y, "Surv")) {
    stop("outcome cannot be survival object.")
  } else {
    type.outcome <- "continuous"
  }
  
  
  ## set outcome.mthd and pscore.mthd
  if (missing(outcome.mthd)) {
    print("Missing outcome.mthd, GBM will be used for fitting outcome model.")
    outcome.mthd <- "GBM"
  }
  
  if (missing(pscore.mthd)) {
    print("Missing pscore.mthd, GLM will be used for fitting propensity score model.")
    pscore.mthd <- "GLM"
  }
  
  
  ## set parameter lists of two nuisance-fitting patterns
  if (modelFitBefore) {
    
    cond.outcome  <- est.outcome(form.outcome, data.train, method = outcome.mthd, folds = 5)
    pscore        <- est.pscore(form.pscore, data.train, method = pscore.mthd, folds = 5)
    
    parms <- list(covariates   = data.covariates,
                  treatment    = data.train$A,
                  outcome      = data.train$Y,
                  est.data     = data.est,
                  pscore       = pscore,
                  cond.outcome = cond.outcome,
                  form.outcome = NULL,
                  form.pscore  = NULL,
                  pscore.mthd  = NULL,
                  outcome.mthd = NULL,
                  minsplit     = minsplit,
                  minsize      = minsize)
    
    
  } else { # modelFitBefore = FALSE,
      
    parms <- list(covariates   = data.covariates,
                  treatment    = data.train$A,
                  outcome      = data.train$Y,
                  est.data     = data.est,
                  pscore       = NULL,
                  cond.outcome = NULL, 
                  form.outcome = form.outcome,  # form must transport in "character" type, not formula
                  form.pscore  = form.pscore,  
                  pscore.mthd  = pscore.mthd,
                  outcome.mthd = outcome.mthd, 
                  minsplit     = minsplit,
                  minsize      = minsize)    
    
  }
  
  # Make CART split using user-written splitting function (Package rpart)
  ulist <- list(eval = etemp, split = stemp, init = itemp)    
  
  # generate initial tree
  init.tree <- rpart(form.rpart, data = data.rownumb,
                     method  = ulist, parms = parms,
                     control = rpart.control(cp = 0, minsplit = minsplit * 2, 
                                             minbucket = minsize, 
                                             maxsurrogate = 0, 
                                             maxcompete = 0))
  
  # add num.leaf, leaf.size and leaf.obs to cptree object. 
  init.tree$num.leaf   <- sum(init.tree$frame$var == "<leaf>")
  init.tree$leaf.size  <- init.tree$frame$n[init.tree$frame$var == "<leaf>"]
  
  for (g in 1:init.tree$num.leaf) {
    init.tree$leaf.obs[[g]] <- names(init.tree$where)[init.tree$where == which(init.tree$frame$var == "<leaf>")[g]]
  }
  
  
  # if init.tree has only one node
  if (nrow(init.tree$frame) == 1) {
    
    is.leaf <- TRUE
    tmp.frame <- init.tree$frame %>% mutate(delta = 0) %>% mutate(g.h = Inf)
    init.tree$frame <- tmp.frame
    init.tree$num.leaf <- 1
    init.tree$leaf.size <- nrow(data.train)
    init.tree$leaf.obs <- 1:nrow(data.train)
    
  } else {
    
    # judge if a node in row of the frame is a leaf 
    is.leaf <- (init.tree$frame$var == "<leaf>")
    
    # add delta(goodness) of each split to the frame
    # terminal nodes can not be splitted, delta = 0 
    init.tree$frame$delta = 0
    init.tree$frame$delta[!is.leaf] = init.tree$splits[, 3]
    
    # Calculate g.h for each non-terminal node
    g.h <- numeric(nrow(init.tree$frame)) 
    child.list <- search.childs(is.leaf)
    
    for(i in 1:nrow(init.tree$frame)) {
      
      childs.i <- child.list[[i]]
      
      if(is.leaf[i]) {
        g.h[i] <- Inf
      } else{
        g.h[i] <- sum(init.tree$frame$delta[childs.i]) / length(!is.leaf(childs.i))
      }
    }
    
    # add g.h of each node to the frame
    init.tree$frame$g.h = g.h
  }

  

  # TODO: estimate treatment effect on each node with data.est  
  # add est.trt.eff and variance into init.tree$frame
  est.trt.eff <- numeric(nrow(init.tree$frame))
  variance    <- numeric(nrow(init.tree$frame))
  n.est    <- numeric(nrow(init.tree$frame))
  
  if (est.modelFitBefore) {
    
    est.cond.outcome  <- est.outcome(form.outcome, data.est, method = outcome.mthd, folds = 5)
    est.prop.score    <- est.pscore(form.pscore, data.est, method = pscore.mthd, folds = 5)
    
    data.est$cond.outcome <- est.cond.outcome
    data.est$pscore       <- est.prop.score
    
  }
  
  right.ind  <- 1
  next.r <- NULL
  
  for (h in 1:nrow(init.tree$frame)) {
    
    # print(paste("node h = ", h))
    
    if (h == 1) { # root
      
      data.node.est <- data.est
      
      # est.cond.outcome  <- est.outcome(form.outcome, data.node.est, method = outcome.mthd, folds = 5)
      # est.prop.score    <- est.pscore(form.pscore, data.node.est, method = pscore.mthd, folds = 5)
      # 
      # v.hat  <- data.node.est$A - est.prop.score
      # Y.star <- data.node.est$Y - est.cond.outcome
      # 
      # trt.eff.root <- sum(v.hat * Y.star) / sum(v.hat ^ 2)
      # 
      # var.root <- sum((Y.star - trt.eff.root * v.hat)^2 * v.hat^2) / (sum(v.hat ^ 2))^2
      #  
      # est.trt.eff[1] <- trt.eff.root
      # variance[1]    <- var.root
      # n.est[1] <- nrow(data.node.est)
      
    } else if (init.tree$frame$var[h] == "<leaf>") { # leaf  
      
      next
      
    } else { # internal nodes
      
      if(!is.null(next.l)) { 
        data.node.est <- next.l[[1]]
      } else {# if left branchs go to the end, starting to est right branches.
        # right index subtract 1, point to current right child node.
        data.node.est <- as.data.frame(next.r[[right.ind - 1]]) 
        
        # update next.r: remove the used element next.r[[right.ind - 1]]
        # if no more element except the used one: empty out next.r to store new next.r 
        if (length(next.r) >  1){  
          next.r[[right.ind - 1]] <- NULL
        } else {
          next.r <- NULL
        }
        
        # data in [[right.ind - 1]] has been used, update right index
        right.ind <- right.ind - 1 
      }
    }  
    

    # the row index of node h in splits, i.e. the times of splits up to node h
    h.ind.splits <- length(which(!is.leaf[1:h]))
    splt.var     <- init.tree$frame$var[h]
    splt.point   <- init.tree$splits[h.ind.splits, 4]
    direction    <- init.tree$splits[h.ind.splits, 2]
    
    splt.ind     <- which(colnames(data.node.est) == as.character(splt.var))
    
    if (class(data.est[, splt.ind]) == "factor") {
      
      level <- levels(data.est[, splt.ind])
      
      ind.est.l  <- which(data.node.est[, splt.ind] %in% level[init.tree$csplit[splt.point, ] == 1])
      ind.est.r  <- which(data.node.est[, splt.ind] %in% level[init.tree$csplit[splt.point, ] == 3])
      
      data.est.l <- data.node.est[ind.est.l, ]
      data.est.r <- data.node.est[ind.est.r, ]
      
    } else {
      
      if (direction > 0) {
        
        ind.est.l  <- which(data.node.est[, splt.ind] >= splt.point)
        ind.est.r  <- which(data.node.est[, splt.ind] <  splt.point)
        
        data.est.l <- data.node.est[ind.est.l, ]       # [data.node.est[, splt.ind] >= splt.point, ]
        data.est.r <- data.node.est[ind.est.r, ]       # [data.node.est[, splt.ind] < splt.point, ]
         
      } else {
        
        ind.est.l  <- which(data.node.est[, splt.ind] <  splt.point)
        ind.est.r  <- which(data.node.est[, splt.ind] >= splt.point)
        
        data.est.l <- data.node.est[ind.est.l, ]     # [data.node.est[, splt.ind] < splt.point, ]
        data.est.r <- data.node.est[ind.est.r, ]     # [data.node.est[, splt.ind] >= splt.point, ]
        
      }
      
    }
    
    
    if (init.tree$frame$var[h+1] == "<leaf>") {
      next.l <- NULL
    } else {
      next.l <- list(data.est.l)
    } 
    
    r.child.ind <- as.numeric(rownames(init.tree$frame)[h+1]) + 1
    r.child.row <- which(as.numeric(rownames(init.tree$frame)) == r.child.ind)
    
    if(init.tree$frame$var[r.child.row] != "<leaf>") {
      next.r[[right.ind]] <- data.est.r
      right.ind <- right.ind + 1
    }
    
    if (!est.modelFitBefore) {
      
      # treated units and control units in each node must larger than minsize
      if (sum(data.node.est$A == 1) < minsize | sum(data.node.est$A == 0) < minsize) {
        next
      } else {
        # remove the covariates columns if they have only one level
        data.node.est.update <- data.node.est[, c(sapply(data.node.est[, !colnames(data.node.est) %in% c("A", "Y")],
                                                             function(col) {length(unique(col))}) > 1, TRUE, TRUE)]
        
        # update formula
        form.outc.r       <- unlist(strsplit(form.outcome, " ~ "))[2]
        form.outc.r.terms <- unlist(strsplit(form.outc.r, " + ", fixed = TRUE))
        form.outc.update  <- paste("Y ~ ", paste0(form.outc.r.terms[form.outc.r.terms %in% colnames(data.node.est.update)], collapse = " + "), sep = "")
        
        form.psc.r        <- unlist(strsplit(form.pscore, " ~ "))[2]
        form.psc.r.terms  <- unlist(strsplit(form.psc.r, " + ", fixed = TRUE))
        form.psc.update   <- paste("A ~ ", paste0(form.psc.r.terms[form.psc.r.terms %in% colnames(data.node.est.update)], collapse = " + "), sep = "")
        
        
        if(nrow(data.node.est.update) >= 2*minsize*5) {
          
          est.prop.score   <- est.pscore(form.psc.update, data.node.est.update, method = pscore.mthd, folds = 5)
          est.cond.outcome <- est.outcome(form.outc.update, data.node.est.update, method = outcome.mthd, folds = 5, type.outcome = "continuous")
          
        } else {
          
          est.prop.score   <- est.pscore(form.psc.update, data.node.est.update, method = pscore.mthd, crossfit = FALSE)
          est.cond.outcome <- est.outcome(form.outc.update, data.node.est.update, method = outcome.mthd, crossfit = FALSE, type.outcome = "continuous")
          
        }
        
        v.hat  <- data.node.est.update$A - est.prop.score
        Y.star <- data.node.est.update$Y - est.cond.outcome 
        
      }
      
    } else {
      
      v.hat  <- data.node.est$A - data.node.est$pscore
      Y.star <- data.node.est$Y - data.node.est$cond.outcome
      
    } 
    
    ## TODO: estimate treatment effect
    if(h==1){ # estimate te and var of root node
      
      trt.eff.root <- sum(v.hat * Y.star) / sum(v.hat ^ 2)
      
      var.root <- sum((Y.star - trt.eff.root * v.hat)^2 * v.hat^2) / (sum(v.hat ^ 2))^2
      
      est.trt.eff[1] <- trt.eff.root
      variance[1]    <- var.root
      n.est[1] <- nrow(data.node.est)
      
    }
    
    trt.eff.l <- sum(v.hat[ind.est.l] * Y.star[ind.est.l]) / sum(v.hat[ind.est.l] ^ 2)
    trt.eff.r <- sum(v.hat[ind.est.r] * Y.star[ind.est.r]) / sum(v.hat[ind.est.r] ^ 2)
    
    est.trt.eff[h+1]         <- trt.eff.l
    est.trt.eff[r.child.row] <- trt.eff.r
    
    if (est.modelFitBefore) {
      n.tar <- nrow(data.est)
    } else {
      n.tar <- nrow(data.node.est)
    }
    
    # the variance of estimated treatment effect
    var.left  <-  sum((Y.star[ind.est.l] - trt.eff.l * v.hat[ind.est.l])^2
                                                    * v.hat[ind.est.l]^2) / (sum(v.hat[ind.est.l] ^ 2))^2
    var.right <-  sum((Y.star[ind.est.r] - trt.eff.r * v.hat[ind.est.r])^2
                                                    * v.hat[ind.est.r]^2) / (sum(v.hat[ind.est.r] ^ 2))^2
    
    variance[h+1]         <- var.left #/n.tar
    variance[r.child.row] <- var.right #/n.tar
 
    n.est[h+1]         <- length(ind.est.l)
    n.est[r.child.row] <- length(ind.est.r)
    
  }
  
  init.tree$frame$yval      <- est.trt.eff
  init.tree$frame$variance  <- variance
  init.tree$frame$n.est     <- n.est
  
  # confidence interval
  init.tree$frame$ci.lower <- est.trt.eff - sqrt(variance) * qnorm(0.975) # qnorm(1-(1-sqrt(0.9))/2) # 
  init.tree$frame$ci.upper <- est.trt.eff + sqrt(variance) * qnorm(0.975) # qnorm(1-(1-sqrt(0.9))/2) # 
  init.tree$frame$widthCI  <- 2 * sqrt(variance) * qnorm(0.975) # init.tree$frame$ci.upper - init.tree$frame$ci.lower
  

  
  
  #### prune
  ## TODO: find minimum lambda, i.e. minimum g.h of all nodes in current tree.
  subtree.seq <- vector("list")
  lambda <- vector("numeric")
  
  # initialize
  subtree.seq[[1]] <- init.tree 
  stop.prune <- FALSE 
  s <- 1
  
  while (stop.prune == FALSE ) {
    
    print(paste("Now starting to prune subtree s = ", s))
    
    tree.to.prune   <- subtree.seq[[s]]
    frame.current   <- tree.to.prune$frame
    is.leaf.current <- (frame.current$var == "<leaf>")
    splits.current  <- tree.to.prune$splits
    
    lambda[s] <- min(frame.current$g.h)
    prune.ind <- which.min(frame.current$g.h)
    
    # find the subtree that root in node h that with minimum g.h
    subtree.h <- search.node.childs(prune.ind, is.leaf.current)
    
    # update the frame 
    childs.h <- subtree.h[-1]
    frame.update <- frame.current[-childs.h, ]
    frame.update$var[prune.ind] <- "<leaf>"
    frame.update$delta[prune.ind] <- 0  
    
    # update the splits
    # the row index of node h ("prune.ind") in splits, i.e. the times of splits up to node h
    h.ind.splits <- length(which(!is.leaf.current[1:prune.ind]))
    
    # the times of splits of subtree.h
    nsplit.h <- length(which(!is.leaf.current[subtree.h]))
    if (prune.ind == 1) { # if prune the root node, new tree have no splits
      splits.update <- NULL
    } else {
      splits.update <- splits.current[-c(h.ind.splits : (h.ind.splits + nsplit.h - 1)), ,drop = FALSE]
    }
    
    
    tree.new <- tree.to.prune
    tree.new$frame <- frame.update
    
    tree.new$splits <- splits.update   

    # when splits have only one line, class(splits.current) = "numeric"
    # if (class(splits.update) == "numeric") { # if prune the root node
    #   tree.new$splits <- as.matrix(splits.update)
    # } else if (class(splits.update) == "matrix") {
    #   tree.new$splits <- splits.update
    # }
    
    # update where
    tree.to.prune$where[tree.to.prune$where %in% subtree.h] <- prune.ind
    tree.to.prune$where[tree.to.prune$where > max(subtree.h)] <- tree.to.prune$where[tree.to.prune$where > max(subtree.h)] - length(childs.h)
    tree.new$where <- tree.to.prune$where
    
    # update num.leaf, leaf.size, leaf.obs
    tree.new$num.leaf <- sum(tree.new$frame$var == "<leaf>")
    tree.new$leaf.size <- tree.new$frame$n[which(tree.new$frame$var != "<leaf>")]
    # leaf.obs  <- which(tree$where == sort(unique(tree.new$where)))
    
    leaf.obs <- list(tree.new$num.leaf)
    for (g in 1:tree.new$num.leaf) {
      tree.new$leaf.obs[[g]] <- names(tree.new$where)[tree.new$where == which(tree.new$frame$var == "<leaf>")[g]]
    }
    # names(leaf.obs) <- paste0("Subgroup", 1:tree.new$num.leaf, sep=" ")
    
    # update g(h) for the new tree
    is.leaf.update     <- (tree.new$frame$var == "<leaf>")
    tree.new$frame$delta[is.leaf.update] <- 0
    tree.new$frame$g.h <- numeric(nrow(tree.new$frame)) 
    child.list.update  <- search.childs(is.leaf.update)
    
    for(i in 1:nrow(tree.new$frame)) {
      if(is.leaf.update[i]) {
        tree.new$frame$g.h[i] <- Inf
      } else{
        subtree.i <- child.list.update[[i]]
        tree.new$frame$g.h[i] <- sum(tree.new$frame$delta[subtree.i]) / length(!is.leaf.update[subtree.i])
      }
    }
    
    s <- s + 1
    subtree.seq[[s]] <- tree.new
    if (sum(tree.new$frame$var == "<leaf>") == 1) {stop.prune = TRUE}
  }
  
  
  
  
  #### Validation: Optimal subtree selection
  #### avaliable: subtree.seq, lambda 
  
  ## TODO: calculate G.T for each subtree
  # if model fit before tree building, estimates need to be input into each node
  if (modelFitBefore) {
    
    valid.cond.outcome  <- est.outcome(form.outcome, data.valid, method = outcome.mthd, folds = 5)
    valid.pscore        <- est.pscore(form.pscore, data.valid, method = pscore.mthd, folds = 5)
    
    # if modelFitBefore, add columns pscore and cond.outcome into data.valid
    data.valid$cond.outcome <- valid.cond.outcome
    data.valid$pscore       <- valid.pscore
    
  }
  
  # initialize g(h), lambda and G(tree)
  G.t <- numeric(length(subtree.seq))
  nsplit <- numeric(length(subtree.seq))
  
  for (t in 1:length(subtree.seq)) { # t in 1:2
    
    print(paste("Now starting to valid subtree t = ", t))
    
    tree.to.valid <- subtree.seq[[t]]
    
    if (nrow(tree.to.valid$frame) == 1) {
      
      is.leaf.valid <- (tree.to.valid$frame$var == "<leaf>")
      valid.goodness <- 0
      G.t[t] <- 0
      
    } else {
      
      is.leaf.valid <- (tree.to.valid$frame$var == "<leaf>")
      valid.goodness <- 0
      
      right.ind  <- 1
      next.r <- NULL
      
      # put the data.valid into subtree
      for (h in 1:nrow(tree.to.valid$frame)) {
        
        # print(paste("node h = ", h))
        # 3 types: root, leaf, internal nodes
        if (h == 1) { # root
          data.node.valid <- data.valid
        } else if (tree.to.valid$frame$var[h] == "<leaf>") { # leaf
          next
        } else { # internal nodes
          if(!is.null(next.l)) {
            data.node.valid <- next.l[[1]]
          } else {   # if left branchs go to the end, starting to valid right branches.
            
            # right index subtract 1, poinit to current right child node.
            data.node.valid <- as.data.frame(next.r[[right.ind - 1]]) 
            
            # update next.r: remove the used element next.r[[right.ind - 1]]
            # if no more element except the used one: empty out next.r to store new next.r 
            if (length(next.r) >  1){  
              next.r[[right.ind - 1]] <- NULL
            } else {
              next.r <- NULL
            }
            
            # data in [[right.ind - 1]] has been used, update right index
            right.ind <- right.ind - 1 
          }
        }   
        
        
        # the row index of node h in splits, i.e. the times of splits up to node h
        h.ind.splits <- length(which(!is.leaf.valid[1:h]))
        splt.var     <- tree.to.valid$frame$var[h]
        splt.point   <- tree.to.valid$splits[h.ind.splits, 4]
        direction    <- tree.to.valid$splits[h.ind.splits, 2]
        
        splt.ind     <- which(colnames(data.node.valid) == splt.var)
        
        if (class(data[, as.character(splt.var)]) == "factor") {
          
          level <- levels(data.valid[, splt.ind])
          
          ind.valid.l  <- which(data.node.valid[, splt.ind] %in% level[tree.to.valid$csplit[splt.point, ] == 1])
          ind.valid.r  <- which(data.node.valid[, splt.ind] %in% level[tree.to.valid$csplit[splt.point, ] == 3])
          
          data.valid.l <- data.node.valid[ind.valid.l, ]
          data.valid.r <- data.node.valid[ind.valid.r, ]
          
        } else{
          
          if (direction > 0) {
            
            ind.valid.l  <- which(data.node.valid[, splt.ind] >= splt.point)
            ind.valid.r  <- which(data.node.valid[, splt.ind] <  splt.point)
            
            data.valid.l <- data.node.valid[ind.valid.l, ]       # [data.node.valid[, splt.ind] >= splt.point, ]
            data.valid.r <- data.node.valid[ind.valid.r, ]       # [data.node.valid[, splt.ind] < splt.point, ]
            
          } else {
            
            ind.valid.l  <- which(data.node.valid[, splt.ind] <  splt.point)
            ind.valid.r  <- which(data.node.valid[, splt.ind] >= splt.point)
            
            data.valid.l <- data.node.valid[ind.valid.l, ]     # [data.node.valid[, splt.ind] < splt.point, ]
            data.valid.r <- data.node.valid[ind.valid.r, ]     # [data.node.valid[, splt.ind] >= splt.point, ]
            
          }
          
        }
        
        
        if (tree.to.valid$frame$var[h+1] == "<leaf>") {
          next.l <- NULL
        } else {
          next.l <- list(data.valid.l)
        } 
        
        r.child.ind <- as.numeric(rownames(tree.to.valid$frame)[h+1]) + 1
        r.child.row <- which(as.numeric(rownames(tree.to.valid$frame)) == r.child.ind)
        
        if(tree.to.valid$frame$var[r.child.row] != "<leaf>") {
          next.r[[right.ind]] <- data.valid.r
          right.ind <- right.ind + 1
        }
        
        
        # treated units and control units in each node must larger than minsize
        if (!modelFitBefore) {
          
          if (sum(data.node.valid$A == 1) < minsize | sum(data.node.valid$A == 0) < minsize) {
            
            next
            
          } else {
            
            # remove the covariates columns if they have only one level
            data.node.valid.update <- data.node.valid[, c(sapply(data.node.valid[, !colnames(data.node.valid) %in% c("A", "Y")],
                                                                 function(col) {length(unique(col))}) > 1, TRUE, TRUE)]
            # data.node.valid.update <- lapply(data.node.valid.update,
            #                                  function(x) if(is.factor(x)) factor(x) else x)
            
            # update formula
            form.outc.r       <- unlist(strsplit(form.outcome, " ~ "))[2]
            form.outc.r.terms <- unlist(strsplit(form.outc.r, " + ", fixed = TRUE))
            form.outc.update  <- paste("Y ~ ", paste0(form.outc.r.terms[form.outc.r.terms %in% colnames(data.node.valid.update)], collapse = " + "), sep = "")
            
            form.psc.r        <- unlist(strsplit(form.pscore, " ~ "))[2]
            form.psc.r.terms  <- unlist(strsplit(form.psc.r, " + ", fixed = TRUE))
            form.psc.update   <- paste("A ~ ", paste0(form.psc.r.terms[form.psc.r.terms %in% colnames(data.node.valid.update)], collapse = " + "), sep = "")
            
            valid.pscore       <- est.pscore(form.psc.update, data.node.valid.update, method = pscore.mthd, crossfit = FALSE, folds = 5)
            
            valid.cond.outcome <- est.outcome(form.outc.update, data.node.valid.update, method = outcome.mthd,
                                              crossfit = FALSE, folds = 5, type.outcome = "continuous")
            
            v.hat  <- data.node.valid.update$A - valid.pscore
            Y.star <- data.node.valid.update$Y - valid.cond.outcome 
            
          }
          
        } else {
          
          if (sum(data.node.valid$A == 1) < minsize | sum(data.node.valid$A == 0) < minsize) {
            
            next
            
          } else {
            
            v.hat  <- data.node.valid$A - data.node.valid$pscore
            Y.star <- data.node.valid$Y - data.node.valid$cond.outcome 
            
          }
        }
        
        
        ## TODO: estimate treatment effect
        trt.eff.l <- sum(v.hat[ind.valid.l] * Y.star[ind.valid.l]) / sum(v.hat[ind.valid.l] ^ 2)
        trt.eff.r <- sum(v.hat[ind.valid.r] * Y.star[ind.valid.r]) / sum(v.hat[ind.valid.r] ^ 2)
        
        delta.trt.eff <- trt.eff.l - trt.eff.r
        
        # calculate the variance estimator
        if (modelFitBefore) {
          n.tar <- nrow(data.valid)
        } else {
          n.tar <- nrow(data.node.valid)
        }
        
        var.left  <- sum((Y.star[ind.valid.l] - trt.eff.l * v.hat[ind.valid.l])^2
                                                 * v.hat[ind.valid.l]^2) / (sum(v.hat[ind.valid.l]^2))^2
        var.right <- sum((Y.star[ind.valid.r] - trt.eff.r * v.hat[ind.valid.r])^2
                                                 * v.hat[ind.valid.r]^2) / (sum(v.hat[ind.valid.r]^2))^2

        valid.goodness <- valid.goodness + (delta.trt.eff) ^ 2  / (var.left + var.right)
         
      } # end h
      
    } # end else (nrow(tree.to.valid$frame) > 1)
    # valid.lambda[t] <- valid.goodness / sum(!is.leaf.valid)
    G.t[t] <- valid.goodness - cp * sum(!is.leaf.valid)  # lambda[t]
    
    # generate cptable, contains complexity parameter lambda, number of split nsplit, validation complexity of heterogeneity.
    nsplit[t] <- sum(!is.leaf.valid)  # nrow(subtree.seq[[t]]$splits) 
  } # end t
  
  res.tab <- data.frame(subtree = 1:length(subtree.seq),
                        nsplit  = nsplit, 
                        val.complx = G.t) 
  
  opt.ind <- which.max(G.t)
  opt.tree <- subtree.seq[[opt.ind]]
  
  # class(opt.tree) <- "causalPartitionTree"
  return(list(optree = opt.tree, 
              res.tab = res.tab, 
              init.tree = init.tree,
              subtree.list = subtree.seq,
              lambda.list = lambda))
  
}



