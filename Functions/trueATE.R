TrueATE <- function(tree, data.test, true.trt.eff) {
  
  
  if (nrow(tree$frame) == 1) {
    
   true.avg.trt.eff <- mean(true.trt.eff, na.rm = TRUE)
   n.test           <- nrow(data.test)
    
  } else{
    
    #initialize
    true.avg.trt.eff <- numeric(nrow(tree$frame))
    n.test           <- numeric(nrow(tree$frame))
    
    right.ind  <- 1
    next.r <- NULL
    
    is.leaf <- (tree$frame$var == "<leaf>")
    data.test <- cbind(data.test, true.trt.eff)
    
    for (h in 1:nrow(tree$frame)) {
      
      # print(paste("node h = ", h))
      
      if (h == 1) { # root
        
        data.node.test <- data.test
        
      } else if (tree$frame$var[h] == "<leaf>") { # leaf  
        
        next
        
      } else { # internal nodes
        
        if(!is.null(next.l)) { 
          data.node.test <- next.l[[1]]
        } else {# if left branchs go to the end, starting to est right branches.
          # right index subtract 1, point to current right child node.
          data.node.test <- as.data.frame(next.r[[right.ind - 1]]) 
          
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
      splt.var     <- tree$frame$var[h]
      splt.point   <- tree$splits[h.ind.splits, 4]
      direction    <- tree$splits[h.ind.splits, 2]
      
      splt.ind     <- which(colnames(data.node.test) == as.character(splt.var))
      
      if (class(data.test[, splt.ind]) == "factor") {
        
        level <- levels(data.test[, splt.ind])
        
        ind.test.l  <- which(data.node.test[, splt.ind] %in% level[tree$csplit[splt.point, ] == 1])
        ind.test.r  <- which(data.node.test[, splt.ind] %in% level[tree$csplit[splt.point, ] == 3])
        
        data.test.l <- data.node.test[ind.test.l, ]
        data.test.r <- data.node.test[ind.test.r, ]
        
      } else {
        
        if (direction > 0) {
          
          ind.test.l  <- which(data.node.test[, splt.ind] >= splt.point)
          ind.test.r  <- which(data.node.test[, splt.ind] <  splt.point)
          
          data.test.l <- data.node.test[ind.test.l, ]       # [data.node.est[, splt.ind] >= splt.point, ]
          data.test.r <- data.node.test[ind.test.r, ]       # [data.node.est[, splt.ind] < splt.point, ]
          
        } else {
          
          ind.test.l  <- which(data.node.test[, splt.ind] <  splt.point)
          ind.test.r  <- which(data.node.test[, splt.ind] >= splt.point)
          
          data.test.l <- data.node.test[ind.test.l, ]     # [data.node.est[, splt.ind] < splt.point, ]
          data.test.r <- data.node.test[ind.test.r, ]     # [data.node.est[, splt.ind] >= splt.point, ]
          
        }
        
      }
      
      
      if (tree$frame$var[h+1] == "<leaf>") {
        next.l <- NULL
      } else {
        next.l <- list(data.test.l)
      } 
      
      r.child.ind <- as.numeric(rownames(tree$frame)[h+1]) + 1
      r.child.row <- which(as.numeric(rownames(tree$frame)) == r.child.ind)
      
      if(tree$frame$var[r.child.row] != "<leaf>") {
        next.r[[right.ind]] <- data.test.r
        right.ind <- right.ind + 1
      }
      
      
      
      ## TODO: calculate the true treatment effect
      
      if(h==1){
        
        true.avg.trt.eff[1] <- mean(data.node.test$true.trt.eff, na.rm = TRUE)
        n.test[1] <- nrow(data.node.test)
        
      } 
      
      trt.eff.l <- mean(data.test.l$true.trt.eff, na.rm = TRUE)
      trt.eff.r <- mean(data.test.r$true.trt.eff, na.rm = TRUE)
      
      true.avg.trt.eff[h+1]         <- trt.eff.l
      true.avg.trt.eff[r.child.row] <- trt.eff.r
      
      n.test[h+1]         <- length(ind.test.l)
      n.test[r.child.row] <- length(ind.test.r)
      
    }
    
  }
  
  
  tree$frame$true.avg.trt.eff <- true.avg.trt.eff
  tree$frame$n.test           <- n.test
  
  return(tree$frame)
  
}
