## search for child nodes of nodes in init.tree$frame.
# If node is a terminal node, return NULL
# If node is a internal node, return the set of child nodes.

rpart.kids <- function(i, is.leaf) {
  if (is.leaf[i]) {
    return(NULL)
  } else {
    return(c(i + 1L, 
                which((cumsum(!is.leaf[-(1L:i)]) + 1L) == cumsum(is.leaf[-(1L:i)]))[1L] + 1L + i))
  }
}




search.childs <- function(is.leaf){
  
  childs <- vector("list", length(is.leaf))
  
  for (node in 1:length(is.leaf)) {
    
    if (is.leaf[node]) {
      childs[[node]] <- node
    } else {
      chds.seq     <- node
      chds.seq.old <- NULL
      while (length(chds.seq.old) < length(chds.seq)) {
        chds.seq.old <- chds.seq
        for (i in chds.seq) {
          chds.seq <- unique(append(chds.seq, rpart.kids(i, is.leaf)))
        }
      }
      childs[[node]] <- chds.seq
    }
  }
  
  childs
}





search.node.childs <- function(node, is.leaf){
  
  if (is.leaf[node]) {
    childs <- node
  } else {
    chds.seq     <- node
    chds.seq.old <- NULL
    while (length(chds.seq.old) < length(chds.seq)) {
      chds.seq.old <- chds.seq
      for (i in chds.seq) {
        chds.seq <- unique(append(chds.seq, rpart.kids(i, is.leaf)))
      }
    }
    childs <- chds.seq
  }
  
  childs
}

