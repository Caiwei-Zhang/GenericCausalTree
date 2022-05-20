est.outcome <- function(form.outcome, data.node, method = c("RF", "GBM", "GLM"), 
                        crossfit, folds, type.outcome = "continuous") {
  
  # estimate conditional mean of outcome l(x)=E(Y|X)
  data.xy <- data.node[, !colnames(data.node) %in% c("A")]
  form.r       <- unlist(strsplit(form.outcome, " ~ "))[2]
  form.r.terms <- unlist(strsplit(form.r, " + ", fixed = TRUE))
  
  Y.hat  <- numeric(nrow(data.node))

  if (missing(crossfit)) {
    if(missing(folds)) { # "crossfit" & "folds" missing, without cross-fitting.
      crossfit <- FALSE
    } else {  # "crossfit" missing, "folds" is given, cross-fitting.
      crossfit <- TRUE
    }
  }
  
  
  if (crossfit == FALSE) {
    
    X.pred <- data.xy[, !colnames(data.xy) %in% c("Y")] #0414?? !colnames(data.xy) %in% c("Y")
      
    if (method == "RF"){
      form.outcome <- as.formula(form.outcome)
      modell <- rfsrc(form.outcome, data = data.xy, ntree = 500)
      Y.hat  <- predict(modell, X.pred)$predicted
    }
    
    if (method == "GBM") {
      
      # gbm
      modell <- gbm.fit(x = data.xy[, colnames(data.xy) %in% form.r.terms],
                        y = data.xy$Y, distribution = "gaussian",
                        n.trees = 100, shrinkage = 0.1,
                        interaction.depth = 1, bag.fraction = 1,
                        n.minobsinnode = 5, 
                        nTrain = nrow(data.xy), verbose = FALSE, keep.data = FALSE)
      
      Y.hat  <- predict(modell, X.pred, type = "response")

    }
    
    if (method == "XGB") {
      # xgb
      ddata <- xgb.DMatrix(as.matrix(data.xy[, !colnames(data.xy) %in% c("Y")]),
                           label = data.xy$Y)
      para <- list(booster = "gbtree", max_depth = 1)
      modell <- xgb.train(data = ddata, params = para, nrounds = 500, nthreads = 4)

      dtest <- xgb.DMatrix(as.matrix(X.pred))
      Y.hat <- predict(modell, dtest)
      
    }
    
    if (method == "GLM") {
      
      form.outcome <- as.formula(form.outcome)
      if (type.outcome == "continuous") {
        modell <- glm(formula = form.outcome, data = data.xy, family = gaussian(link = "identity"))
      } else if (type.outcome == "binary") {
        modell <- glm(formula = form.outcome, data = data.xy, family = binomial(link = "logit"))
      }
      
      Y.hat  <- predict(modell, X.pred, type = "response")
      
    }
    
  } else {
    
    if (missing(folds)) {
      folds <- 5
    }
    
    idx <- createFolds(data.node$A, k = folds, list = TRUE, returnTrain = FALSE)
    
    for (fold in 1:folds) {
    
      main.idx <- idx[[fold]]
      aux.xy   <- data.xy[-main.idx, ]
      X.main   <- data.xy[main.idx, !colnames(data.xy) %in% c("Y")]
      
      if (method == "RF"){
        form.outcome <- as.formula(form.outcome)
        modell <- rfsrc(form.outcome, data = aux.xy, ntree = 500)
        Y.hat[main.idx]  <- predict(modell, X.main)$predicted
      }
      
      if (method == "GBM") {
        
        # gbm
        modell <- gbm.fit(x = aux.xy[, colnames(aux.xy) %in% form.r.terms],
                          y = aux.xy$Y, distribution = "gaussian",
                          n.trees = 100, shrinkage = 0.1,
                          interaction.depth = 1, bag.fraction = 1,
                          n.minobsinnode = 5, 
                          nTrain = nrow(aux.xy), verbose = FALSE, keep.data = FALSE)
        # n.minobsinnode = 1, 
        Y.hat[main.idx]  <- predict(modell, X.main, type = "response")
        
      }
      
      if (method == "XGB") {
        # xgb
        ddata <- xgb.DMatrix(as.matrix(aux.xy[, !colnames(aux.xy) %in% c("Y")]),
                             label = aux.xy$Y)
        para <- list(booster = "gbtree", max_depth = 1)
        modell <- xgb.train(data = ddata, params = para, nrounds = 500, nthreads = 4)

        dtest <- xgb.DMatrix(as.matrix(X.main))
        Y.hat[main.idx]  <- predict(modell, dtest)
        
      }
      
      if (method == "GLM") {
        form.outcome <- as.formula(form.outcome)
        if (type.outcome == "continuous") {
          modell <- glm(formula = form.outcome, data = aux.xy, family = gaussian(link = "identity"))
        } else if (type.outcome == "binary") {
          modell <- glm(formula = form.outcome, data = aux.xy, family = binomial(link = "logit"))
        }
        
        Y.hat[main.idx]  <- predict(modell, X.main, type = "response")
        
      }
    
    }
    
  }
  
  return(cond.outcome = Y.hat)
  
}




est.pscore <- function(form.pscore, data.node, method = c("RF", "GBM", "GLM"), crossfit, folds){
  
  # estimate the propensity scores e(x)
  data.xa <- data.node[, !colnames(data.node) %in% c("Y")]
  form.r       <- unlist(strsplit(form.pscore, " ~ "))[2]
  form.r.terms <- unlist(strsplit(form.r, " + ", fixed = TRUE))
  
  pscore <- numeric(nrow(data.node))

  if (missing(crossfit)) {
    if(missing(folds)) { # "crossfit" & "folds" missing, without cross-fitting.
      crossfit <- FALSE
    } else {  # "crossfit" missing, "folds" is given, cross-fitting.
      crossfit <- TRUE
    }
  }
  
  if (crossfit == FALSE) {
    
    X <- data.xa[, !colnames(data.xa) %in% c("A")]
    
    if (method == "RF") {
      
      form.pscore <- as.formula(form.pscore)
      modele <- rfsrc(form.pscore, data = data.xa, ntree = 500)
      pscore <- predict(modele, X)$predicted
      
    }
    
    if (method == "GBM") {

      modele <- gbm.fit(x = data.xa[, colnames(data.xa) %in% form.r.terms], 
                        y = data.xa$A, distribution = "bernoulli",
                        n.trees = 100, shrinkage = 0.1,
                        interaction.depth = 1, bag.fraction = 1,
                        n.minobsinnode = 5, 
                        nTrain = nrow(data.xa), verbose = FALSE, keep.data = FALSE)
      
      pscore <- predict(modele, X, type = "response")
      
    }
    
    if (method == "GLM") {
      
      form.pscore <- as.formula(form.pscore)
      modele <- glm(formula = form.pscore, data = data.xa, family = binomial(link = "logit"))
      pscore <- predict(modele, X, type = "response")
      
    }
    
  } else {
    
    if (missing(folds)) {
      folds <- 5
    }
    
    idx <- createFolds(data.node$A, k = folds, list = TRUE, returnTrain = FALSE)
    
    for (fold in 1:folds) {
      
      main.idx <- idx[[fold]]
      aux.xa   <- data.xa[-main.idx, ]
      X.main   <- data.xa[main.idx, colnames(data.xa) != c("A")]
      
      if (method == "RF"){
        form.pscore <- as.formula(form.pscore)
        modele <- rfsrc(form.pscore, data = aux.xa, ntree = 500)
        pscore[main.idx] <- predict(modele, X.main)$predicted
      }
      
      if (method == "GBM") {
        
        modele <- gbm.fit(x = aux.xa[, colnames(aux.xa) %in% form.r.terms], 
                          y = aux.xa$A, distribution = "bernoulli",
                          n.trees = 100, shrinkage = 0.1,
                          interaction.depth = 1, bag.fraction = 1,
                          n.minobsinnode = 5, 
                          nTrain = nrow(aux.xa), verbose = FALSE, keep.data = FALSE)
        
        pscore[main.idx] <- predict(modele, X.main, type = "response")
        
      }
      
      if (method == "GLM") {
        form.pscore <- as.formula(form.pscore)
        modele <- glm(formula = form.pscore, data = aux.xa, family = binomial(link = "logit"))
        pscore[main.idx] <- predict(modele, X.main, type = "response")
      }
      
    }
    
  }
  
  
  return(pscore = pscore)
  
}
