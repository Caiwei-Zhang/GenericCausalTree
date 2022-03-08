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
    
    X.pred <- data.xy[, !colnames(data.xy) %in% c("A")]
      
    if (method == "RF"){
      form.outcome <- as.formula(form.outcome)
      modell <- rfsrc(form.outcome, data = data.xy, ntree = 500)
      Y.hat  <- predict(modell, X.pred)$predicted
    }
    
    if (method == "GBM") {
      
      # form.outcome <- as.formula(form.outcome)
      # modell <- gbm(form.outcome, data = data.xy, distribution = "gaussian",
      #               n.trees = 300, interaction.depth = 2,  # bag.fraction = 0.5, train.fraction = 0.8
      #               n.minobsinnode = 5, shrinkage = 0.05)
      modell <- gbm.fit(x = data.xy[, colnames(data.xy) %in% form.r.terms],
                        y = data.xy$Y, distribution = "gaussian",
                        n.trees = 300, interaction.depth = 1,  # bag.fraction = 0.5,
                        n.minobsinnode = 5, shrinkage = 0.05,
                        nTrain = nrow(data.xy), verbose = FALSE, keep.data = FALSE)
      Y.hat  <- predict(modell, X.pred, type = "response")


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
        
        modell <- gbm.fit(x = aux.xy[, colnames(aux.xy) %in% form.r.terms],
                          y = aux.xy$Y, distribution = "gaussian",
                          n.trees = 300, interaction.depth = 1,  # bag.fraction = 0.5,
                          n.minobsinnode = 5, shrinkage = 0.05,
                          nTrain = nrow(aux.xy), verbose = FALSE, keep.data = FALSE)
        # gbmGrid <- expand.grid(n.trees = 300, interaction.depth = 2, n.minobsinnode = 5, shrinkage = 0.05)
        # modell <- train(form.outcome, data = aux.xy, method = "gbm",
        #                 preProcess = c('center', 'scale'), tuneGrid = gbmGrid)
        # Y.hat[main.idx]  <- predict(modell, X.main)
        # form.outcome <- as.formula(form.outcome)
        # modell <- gbm(form.outcome, data = data.xy, distribution = "gaussian",
        #               n.trees = 300, interaction.depth = 2,  # bag.fraction = 0.5, train.fraction = 0.8
        #               n.minobsinnode = 5, shrinkage = 0.05)#, cv.folds = 10

        Y.hat[main.idx]  <- predict(modell, X.main, type = "response")
        
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
                        n.trees = 200, interaction.depth = 1, # bag.fraction = 0.5,
                        n.minobsinnode = 5, shrinkage = 0.05,
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
                          n.trees = 300, interaction.depth = 1, # bag.fraction = 0.5,
                          n.minobsinnode = 5, shrinkage = 0.05,
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





# model.fit <- function(form.outcome, form.pscore, data.node, method = c("RF", "GBM", "GLM"), 
#                       folds = 2, type.outcome = "continuous"){
#   
#   # estimate conditional mean of outcome l(x)=E(Y|X) and propensity score e(x)
#   # data.node:    a data frame contains covariates X, outcome Y and treatment A.
#   # method:       one of "RF", "GBM", or "GLM". Default value is "GBM".
#   # folds:        number of data splitting
#   # type.outcome: type of outcome variable, if continuous, give predicted value; 
#   #               if binary, give the probability of Y=1. 
#   
#   method <- match.arg(method)
#   
#   idx <- createFolds(data.node$A, k = folds, list = TRUE, returnTrain = FALSE)
#   
#   data.xy <- data.node[, !colnames(data.node) %in% c("A")]
#   data.xa <- data.node[, !colnames(data.node) %in% c("Y")]
#   
#   Y.hat  <- numeric(nrow(data.node)) 
#   pscore <- numeric(nrow(data.node))
#   
#   for (fold in 1:folds) {
#     
#     main.idx <- idx[[fold]]
#     aux.xy   <- data.xy[-main.idx, ]
#     aux.xa   <- data.xa[-main.idx, ]
#     X.main   <- data.xy[main.idx, colnames(data.xy) != c("Y")]
#     
#     if (method == "RF"){
#       
#       modell <- rfsrc(as.formula(form.outcome), data = aux.xy, ntree = 500)
#       modele <- rfsrc(as.formula(form.pscore), data = aux.xa, ntree = 500)
#       
#       Y.hat[main.idx]  <- predict(modell, X.main)$predicted
#       pscore[main.idx] <- predict(modele, X.main)$predicted
#     }
#     
#     if (method == "GBM") {
#       
#       form.outc.r       <- unlist(strsplit(form.outcome, " ~ "))[2]
#       form.outc.r.terms <- unlist(strsplit(form.outc.r, " + ", fixed = TRUE))
#       
#       form.psc.r       <- unlist(strsplit(form.pscore, " ~ "))[2]
#       form.psc.r.terms <- unlist(strsplit(form.psc.r, " + ", fixed = TRUE))
#       
#       modell <- gbm.fit(x = aux.xy[, colnames(aux.xy) %in% form.outc.r.terms], 
#                         y = aux.xy$Y, distribution = "gaussian",
#                         n.trees = 100, interaction.depth = 2, bag.fraction = 0.8,
#                         n.minobsinnode = 5, shrinkage = 0.05, 
#                         nTrain = nrow(aux.xy))
#       
#       modele <- gbm.fit(x = aux.xa[, colnames(aux.xa) %in% form.psc.r.terms], 
#                         y = aux.xa$A, distribution = "bernoulli",
#                         n.trees = 100, interaction.depth = 1, bag.fraction = 0.8,
#                         n.minobsinnode = 5, shrinkage = 0.05,
#                         nTrain = nrow(aux.xa))
#       
#       Y.hat[main.idx]  <- predict(modell, X.main, type = "response")
#       pscore[main.idx] <- predict(modele, X.main, type = "response")
#       
#     }
#     
#     if (method == "GLM") {
#       
#       if (type.outcome == "continuous") {
#         modell <- glm(as.formula(form.outcome), data = aux.xy, family = gaussian(link = "identity"))
#       } else if (type.outcome == "binary") {
#         modell <- glm(as.formula(form.outcome), data = aux.xy, family = binomial(link = "logit"))
#       }
#       
#       modele <- glm(as.formula(form.pscore), data = aux.xa, family = binomial(link = "logit"))
#       
#       Y.hat[main.idx]  <- predict(modell, X.main, type = "response")
#       pscore[main.idx] <- predict(modele, X.main, type = "response")
#       
#     }
#     
#   }
#   
#   # truncate extreme value
#   pscore[pscore < 0.1] <- 0.1
#   pscore[pscore > 0.9] <- 0.9
#   
#   return(list(cond.outcome = Y.hat, 
#               pscore       = pscore))
#   
# }


