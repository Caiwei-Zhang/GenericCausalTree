CTGrid <- expand.grid(split_rule = c("TOT", "CT", "fit", "tstats"), split_honest = c(T, F),
                      cv_option  = c("TOT", "CT", "fit", "matching"), cv_honest   = c(T, F))

CTGrid <- CTGrid %>%
  mutate(split_honest = ifelse(split_rule %in% c("TOT"), NA, split_honest)) %>%
  mutate(cv_honest    = ifelse(cv_option %in% c("TOT", "matching"), NA, cv_honest))
CTGrid <- CTGrid[!duplicated(CTGrid),] # quchong
rownames(CTGrid) <- 1:nrow(CTGrid)


set.seed(999)

pscore.coeff <- matrix(c(0.5, -0.5, 0.5, -0.5, 0.5), 5, 1)
################################################################################################
######################### 1. tuning CT: heterogeneous, nonhonest ###############################
################################################################################################

tune.ct.hetero.nonhonest.simu.func <- function(ind) {
  
  data.lst <- makeData.cont.eff.cont(N = 1000, n.test = 1000, p = 5)
  data <- data.lst$data.used
  # data.noy <- data %>% dplyr::select(-Y)
  
  form.term <- colnames(data)[!colnames(data) %in% c("A", "Y")]
  form <- paste("A ~", paste0(form.term, collapse = " + "))
  
  t0 <- Sys.time()
  tmp.propsc <- est.pscore(form.pscore = form, data.node = data, method = "GLM", crossfit = FALSE)
  pw <- ifelse(data$A == 1, tmp.propsc, 1 - tmp.propsc)
  
  if (is.na(CTGrid$split_honest[ind])) {
    if (is.na(CTGrid$cv_honest[ind])) {
      tmp.ct.nonhonest <- causalTree::causalTree(Y ~ ., data = data, weights = 1/pw,
                                                 treatment   = data$A, 
                                                 split.Rule  = CTGrid$split_rule[ind], 
                                                 cv.option   = CTGrid$cv_option[ind], 
                                                 split.Bucket = T, xval = 10,
                                                 cp = 0, minsize = 20) 
    } else { # !is.na(cv_honest)
      tmp.ct.nonhonest <- causalTree::causalTree(Y ~ ., data = data, weights = 1/pw,
                                                 treatment   = data$A, 
                                                 split.Rule  = CTGrid$split_rule[ind], 
                                                 cv.option   = CTGrid$cv_option[ind], 
                                                 cv.Honest   = CTGrid$cv_honest[ind], 
                                                 split.Bucket = T, xval = 10,
                                                 cp = 0, minsize = 20) 
    }
    
  } else { # !is.na(split_honest)
    
    if (is.na(CTGrid$cv_honest[ind])) {
      tmp.ct.nonhonest <- causalTree::causalTree(Y ~ ., data = data, weights = 1/pw,
                                                 treatment   = data$A, 
                                                 split.Rule  = CTGrid$split_rule[ind], 
                                                 cv.option   = CTGrid$cv_option[ind], 
                                                 split.Honest = CTGrid$split_honest[ind],
                                                 split.Bucket = T, xval = 10,
                                                 cp = 0, minsize = 20) 
    } else { # !is.na(cv_honest)
      tmp.ct.nonhonest <- causalTree::causalTree(Y ~ ., data = data, weights = 1/pw,
                                                 treatment   = data$A, 
                                                 split.Rule  = CTGrid$split_rule[ind], 
                                                 cv.option   = CTGrid$cv_option[ind],
                                                 split.Honest = CTGrid$split_honest[ind],
                                                 cv.Honest   = CTGrid$cv_honest[ind], 
                                                 split.Bucket = T, xval = 10,
                                                 cp = 0, minsize = 20) 
    }
    
  }
  
  opcp  <- tmp.ct.nonhonest$cptable[,1][which.min(tmp.ct.nonhonest$cptable[,4])]
  final.ct.nonhonest <- prune(tmp.ct.nonhonest, cp = opcp)
  
  t1 <- Sys.time()
  
  
  X.test <- as.matrix(data.lst$test.data)
  e.test <- exp(X.test[, 1:5] %*% pscore.coeff)/(1 + exp(X.test[, 1:5] %*% pscore.coeff))
  
  tmp.eval.ct <- eval.measures.eff(final.tree   = final.ct.nonhonest,
                                   test.data    = data.frame(data.lst$test.data, A = rbinom(1000, 1, e.test)),
                                   true.trt.eff = data.lst$true.trt.eff,
                                   noise.var    = data.lst$noise.splt.var,
                                   corr.split   = data.lst$eff.modifier,
                                   type.eff.mod = data.lst$type.eff.mod,
                                   where.split  = data.lst$where.split,
                                   dir.split    = data.lst$dir.split)
  tmp.eval.ct$time <- as.numeric(difftime(t1, t0, units = "secs"))
  tmp.eval.ct
  
}


# tuning non-honest CT
tmp.all <- data.frame(matrix(0, nrow(CTGrid), 13))
for (iter in 1:50) {
  tmp <- foreach(i = 1:nrow(CTGrid), .combine = rbind, .export = ls(), .packages = pack) %dopar% tune.ct.hetero.nonhonest.simu.func(ind = i)
  tmp.all <- tmp.all + tmp
}
## select the optimal parameters of ct.hetero.nonhonest:
# 1st. split_rule = "tstats", split_honest = TRUE,
#      cv_option = "matching", cv_honest = NA
# 2nd. split_rule = "tstats", split_honest = TRUE,
#      cv_option = "CT", cv_honest = FALSE


################################################################################################
################################ 2.tuning CT: homogeneous, nonhonest ###########################
################################################################################################
tune.ct.homo.nonhonest.simu.func <- function(ind) {
  
  data.lst <- makeData.cont.noeff.cont(N = 1000, n.test = 1000, p = 5)
  data <- data.lst$data.used
  # data.noy <- data %>% dplyr::select(-Y)
  
  form.term <- colnames(data)[!colnames(data) %in% c("A", "Y")]
  form <- paste("A ~", paste0(form.term, collapse = " + "))
  
  t0 <- Sys.time()
  tmp.propsc <- est.pscore(form.pscore = form, data.node = data, method = "GLM", crossfit = FALSE)
  pw <- ifelse(data$A == 1, tmp.propsc, 1 - tmp.propsc)
  
  if (is.na(CTGrid$split_honest[ind])) {
    if (is.na(CTGrid$cv_honest[ind])) {
      tmp.ct.nonhonest <- causalTree::causalTree(Y ~ ., data = data, weights = 1/pw,
                                                 treatment   = data$A, 
                                                 split.Rule  = CTGrid$split_rule[ind], 
                                                 cv.option   = CTGrid$cv_option[ind], 
                                                 split.Bucket = T, xval = 10,
                                                 cp = 0, minsize = 20) 
    } else { # !is.na(cv_honest)
      tmp.ct.nonhonest <- causalTree::causalTree(Y ~ ., data = data, weights = 1/pw,
                                                 treatment   = data$A, 
                                                 split.Rule  = CTGrid$split_rule[ind], 
                                                 cv.option   = CTGrid$cv_option[ind], 
                                                 cv.Honest   = CTGrid$cv_honest[ind], 
                                                 split.Bucket = T, xval = 10,
                                                 cp = 0, minsize = 20) 
    }
    
  } else { # !is.na(split_honest)
    
    if (is.na(CTGrid$cv_honest[ind])) {
      tmp.ct.nonhonest <- causalTree::causalTree(Y ~ ., data = data, weights = 1/pw,
                                                 treatment   = data$A, 
                                                 split.Rule  = CTGrid$split_rule[ind], 
                                                 cv.option   = CTGrid$cv_option[ind], 
                                                 split.Honest = CTGrid$split_honest[ind],
                                                 split.Bucket = T, xval = 10,
                                                 cp = 0, minsize = 20) 
    } else { # !is.na(cv_honest)
      tmp.ct.nonhonest <- causalTree::causalTree(Y ~ ., data = data, weights = 1/pw,
                                                 treatment   = data$A, 
                                                 split.Rule  = CTGrid$split_rule[ind], 
                                                 cv.option   = CTGrid$cv_option[ind],
                                                 split.Honest = CTGrid$split_honest[ind],
                                                 cv.Honest   = CTGrid$cv_honest[ind], 
                                                 split.Bucket = T, xval = 10,
                                                 cp = 0, minsize = 20) 
    }
    
  }
  
  opcp  <- tmp.ct.nonhonest$cptable[,1][which.min(tmp.ct.nonhonest$cptable[,4])]
  final.ct.nonhonest <- prune(tmp.ct.nonhonest, cp = opcp)
  
  t1 <- Sys.time()
  
  X.test <- as.matrix(data.lst$test.data)
  e.test <- exp(X.test[, 1:5] %*% pscore.coeff)/(1 + exp(X.test[, 1:5] %*% pscore.coeff))
  
  tmp.eval.ct <- eval.measures.eff(final.tree   = final.ct.nonhonest,
                                   test.data    = data.frame(data.lst$test.data, A = rbinom(1000, 1, e.test)),
                                   true.trt.eff = data.lst$true.trt.eff,
                                   noise.var    = data.lst$noise.splt.var,
                                   corr.split   = data.lst$eff.modifier,
                                   type.eff.mod = data.lst$type.eff.mod,
                                   where.split  = data.lst$where.split,
                                   dir.split    = data.lst$dir.split)
  tmp.eval.ct$time <- as.numeric(difftime(t1, t0, units = "secs"))
  tmp.eval.ct
  
}


# tuning non-honest CT (foreach)
tmp.all <- data.frame(matrix(0, nrow(CTGrid), 13))
for (iter in 1:50) {
  tmp <- foreach(i = 1:nrow(CTGrid), .combine = rbind, .export = ls(), .packages = pack) %dopar% tune.ct.homo.nonhonest.simu.func(ind = i)
  tmp.all <- tmp.all + tmp
}

res.tune.ct.homo.nonhonest <- tmp.all/50
colnames(res.tune.ct.homo.nonhonest) <- colnames(tmp)
# c("fir.splt.ind","abs.err.sp","corr.tree","size.tree", "num.splt",
#   "numb.noise","pps", "mse", "prop.noise.all", "num.corr.splt.frst", 
#   "num.corr.splt.all", "prop.corr.all", "time")

# select optimal parameters for homo causal Tree
# 1st 28:
#     split_rule split_honest cv_option cv_honest
#     tstats      FALSE      matching     NA

# 2nd 25
# 3rd 5/23


################################################################################################
########################### 3. tuning CT: heterogeneous, honest ################################
################################################################################################
tune.ct.hetero.honest.simu.func <- function(ind) {
  
  data.lst <- makeData.cont.eff.cont(N = 1000, n.test = 1000, p = 5)
  data <- data.lst$data.used
  # data.noy <- data %>% dplyr::select(-Y)
  
  form.term <- colnames(data)[!colnames(data) %in% c("A", "Y")]
  form <- paste("A ~", paste0(form.term, collapse = " + "))
  
  t0 <- Sys.time()
  tmp.propsc <- est.pscore(form.pscore = form, data.node = data, method = "GLM", crossfit = FALSE)
  pw <- ifelse(data$A == 1, tmp.propsc, 1 - tmp.propsc)
  
  if (is.na(CTGrid$split_honest[ind])) {
    if (is.na(CTGrid$cv_honest[ind])) {
      tmp.ct.honest <- causalTree::honest.causalTree(Y ~ ., data = train.data, 
                                                     weights = 1/pw[train.idx, ], 
                                                     treatment = train.data$A,
                                                     est_data = est.data, 
                                                     est_weights = 1/pw[-train.idx, ], 
                                                     est_treatment = est.data$A,
                                                     split.Rule = CTGrid$split_rule[ind], 
                                                     cv.option = CTGrid$cv_option[ind], 
                                                     HonestSampleSize = nrow(est.data),
                                                     split.Bucket = T)
    } else { # !is.na(cv_honest)
      tmp.ct.honest <- causalTree::honest.causalTree(Y ~ ., data = train.data, 
                                                     weights = 1/pw[train.idx, ], 
                                                     treatment = train.data$A,
                                                     est_data = est.data, 
                                                     est_weights = 1/pw[-train.idx, ], 
                                                     est_treatment = est.data$A,
                                                     split.Rule = CTGrid$split_rule[ind], 
                                                     cv.option = CTGrid$cv_option[ind], 
                                                     cv.Honest = CTGrid$cv_honest[ind],
                                                     HonestSampleSize = nrow(est.data),
                                                     split.Bucket = T)
    }
    
  } else { # !is.na(split_honest)
    
    if (is.na(CTGrid$cv_honest[ind])) {
      tmp.ct.honest <- causalTree::honest.causalTree(Y ~ ., data = train.data, 
                                                     weights = 1/pw[train.idx, ], 
                                                     treatment = train.data$A,
                                                     est_data = est.data, 
                                                     est_weights = 1/pw[-train.idx, ], 
                                                     est_treatment = est.data$A,
                                                     split.Rule = CTGrid$split_rule[ind], 
                                                     split.Honest = CTGrid$split_honest[ind],
                                                     cv.option = CTGrid$cv_option[ind], 
                                                     HonestSampleSize = nrow(est.data),
                                                     split.Bucket = T)
    } else { # !is.na(cv_honest)
      tmp.ct.honest <- causalTree::honest.causalTree(Y ~ ., data = train.data, 
                                                     weights = 1/pw[train.idx, ], 
                                                     treatment = train.data$A,
                                                     est_data = est.data, 
                                                     est_weights = 1/pw[-train.idx, ], 
                                                     est_treatment = est.data$A,
                                                     split.Rule = CTGrid$split_rule[ind], 
                                                     split.Honest = CTGrid$split_honest[ind],
                                                     cv.option = CTGrid$cv_option[ind],
                                                     cv.Honest= CTGrid$cv_honest[ind],
                                                     HonestSampleSize = nrow(est.data),
                                                     split.Bucket = T)
    }
    
  }
  
  opcp  <- tmp.ct.honest$cptable[,1][which.min(tmp.ct.honest$cptable[,4])]
  final.ct.honest <- prune(tmp.ct.honest, cp = opcp)
  
  t1 <- Sys.time()
  
  X.test <- as.matrix(data.lst$test.data)
  e.test <- exp(X.test[, 1:5] %*% pscore.coeff)/(1 + exp(X.test[, 1:5] %*% pscore.coeff))
  
  tmp.eval.ct <- eval.measures.eff(final.tree   = final.ct.honest,
                                   test.data    = data.frame(data.lst$test.data, A = rbinom(1000, 1, e.test)),
                                   true.trt.eff = data.lst$true.trt.eff,
                                   noise.var    = data.lst$noise.splt.var,
                                   corr.split   = data.lst$eff.modifier,
                                   type.eff.mod = data.lst$type.eff.mod,
                                   where.split  = data.lst$where.split,
                                   dir.split    = data.lst$dir.split)
  tmp.eval.ct$time <- as.numeric(difftime(t1, t0, units = "secs"))
  tmp.eval.ct
  
}


# tuning honest CT 









################################################################################################
################################### 4.tuning CT: homogeneous, honest ###########################
################################################################################################
tune.ct.homo.honest.simu.func <- function(ind) {
  
  data.lst <- makeData.cont.eff.cont(N = 1000, n.test = 1000, p = 5)
  data <- data.lst$data.used
  # data.noy <- data %>% dplyr::select(-Y)
  
  form.term <- colnames(data)[!colnames(data) %in% c("A", "Y")]
  form <- paste("A ~", paste0(form.term, collapse = " + "))
  
  t0 <- Sys.time()
  tmp.propsc <- est.pscore(form.pscore = form, data.node = data, method = "GLM", crossfit = FALSE)
  pw <- ifelse(data$A == 1, tmp.propsc, 1 - tmp.propsc)
  
  if (is.na(CTGrid$split_honest[ind])) {
    if (is.na(CTGrid$cv_honest[ind])) {
      tmp.ct.honest <- causalTree::honest.causalTree(Y ~ ., data = train.data, 
                                                     weights = 1/pw[train.idx, ], 
                                                     treatment = train.data$A,
                                                     est_data = est.data, 
                                                     est_weights = 1/pw[-train.idx, ], 
                                                     est_treatment = est.data$A,
                                                     split.Rule = CTGrid$split_rule[ind], 
                                                     cv.option = CTGrid$cv_option[ind], 
                                                     HonestSampleSize = nrow(est.data),
                                                     split.Bucket = T)
    } else { # !is.na(cv_honest)
      tmp.ct.honest <- causalTree::honest.causalTree(Y ~ ., data = train.data, 
                                                     weights = 1/pw[train.idx, ], 
                                                     treatment = train.data$A,
                                                     est_data = est.data, 
                                                     est_weights = 1/pw[-train.idx, ], 
                                                     est_treatment = est.data$A,
                                                     split.Rule = CTGrid$split_rule[ind], 
                                                     cv.option = CTGrid$cv_option[ind], 
                                                     cv.Honest = CTGrid$cv_honest[ind],
                                                     HonestSampleSize = nrow(est.data),
                                                     split.Bucket = T)
    }
    
  } else { # !is.na(split_honest)
    
    if (is.na(CTGrid$cv_honest[ind])) {
      tmp.ct.honest <- causalTree::honest.causalTree(Y ~ ., data = train.data, 
                                                     weights = 1/pw[train.idx, ], 
                                                     treatment = train.data$A,
                                                     est_data = est.data, 
                                                     est_weights = 1/pw[-train.idx, ], 
                                                     est_treatment = est.data$A,
                                                     split.Rule = CTGrid$split_rule[ind], 
                                                     split.Honest = CTGrid$split_honest[ind],
                                                     cv.option = CTGrid$cv_option[ind], 
                                                     HonestSampleSize = nrow(est.data),
                                                     split.Bucket = T)
    } else { # !is.na(cv_honest)
      tmp.ct.honest <- causalTree::honest.causalTree(Y ~ ., data = train.data, 
                                                     weights = 1/pw[train.idx, ], 
                                                     treatment = train.data$A,
                                                     est_data = est.data, 
                                                     est_weights = 1/pw[-train.idx, ], 
                                                     est_treatment = est.data$A,
                                                     split.Rule = CTGrid$split_rule[ind], 
                                                     split.Honest = CTGrid$split_honest[ind],
                                                     cv.option = CTGrid$cv_option[ind],
                                                     cv.Honest= CTGrid$cv_honest[ind],
                                                     HonestSampleSize = nrow(est.data),
                                                     split.Bucket = T)
    }
    
  }
  
  opcp  <- tmp.ct.honest$cptable[,1][which.min(tmp.ct.honest$cptable[,4])]
  final.ct.honest <- prune(tmp.ct.honest, cp = opcp)
  
  t1 <- Sys.time()
  
  X.test <- as.matrix(data.lst$test.data)
  e.test <- exp(X.test[, 1:5] %*% pscore.coeff)/(1 + exp(X.test[, 1:5] %*% pscore.coeff))
  
  tmp.eval.ct <- eval.measures.eff(final.tree   = final.ct.honest,
                                   test.data    = data.frame(data.lst$test.data, A = rbinom(1000, 1, e.test)),
                                   true.trt.eff = data.lst$true.trt.eff,
                                   noise.var    = data.lst$noise.splt.var,
                                   corr.split   = data.lst$eff.modifier,
                                   type.eff.mod = data.lst$type.eff.mod,
                                   where.split  = data.lst$where.split,
                                   dir.split    = data.lst$dir.split)
  tmp.eval.ct$time <- as.numeric(difftime(t1, t0, units = "secs"))
  tmp.eval.ct
  
}


# tuning honest CT
# res.tune.ct.homo.honest <- lapply(1:nrow(CTGrid), function(i) {
#   
#   tmp <- sapply(1:50, function(rr) {tune.ct.homo.honest.simu.func(ind = i)})
#   tmp <- apply(tmp, 1, as.numeric)
#   tmp.eval.ct <- colMeans(tmp, na.rm = TRUE)
#   
# })