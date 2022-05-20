############### CausalTree [Athey et al.(2016)] #################
setwd("..")
load("./seed.RData")
source("./Simulations/Study1_MethodComparison/CIT-master/Functions/EvalMeas.R")
source("./Functions/dgp.R")



################################################################################################
############################### 1. CT: heterogenous, nonhonest #################################
CT.hetero.nonhonest.simu.func <- function(nn, dim = 5, unmeas.conf = F, seed = NULL) {
  
  if(!is.null(seed)){
    set.seed(seed)
  }
  
  data.lst <- hetero.DataGen(N = nn, p = dim, dgp = "cov", N.test = 1000)
  data.lst$where.split <- list(c(1))   
  data.lst$dir.split  <- list(c(NULL))
  
  if (unmeas.conf == F) {
    data <- data.lst$data
    data.test <- data.lst$data.test
  } else {
    data <- data.lst$data %>% dplyr::select(-X2)
    data.test <- data.lst$data.test %>% dplyr::select(-X2)
  }
  
  form.term <- colnames(data)[!colnames(data) %in% c("A", "Y")]
  form <- paste("A ~", paste0(form.term, collapse = " + "))
  
  t0 <- Sys.time()
  tmp.propsc <- est.pscore(form.pscore = form, data.node = data, method = "GLM", crossfit = FALSE)
  pw <- ifelse(data$A == 1, tmp.propsc, 1 - tmp.propsc)
  
  init.ct <- causalTree::causalTree(Y ~ ., data = data, treatment = data$A, weights = 1 / pw,
                                    split.Rule  = "CT", cv.option = "CT", 
                                    split.Honest = T, cv.Honest = T, 
                                    split.Bucket = T, xval = 10,
                                    cp = 0, minsize = 10)# nn/200
  opcp  <- init.ct$cptable[,1][which.min(init.ct$cptable[,4])]
  final.ct <- prune(init.ct, cp = opcp)
  t1 <- Sys.time()
  
  X.test <- as.matrix(data.test)
  pscore.coeff <- matrix(c(0.5, -0.5, 0.5, -0.5, 0.5), 5, 1)
  if (unmeas.conf == F) {
    e.test <- exp(X.test[, 1:dim] %*% pscore.coeff)/(1 + exp(X.test[, 1:dim] %*% pscore.coeff))
  } else {
    e.test <- exp(X.test[, 1:(dim-1)] %*% pscore.coeff[-2, , drop = FALSE])/(1 + exp(X.test[, 1:(dim-1)] %*% pscore.coeff[-2, , drop = FALSE]))
  }
  
  
  eval.ct <- eval.measures.eff(final.tree   = final.ct,
                               test.data    = data.frame(data.test, A = rbinom(1000, 1, e.test)),
                               true.trt.eff = data.lst$true.trt.eff,
                               noise.var    = data.lst$noise.splt.var,
                               corr.split   = data.lst$eff.modifier,
                               where.split  = data.lst$where.split,
                               dir.split    = data.lst$dir.split)
  
  eval.ct$t <- as.numeric(difftime(t1, t0, units = "secs"))
  
  return(unlist(eval.ct))
  
}


################################################################################################
############################### 2. CT: homogeneous, nonhonest ##################################
################################################################################################
CT.homo.nonhonest.simu.func <- function(nn, dim = 5, unmeas.conf = F, seed = NULL) {
  
  if(!is.null(seed)){
    set.seed(seed)
  }
  
  data.lst <- homo.DataGen(N = nn, p = dim, dgp = "cov", N.test = 1000)
  data.lst$where.split <- list(c(1))   
  data.lst$dir.split  <- list(c(NULL))
  
  if (unmeas.conf == F) {
    data <- data.lst$data
    data.test <- data.lst$data.test
  } else {
    data <- data.lst$data %>% dplyr::select(-X2)
    data.test <- data.lst$data.test %>% dplyr::select(-X2)
  }
  
  form.term <- colnames(data)[!colnames(data) %in% c("A", "Y")]
  form <- paste("A ~", paste0(form.term, collapse = " + "))
  
  t0 <- Sys.time()
  tmp.propsc <- est.pscore(form.pscore = form, data.node = data, method = "GLM", crossfit = FALSE)
  pw <- ifelse(data$A == 1, tmp.propsc, 1 - tmp.propsc)
  
  init.ct <- causalTree::causalTree(Y ~ ., data = data, treatment = data$A, weights = 1 / pw,
                                    split.Rule  = "CT", cv.option = "CT", 
                                    split.Honest = T, cv.Honest = T, 
                                    split.Bucket = T, xval = 10,
                                    cp = 0, minsize = 50)
  opcp  <- init.ct$cptable[,1][which.min(init.ct$cptable[,4])]
  final.ct <- prune(init.ct, cp = opcp)
  t1 <- Sys.time()
  
  X.test <- as.matrix(data.test)
  pscore.coeff <- matrix(c(0.5, -0.5, 0.5, -0.5, 0.5), 5, 1)
  if (unmeas.conf == F) {
    e.test <- exp(X.test[, 1:dim] %*% pscore.coeff)/(1 + exp(X.test[, 1:dim] %*% pscore.coeff))
  } else {
    e.test <- exp(X.test[, 1:(dim-1)] %*% pscore.coeff[-2, , drop = FALSE])/(1 + exp(X.test[, 1:(dim-1)] %*% pscore.coeff[-2, , drop = FALSE]))
  }
  
  eval.ct <- eval.measures.eff(final.tree   = final.ct,
                               test.data    = data.frame(data.test, A = rbinom(1000, 1, e.test)),
                               true.trt.eff = data.lst$true.trt.eff,
                               noise.var    = data.lst$noise.splt.var,
                               corr.split   = data.lst$eff.modifier,
                               where.split  = data.lst$where.split,
                               dir.split    = data.lst$dir.split)
  
  eval.ct$t <- as.numeric(difftime(t1, t0, units = "secs"))
  
  return(unlist(eval.ct))
  
}


# ################################################################################################
# ################################ 3. CT: heterogeneous, honest ##################################
# ################################################################################################
# CT.hetero.honest.simu.func <- function(nn, dim = 5, unmeas.conf = F) {
#   
#   data.lst <- hetero.DataGen(N = nn, p = dim, dgp = "cov", N.test = 1000)
#   data.lst$where.split <- list(c(1))   
#   data.lst$dir.split  <- list(c(NULL))
#   
#   if (unmeas.conf == F) {
#     data <- data.lst$data
#     data.test <- data.lst$data.test
#   } else {
#     data <- data.lst$data %>% dplyr::select(-X2)
#     data.test <- data.lst$data.test %>% dplyr::select(-X2)
#   }
#   
#   train.idx <- createDataPartition(data$A, times = 1, p = 0.5)[[1]]
#   train.data <- data[train.idx, ]
#   est.data   <- data[-train.idx, ]
#   
#   form.term <- colnames(data)[!colnames(data) %in% c("A", "Y")]
#   form <- paste("A ~", paste0(form.term, collapse = " + "))
#   
#   t0 <- Sys.time()
#   tmp.propsc <- est.pscore(form.pscore = form, data.node = data, method = "GLM", crossfit = FALSE)
#   pw <- ifelse(data$A == 1, tmp.propsc, 1 - tmp.propsc)
#   
#   init.ct <- causalTree::honest.causalTree(Y ~ ., data = train.data, 
#                                            weights = 1/pw[train.idx, ], 
#                                            treatment = train.data$A,
#                                            est_data = est.data, 
#                                            est_weights = 1/pw[-train.idx, ], 
#                                            est_treatment = est.data$A,
#                                            split.Rule = "CT", split.Honest = T, 
#                                            cv.option = "CT", cv.Honest = T,
#                                            HonestSampleSize = nrow(est.data),
#                                            split.Bucket = T)
#   
#   opcp  <- init.ct$cptable[,1][which.min(init.ct$cptable[,4])]
#   final.ct <- prune(init.ct, cp = opcp)
#   
#   t1 <- Sys.time()
#   
#   X.test <- as.matrix(data.lst$test.data)
#   pscore.coeff <- matrix(c(0.5, -0.5, 0.5, -0.5, 0.5), 5, 1)
#   e.test <- exp(X.test[, 1:5] %*% pscore.coeff)/(1 + exp(X.test[, 1:5] %*% pscore.coeff))
#   
#   eval.ct <- eval.measures.eff(final.tree   = final.ct,
#                                test.data    = data.frame(test.data, A = rbinom(1000, 1, e.test)),
#                                true.trt.eff = data.lst$true.trt.eff,
#                                noise.var    = data.lst$noise.splt.var,
#                                corr.split   = data.lst$eff.modifier,
#                                where.split  = data.lst$where.split,
#                                dir.split    = data.lst$dir.split)
#   
#   eval.ct$time <- as.numeric(difftime(t1, t0, units = "secs"))
#   
#   eval.ct
#   
# }
# 
# 
# ################################################################################################
# ################################ 4. CT: homogeneous, honest ####################################
# ################################################################################################
# CT.homo.honest.simu.func <- function(nn, dim = 5, unmeas.conf = F) {
#   
#   data.lst <- makeData.cont.noeff.cont(N = nn, n.test = 1000, p = dim)
#   if (unmeas.conf == F) {
#     data <- data.lst$data.used
#     test.data <- data.lst$test.data
#   } else {
#     data <- data.lst$data.used %>% dplyr::select(-X2)
#     test.data <- data.lst$test.data %>% dplyr::select(-X2)
#   }
#   
#   train.idx <- createDataPartition(data$A, times = 1, p = 0.5)[[1]]
#   train.data <- data[train.idx, ]
#   est.data   <- data[-train.idx, ]
#   
#   form.term <- colnames(data)[!colnames(data) %in% c("A", "Y")]
#   form <- paste("A ~", paste0(form.term, collapse = " + "))
#   
#   t0 <- Sys.time()
#   tmp.propsc <- est.pscore(form.pscore = form, data.node = data, method = "GLM", crossfit = FALSE)
#   pw <- ifelse(data$A == 1, tmp.propsc, 1 - tmp.propsc)
#   
#   init.ct <- causalTree::honest.causalTree(Y ~ ., data = train.data, 
#                                            weights = 1/pw[train.idx, ], 
#                                            treatment = train.data$A,
#                                            est_data = est.data, 
#                                            est_weights = 1/pw[-train.idx, ], 
#                                            est_treatment = est.data$A,
#                                            split.Rule = "CT", split.Honest = T, 
#                                            cv.option = "CT", cv.Honest = T,
#                                            HonestSampleSize = nrow(est.data),
#                                            split.Bucket = T)
#   
#   opcp  <- init.ct$cptable[,1][which.min(init.ct$cptable[,4])]
#   final.ct <- prune(init.ct, cp = opcp)
#   
#   t1 <- Sys.time()
#   
#   X.test <- as.matrix(data.lst$test.data)
#   pscore.coeff <- matrix(c(0.5, -0.5, 0.5, -0.5, 0.5), 5, 1)
#   e.test <- exp(X.test[, 1:5] %*% pscore.coeff)/(1 + exp(X.test[, 1:5] %*% pscore.coeff))
#   
#   eval.ct <- eval.measures.eff(final.tree   = final.ct,
#                                test.data    = data.frame(test.data, A = rbinom(1000, 1, e.test)),
#                                true.trt.eff = data.lst$true.trt.eff,
#                                noise.var    = data.lst$noise.splt.var,
#                                corr.split   = data.lst$eff.modifier,
#                                where.split  = data.lst$where.split,
#                                dir.split    = data.lst$dir.split)
#   
#   eval.ct$time <- as.numeric(difftime(t1, t0, units = "secs"))
#   
#   eval.ct
#   
# }

# foreach
clnum <- detectCores()
cl <- makeCluster(getOption("cl.cores", clnum - 4))
registerDoParallel(cl)
packCT <- c("caret", "devtools", "causalTree", "dplyr", "gbm", "MASS", "rpart")

CT.hetero.res.nonhonest <- foreach(i = 1:sim.reps, .combine = rbind, .export = ls(), .packages = packCT) %dopar% CT.hetero.nonhonest.simu.func(nn = 1000, seed = seed[i])
CT.homo.res.nonhonest   <- foreach(i = 1:sim.reps, .combine = rbind, .export = ls(), .packages = packCT) %dopar% CT.homo.nonhonest.simu.func(nn = 1000, seed = seed[i+sim.reps])

CT.hetero.nonhonest <- colMeans(CT.hetero.res.nonhonest, na.rm = TRUE)
CT.homo.nonhonest   <- colMeans(CT.homo.res.nonhonest, na.rm = TRUE)


save(CT.hetero.res.nonhonest,
     CT.homo.res.nonhonest,  
     CT.hetero.nonhonest, 
     CT.homo.nonhonest,
     file = "./Simulations/Study1_MethodComparison/res/CT.1000.RData")


#stopImplicitCluster()















################################################################################################
############################## 1. opt.CT: heterogeneous, nonhonest #############################
################################################################################################
# split_rule = "tstats, split_honest = TRUE,
# cv_option = "matching", cv_honest = NA
opt.ct.hetero.nonhonest.simu.func <- function(nn, dim = 5, unmeas.conf = F, seed = NULL) {
  
  if(!is.null(seed)){
    set.seed(seed)
  }
  
  data.lst <- hetero.DataGen(N = nn, p = dim, dgp = "cov", N.test = 1000)
  data.lst$where.split <- list(c(1))   
  data.lst$dir.split  <- list(c(NULL))
  
  if (unmeas.conf == F) {
    data <- data.lst$data
    data.test <- data.lst$data.test
  } else {
    data <- data.lst$data %>% dplyr::select(-X2)
    data.test <- data.lst$data.test %>% dplyr::select(-X2)
  }
  
  form.term <- colnames(data)[!colnames(data) %in% c("A", "Y")]
  form <- paste("A ~", paste0(form.term, collapse = " + "))
  
  t0 <- Sys.time()
  tmp.propsc <- est.pscore(form.pscore = form, data.node = data, method = "GLM", crossfit = FALSE)
  pw <- ifelse(data$A == 1, tmp.propsc, 1 - tmp.propsc)
  
  tmp.ct.nonhonest <- causalTree::causalTree(Y ~ ., data = data,  treatment = data$A, weights = 1 / pw,
                                             split.Rule  = "tstats", cv.option = "matching", 
                                             split.Honest = T, # cv.Honest = FALSE, 
                                             split.Bucket = T, xval = 10,
                                             cp = 0, minsize = 10) # nn/200
  
  opcp  <- tmp.ct.nonhonest$cptable[,1][which.min(tmp.ct.nonhonest$cptable[,4])]
  final.ct.nonhonest <- prune(tmp.ct.nonhonest, cp = opcp)
  
  t1 <- Sys.time()
  
  X.test <- as.matrix(data.test)
  pscore.coeff <- matrix(c(0.5, -0.5, 0.5, -0.5, 0.5), 5, 1)
  if (unmeas.conf == F) {
    e.test <- exp(X.test[, 1:dim] %*% pscore.coeff)/(1 + exp(X.test[, 1:dim] %*% pscore.coeff))
  } else {
    e.test <- exp(X.test[, 1:(dim-1)] %*% pscore.coeff[-2, , drop = FALSE])/(1 + exp(X.test[, 1:(dim-1)] %*% pscore.coeff[-2, , drop = FALSE]))
  }
  
  tmp.eval <- eval.measures.eff(final.tree   = final.ct.nonhonest,
                                test.data    = data.frame(data.test, A = rbinom(1000, 1, e.test)),
                                true.trt.eff = data.lst$true.trt.eff,
                                noise.var    = data.lst$noise.splt.var,
                                corr.split   = data.lst$eff.modifier,
                                where.split  = data.lst$where.split,
                                dir.split    = data.lst$dir.split)
  tmp.eval$time <- as.numeric(difftime(t1, t0, units = "secs"))
  
  return(unlist(tmp.eval))
  
}




################################################################################################
############################## 2. opt.CT: homogeneous, nonhonest ###############################
################################################################################################
opt.ct.homo.nonhonest.simu.func <- function(nn, dim = 5, unmeas.conf = F, seed = NULL) {
  
  if(!is.null(seed)){
    set.seed(seed)
  }
  
  data.lst <- homo.DataGen(N = nn, p = dim, dgp = "cov", N.test = 1000)
  data.lst$where.split <- list(c(1))   
  data.lst$dir.split  <- list(c(NULL))
  
  if (unmeas.conf == F) {
    data <- data.lst$data
    data.test <- data.lst$data.test
  } else {
    data <- data.lst$data %>% dplyr::select(-X2)
    data.test <- data.lst$data.test %>% dplyr::select(-X2)
  }
  
  
  form.term <- colnames(data)[!colnames(data) %in% c("A", "Y")]
  form <- paste("A ~", paste0(form.term, collapse = " + "))
  
  t0 <- Sys.time()
  
  tmp.propsc <- est.pscore(form.pscore = form, data.node = data, method = "GLM", crossfit = FALSE)
  pw <- ifelse(data$A == 1, tmp.propsc, 1 - tmp.propsc)
  
  tmp.ct.nonhonest <- causalTree::causalTree(Y ~ ., data = data, treatment = data$A, weights = 1 / pw,
                                             split.Rule  = "tstats", cv.option = "matching", 
                                             split.Honest = F, split.Bucket = T, xval = 10,
                                             cp = 0, minsize = 50)
  
  opcp  <- tmp.ct.nonhonest$cptable[,1][which.min(tmp.ct.nonhonest$cptable[,4])]
  final.ct.nonhonest <- prune(tmp.ct.nonhonest, cp = opcp)
  
  t1 <- Sys.time()
  
  X.test <- as.matrix(data.test)
  pscore.coeff <- matrix(c(0.5, -0.5, 0.5, -0.5, 0.5), 5, 1)
  if (unmeas.conf == F) {
    e.test <- exp(X.test[, 1:dim] %*% pscore.coeff)/(1 + exp(X.test[, 1:dim] %*% pscore.coeff))
  } else {
    e.test <- exp(X.test[, 1:(dim-1)] %*% pscore.coeff[-2, , drop = FALSE])/(1 + exp(X.test[, 1:(dim-1)] %*% pscore.coeff[-2, , drop = FALSE]))
  }
  
  tmp.eval <- eval.measures.eff(final.tree   = final.ct.nonhonest,
                                test.data    = data.frame(data.test, A = rbinom(1000, 1, e.test)),
                                true.trt.eff = data.lst$true.trt.eff,
                                noise.var    = data.lst$noise.splt.var,
                                corr.split   = data.lst$eff.modifier,
                                where.split  = data.lst$where.split,
                                dir.split    = data.lst$dir.split)
  
  tmp.eval$time <- as.numeric(difftime(t1, t0, units = "secs"))
  
  return(unlist(tmp.eval))
  
}



# foreach
clnum <- detectCores()
cl <- makeCluster(getOption("cl.cores", clnum - 1))
registerDoParallel(cl)

opt.ct.hetero.res.nonhonest <- foreach(i = 1:sim.reps, .combine = rbind, .export = ls(), .packages = packCT) %dopar% 
  opt.ct.hetero.nonhonest.simu.func(nn = 1000, seed = seed[i])
opt.ct.homo.res.nonhonest   <- foreach(i = 1:sim.reps, .combine = rbind, .export = ls(), .packages = packCT) %dopar% 
  opt.ct.homo.nonhonest.simu.func(nn = 1000, seed = seed[i + sim.reps])

opt.CT.hetero.nonhonest <- colMeans(opt.ct.hetero.res.nonhonest, na.rm = TRUE)
opt.CT.homo.nonhonest   <- colMeans(opt.ct.homo.res.nonhonest, na.rm = TRUE)


save(opt.CT.hetero.nonhonest, opt.CT.homo.nonhonest, 
     opt.ct.hetero.res.nonhonest, opt.ct.homo.res.nonhonest,
     file = "./Simulations/Study1_MethodComparison/res/opt.CT.1000.RData")

stopImplicitCluster()




