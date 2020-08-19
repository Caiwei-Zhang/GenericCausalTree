setwd("E:/DownloadCode/CIT-master")  #### setwd("..")
folder <- paste(getwd(), "/Functions/", sep="")
functions <- list.files(folder)
functions <- paste(folder, functions, sep = "")
for (i in functions){
  source(i)
}

## adapt primary  code for fitting outcome model by gbm
sim.reps <- 1000
set.seed(999)
pscore.coeff <- matrix(c(0.5, -0.5, 0.5, -0.5, 0.5), 5, 1)

## non-honest
################################################################################################
############################### 1. CT: heterogenous, nonhonest #################################
################################################################################################
CT.hetero.nonhonest.simu.func <- function(nn, dim = 5, unmeas.conf = F) {
  
  data.lst <- makeData.cont.eff.cont(N = nn, n.test = 1000, p = dim)

  if (unmeas.conf == F) {
    data <- data.lst$data.used
    test.data <- data.lst$test.data
  } else {
    data <- data.lst$data.used %>% dplyr::select(-X2)
    test.data <- data.lst$test.data %>% dplyr::select(-X2)
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
                                    cp = 0, minsize = 30)
  opcp  <- init.ct$cptable[,1][which.min(init.ct$cptable[,4])]
  final.ct <- prune(init.ct, cp = opcp)
  t1 <- Sys.time()
  
  # predict.rpart(): test.data must contains "A"
  X.test <- as.matrix(data.lst$test.data)
  e.test <- exp(X.test[, 1:5] %*% pscore.coeff)/(1 + exp(X.test[, 1:5] %*% pscore.coeff))
  
  eval.ct <- eval.measures.eff(final.tree   = final.ct,
                               test.data    = data.frame(test.data, A = rbinom(1000, 1, e.test)),
                               true.trt.eff = data.lst$true.trt.eff,
                               noise.var    = data.lst$noise.splt.var,
                               corr.split   = data.lst$eff.modifier,
                               type.eff.mod = data.lst$type.eff.mod,
                               where.split  = data.lst$where.split,
                               dir.split    = data.lst$dir.split)

  eval.ct$t <- as.numeric(difftime(t1, t0, units = "secs"))
  
  eval.ct
  
}


################################################################################################
############################### 2. CT: homogeneous, nonhonest ##################################
################################################################################################
CT.homo.nonhonest.simu.func <- function(nn, dim = 5, unmeas.conf = F) {
  
  data.lst <- makeData.cont.noeff.cont(N = nn, n.test = 1000, p = dim)
  
  if (unmeas.conf == F) {
    data <- data.lst$data.used
    test.data <- data.lst$test.data
  } else {
    data <- data.lst$data.used %>% dplyr::select(-X2)
    test.data <- data.lst$test.data %>% dplyr::select(-X2)
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
                                    cp = 0, minsize = 30)
  opcp  <- init.ct$cptable[,1][which.min(init.ct$cptable[,4])]
  final.ct <- prune(init.ct, cp = opcp)
  t1 <- Sys.time()
  
  X.test <- as.matrix(data.lst$test.data)
  e.test <- exp(X.test[, 1:5] %*% pscore.coeff)/(1 + exp(X.test[, 1:5] %*% pscore.coeff))
  
  eval.ct <- eval.measures.eff(final.tree   = final.ct,
                               test.data    = data.frame(test.data, A = rbinom(1000, 1, e.test)),
                               true.trt.eff = data.lst$true.trt.eff,
                               noise.var    = data.lst$noise.splt.var,
                               corr.split   = data.lst$eff.modifier,
                               type.eff.mod = data.lst$type.eff.mod,
                               where.split  = data.lst$where.split,
                               dir.split    = data.lst$dir.split)
  
  eval.ct$t <- as.numeric(difftime(t1, t0, units = "secs"))
  
  eval.ct
  
}


################################################################################################
################################ 3. CT: heterogeneous, honest ##################################
################################################################################################
CT.hetero.honest.simu.func <- function(nn, dim = 5, unmeas.conf = F) {
  
  data.lst <- makeData.cont.eff.cont(N = nn, n.test = 1000, p = dim)
  if (unmeas.conf == F) {
    data <- data.lst$data.used
    test.data <- data.lst$test.data
  } else {
    data <- data.lst$data.used %>% dplyr::select(-X2)
    test.data <- data.lst$test.data %>% dplyr::select(-X2)
  }
  
  train.idx <- createDataPartition(data$A, times = 1, p = 0.5)[[1]]
  train.data <- data[train.idx, ]
  est.data   <- data[-train.idx, ]
  
  form.term <- colnames(data)[!colnames(data) %in% c("A", "Y")]
  form <- paste("A ~", paste0(form.term, collapse = " + "))
  
  t0 <- Sys.time()
  tmp.propsc <- est.pscore(form.pscore = form, data.node = data, method = "GLM", crossfit = FALSE)
  pw <- ifelse(data$A == 1, tmp.propsc, 1 - tmp.propsc)
  
  init.ct <- causalTree::honest.causalTree(Y ~ ., data = train.data, 
                                           weights = 1/pw[train.idx, ], 
                                           treatment = train.data$A,
                                           est_data = est.data, 
                                           est_weights = 1/pw[-train.idx, ], 
                                           est_treatment = est.data$A,
                                           split.Rule = "CT", split.Honest = T, 
                                           cv.option = "CT", cv.Honest = T,
                                           HonestSampleSize = nrow(est.data),
                                           split.Bucket = T)
  
  opcp  <- init.ct$cptable[,1][which.min(init.ct$cptable[,4])]
  final.ct <- prune(init.ct, cp = opcp)
  
  t1 <- Sys.time()
  
  X.test <- as.matrix(data.lst$test.data)
  e.test <- exp(X.test[, 1:5] %*% pscore.coeff)/(1 + exp(X.test[, 1:5] %*% pscore.coeff))
  
  eval.ct <- eval.measures.eff(final.tree   = final.ct,
                               test.data    = data.frame(test.data, A = rbinom(1000, 1, e.test)),
                               true.trt.eff = data.lst$true.trt.eff,
                               noise.var    = data.lst$noise.splt.var,
                               corr.split   = data.lst$eff.modifier,
                               type.eff.mod = data.lst$type.eff.mod,
                               where.split  = data.lst$where.split,
                               dir.split    = data.lst$dir.split)
  
  eval.ct$time <- as.numeric(difftime(t1, t0, units = "secs"))
  
  eval.ct
  
}


################################################################################################
################################ 4. CT: homogeneous, honest ####################################
################################################################################################
CT.homo.honest.simu.func <- function(nn, dim = 5, unmeas.conf = F) {
  
  data.lst <- makeData.cont.noeff.cont(N = nn, n.test = 1000, p = dim)
  if (unmeas.conf == F) {
    data <- data.lst$data.used
    test.data <- data.lst$test.data
  } else {
    data <- data.lst$data.used %>% dplyr::select(-X2)
    test.data <- data.lst$test.data %>% dplyr::select(-X2)
  }
  
  train.idx <- createDataPartition(data$A, times = 1, p = 0.5)[[1]]
  train.data <- data[train.idx, ]
  est.data   <- data[-train.idx, ]
  
  form.term <- colnames(data)[!colnames(data) %in% c("A", "Y")]
  form <- paste("A ~", paste0(form.term, collapse = " + "))
  
  t0 <- Sys.time()
  tmp.propsc <- est.pscore(form.pscore = form, data.node = data, method = "GLM", crossfit = FALSE)
  pw <- ifelse(data$A == 1, tmp.propsc, 1 - tmp.propsc)
  
  init.ct <- causalTree::honest.causalTree(Y ~ ., data = train.data, 
                                           weights = 1/pw[train.idx, ], 
                                           treatment = train.data$A,
                                           est_data = est.data, 
                                           est_weights = 1/pw[-train.idx, ], 
                                           est_treatment = est.data$A,
                                           split.Rule = "CT", split.Honest = T, 
                                           cv.option = "CT", cv.Honest = T,
                                           HonestSampleSize = nrow(est.data),
                                           split.Bucket = T)
  
  opcp  <- init.ct$cptable[,1][which.min(init.ct$cptable[,4])]
  final.ct <- prune(init.ct, cp = opcp)
  
  t1 <- Sys.time()
  
  X.test <- as.matrix(data.lst$test.data)
  e.test <- exp(X.test[, 1:5] %*% pscore.coeff)/(1 + exp(X.test[, 1:5] %*% pscore.coeff))
  
  eval.ct <- eval.measures.eff(final.tree   = final.ct,
                               test.data    = data.frame(test.data, A = rbinom(1000, 1, e.test)),
                               true.trt.eff = data.lst$true.trt.eff,
                               noise.var    = data.lst$noise.splt.var,
                               corr.split   = data.lst$eff.modifier,
                               type.eff.mod = data.lst$type.eff.mod,
                               where.split  = data.lst$where.split,
                               dir.split    = data.lst$dir.split)
  
  eval.ct$time <- as.numeric(difftime(t1, t0, units = "secs"))
  
  eval.ct
  
}


# sapply 
CT.hetero.res.nonhonest <- apply(sapply(1:sim.reps, function(iter) CT.hetero.nonhonest.simu.func(n = 1000)), 1, as.numeric)
CT.homo.res.nonhonest   <- apply(sapply(1:sim.reps, function(iter) CT.hetero.nonhonest.simu.func(n = 1000)), 1, as.numeric)

CT.hetero.fitBefore <- colMeans(CT.hetero.res.nonhonest, na.rm = TRUE)
CT.homo.fitBefore   <- colMeans(CT.homo.res.nonhonest, na.rm = TRUE)

###################################################################################################
####################   Task 1:  foreach: CausalTree [Athey et al.(2016)]    #######################
###################################################################################################
clnum <- detectCores()
cl <- makeCluster(getOption("cl.cores", clnum - 1))
registerDoParallel(cl)
pack <- c("caret", "devtools", "causalTree", "causalForest", "dplyr", "grf",
          "gbm", "MASS", "rpart", "randomForestSRC", "rattle", "rpart.plot")

CT.hetero.res.nonhonest <- foreach(1:sim.reps, .combine = rbind, 
                                   .export = ls(), .packages = pack) %dopar% CT.hetero.nonhonest.simu.func(nn = 1000)


CT.homo.res.nonhonest <- foreach(1:sim.reps, .combine = rbind, 
                                 .export = ls(), .packages = pack) %dopar% CT.homo.nonhonest.simu.func(nn = 1000)


CT.hetero.nonhonest <- colMeans(CT.hetero.res.nonhonest, na.rm = TRUE)
CT.homo.nonhonest   <- colMeans(CT.homo.res.nonhonest, na.rm = TRUE)


stopImplicitCluster()









################################################################################################
############################## 1. opt.CT: heterogeneous, nonhonest #############################
################################################################################################
# split_rule = "tstats, split_honest = TRUE,
# cv_option = "matching", cv_honest = NA
opt.ct.hetero.nonhonest.simu.func <- function(nn, dim = 5, unmeas.conf = F) {
  
  data.lst <- makeData.cont.eff.cont(N = nn, n.test = 1000, p = dim)
  
  if (unmeas.conf == F) {
    data <- data.lst$data.used
    test.data <- data.lst$test.data
  } else {
    data <- data.lst$data.used %>% dplyr::select(-X2)
    test.data <- data.lst$test.data %>% dplyr::select(-X2)
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
                                             cp = 0, minsize = 30)
  
  opcp  <- tmp.ct.nonhonest$cptable[,1][which.min(tmp.ct.nonhonest$cptable[,4])]
  final.ct.nonhonest <- prune(tmp.ct.nonhonest, cp = opcp)
  
  t1 <- Sys.time()
  
  X.test <- as.matrix(data.lst$test.data)
  e.test <- exp(X.test[, 1:5] %*% pscore.coeff)/(1 + exp(X.test[, 1:5] %*% pscore.coeff))
  
  tmp.eval <- eval.measures.eff(final.tree   = final.ct.nonhonest,
                               test.data    = data.frame(test.data, A = rbinom(1000, 1, e.test)),
                               true.trt.eff = data.lst$true.trt.eff,
                               noise.var    = data.lst$noise.splt.var,
                               corr.split   = data.lst$eff.modifier,
                               type.eff.mod = data.lst$type.eff.mod,
                               where.split  = data.lst$where.split,
                                dir.split    = data.lst$dir.split)
  tmp.eval$time <- as.numeric(difftime(t1, t0, units = "secs"))
  
  tmp.eval
  
}




################################################################################################
############################## 2. opt.CT: homogeneous, nonhonest ###############################
################################################################################################
# opt.option.homo.nonhonest <- CTGrid[opt.ind.homo.nonhonest, ]
opt.ct.homo.nonhonest.simu.func <- function(nn, dim = 5, unmeas.conf = F) {
  
  data.lst <- makeData.cont.noeff.cont(N = nn, n.test = 1000, p = dim)

  if (unmeas.conf == F) {
    data <- data.lst$data.used
    test.data <- data.lst$test.data
  } else {
    data <- data.lst$data.used %>% dplyr::select(-X2)
    test.data <- data.lst$test.data %>% dplyr::select(-X2)
  }
  
  form.term <- colnames(data)[!colnames(data) %in% c("A", "Y")]
  form <- paste("A ~", paste0(form.term, collapse = " + "))
  
  t0 <- Sys.time()

  tmp.propsc <- est.pscore(form.pscore = form, data.node = data, method = "GLM", crossfit = FALSE)
  pw <- ifelse(data$A == 1, tmp.propsc, 1 - tmp.propsc)
  
  tmp.ct.nonhonest <- causalTree::causalTree(Y ~ ., data = data, treatment = data$A, weights = 1 / pw,
                                             split.Rule  = "tstats", cv.option = "matching", 
                                             split.Honest = F, split.Bucket = T, xval = 10,
                                             cp = 0, minsize = 30)
  
  opcp  <- tmp.ct.nonhonest$cptable[,1][which.min(tmp.ct.nonhonest$cptable[,4])]
  final.ct.nonhonest <- prune(tmp.ct.nonhonest, cp = opcp)
  
  t1 <- Sys.time()
  
  X.test <- as.matrix(data.lst$test.data)
  e.test <- exp(X.test[, 1:5] %*% pscore.coeff)/(1 + exp(X.test[, 1:5] %*% pscore.coeff))
  
  tmp.eval <- eval.measures.eff(final.tree   = final.ct.nonhonest,
                                test.data    = data.frame(test.data, A = rbinom(1000, 1, e.test)),
                                true.trt.eff = data.lst$true.trt.eff,
                                noise.var    = data.lst$noise.splt.var,
                                corr.split   = data.lst$eff.modifier,
                                type.eff.mod = data.lst$type.eff.mod,
                                where.split  = data.lst$where.split,
                                dir.split    = data.lst$dir.split)
  
  tmp.eval$time <- as.numeric(difftime(t1, t0, units = "secs"))
  
  tmp.eval
  
}




################################################################################################
############################## 3. opt.CT: heterogeneous, honest ################################
################################################################################################
opt.ct.hetero.honest.simu.func <- function(nn, dim = 5, unmeas.conf = F) {
  
  data.lst <- makeData.cont.eff.cont(N = nn, n.test = 1000, p = dim)

  if (unmeas.conf == F) {
    data <- data.lst$data.used
    test.data <- data.lst$test.data
  } else {
    data <- data.lst$data.used %>% dplyr::select(-X2)
    test.data <- data.lst$test.data %>% dplyr::select(-X2)
  }
  
  train.idx  <- createDataPartition(data$A, times = 1, p = 0.5)[[1]]
  train.data <- data[train.idx, ]
  est.data   <- data[-train.idx, ]
  
  form.term <- colnames(data)[!colnames(data) %in% c("A", "Y")]
  form <- paste("A ~", paste0(form.term, collapse = " + "))
  
  t0 <- Sys.time()

  tmp.propsc <- est.pscore(form.pscore = form, data.node = data, method = "GBM", crossfit = FALSE)
  pw <- ifelse(data$A == 1, tmp.propsc, 1 - tmp.propsc)
  
  tmp.ct.honest <- honest.causalTree(Y ~ ., data = train.data, weights = 1/pw[train.idx, ], treatment = train.data$A,
                                     est_data = est.data, est_weights = 1/pw[-train.idx, ], est_treatment = est.data$A,
                                     split.Rule = "tstats", split.Honest = T, 
                                     cv.option = "matching", # cv.Honest = F,
                                     HonestSampleSize = nrow(est.data),
                                     split.Bucket = F)
  
  opcp  <- tmp.ct.honest$cptable[,1][which.min(tmp.ct.honest$cptable[,4])]
  final.ct.honest <- prune(tmp.ct.honest, cp = opcp)
  
  t1 <- Sys.time()
  
  X.test <- as.matrix(data.lst$test.data)
  e.test <- exp(X.test[, 1:5] %*% pscore.coeff)/(1 + exp(X.test[, 1:5] %*% pscore.coeff))
  
  tmp.eval <- eval.measures.eff(final.tree   = final.ct.honest,
                                test.data    = data.frame(test.data, A = rbinom(1000, 1, e.test)),
                                true.trt.eff = data.lst$true.trt.eff,
                                noise.var    = data.lst$noise.splt.var,
                                corr.split   = data.lst$eff.modifier,
                                type.eff.mod = data.lst$type.eff.mod,
                                where.split  = data.lst$where.split,
                                dir.split    = data.lst$dir.split)
  
  tmp.eval$time <- as.numeric(difftime(t1, t0, units = "secs"))
  
  tmp.eval
  
}




################################################################################################
############################## 4. opt.CT: homogeneous, honest ##################################
################################################################################################

opt.ct.homo.honest.simu.func <- function(nn, dim = 5, unmeas.conf = F) {
  
  data.lst <- makeData.cont.noeff.cont(N = nn, n.test = 1000, p = dim)

  if (unmeas.conf == F) {
    data <- data.lst$data.used
    test.data <- data.lst$test.data
  } else {
    data <- data.lst$data.used %>% dplyr::select(-X2)
    test.data <- data.lst$test.data %>% dplyr::select(-X2)
  }
  
  train.idx  <- createDataPartition(data$A, times = 1, p = 0.5)[[1]]
  train.data <- data[train.idx, ]
  est.data   <- data[-train.idx, ]
  
  form.term <- colnames(data)[!colnames(data) %in% c("A", "Y")]
  form <- paste("A ~", paste0(form.term, collapse = " + "))
                
  t0 <- Sys.time()
  
  tmp.propsc <- est.pscore(form.pscore = form, data.node = data, method = "GBM", crossfit = FALSE)
  pw <- ifelse(data$A == 1, tmp.propsc, 1 - tmp.propsc)
  
  tmp.ct.honest <- honest.causalTree(Y ~ ., data = train.data, weights = 1/pw[train.idx, ], treatment = train.data$A,
                                     est_data = est.data, est_weights = 1/pw[-train.idx, ], est_treatment = est.data$A,
                                     split.Rule = "tstats", split.Honest = F, 
                                     cv.option = "matching", # cv.Honest = F,
                                     HonestSampleSize = nrow(est.data),
                                     split.Bucket = F)
  
  opcp  <- tmp.ct.honest$cptable[,1][which.min(tmp.ct.honest$cptable[,4])]
  final.ct.honest <- prune(tmp.ct.honest, cp = opcp)
  
  t1 <- Sys.time()
  
  X.test <- as.matrix(data.lst$test.data)
  e.test <- exp(X.test[, 1:5] %*% pscore.coeff)/(1 + exp(X.test[, 1:5] %*% pscore.coeff))
  
  tmp.eval <- eval.measures.eff(final.tree   = final.ct.honest,
                                test.data    = data.frame(test.data, A = rbinom(1000, 1, e.test)),
                                true.trt.eff = data.lst$true.trt.eff,
                                noise.var    = data.lst$noise.splt.var,
                                corr.split   = data.lst$eff.modifier,
                                type.eff.mod = data.lst$type.eff.mod,
                                where.split  = data.lst$where.split,
                                dir.split    = data.lst$dir.split)
  
  tmp.eval$time <- as.numeric(difftime(t1, t0, units = "secs"))
  
  tmp.eval
  
}



set.seed(999)
opt.CT.hetero.res.nonhonest <- apply(sapply(1:sim.reps, function(iter) opt.ct.hetero.nonhonest.simu.func(n = 1000)), 1, as.numeric)
opt.CT.homo.res.nonhonest   <- apply(sapply(1:sim.reps, function(iter) opt.ct.homo.nonhonest.simu.func(n = 1000)), 1, as.numeric)

opt.CT.hetero.nonhonest   <- colMeans(opt.CT.hetero.res.nonhonest, na.rm = TRUE)
opt.CT.homo.nonhonest     <- colMeans(opt.CT.homo.res.nonhonest, na.rm = TRUE)



# foreach
clnum <- detectCores()
cl <- makeCluster(getOption("cl.cores", clnum - 1))
registerDoParallel(cl)
pack <- c("caret", "devtools", "causalTree", "causalForest", "dplyr", "grf",
          "gbm", "MASS", "rpart", "randomForestSRC", "rattle", "rpart.plot")

opt.ct.hetero.res.nonhonest <- foreach(1:sim.reps, .combine = rbind, 
                                    .export = ls(), .packages = pack) %dopar% opt.ct.hetero.nonhonest.simu.func(nn = 1000)


opt.ct.homo.res.nonhonest <- foreach(1:sim.reps, .combine = rbind, 
                                     .export = ls(), .packages = pack) %dopar% opt.ct.homo.nonhonest.simu.func(nn = 1000)


opt.ct.hetero.res.honest <- foreach(1:sim.reps, .combine = rbind, 
                                        .export = ls(), .packages = pack) %dopar% opt.ct.hetero.honest.simu.func(nn = 1000)


opt.ct.homo.res.honest <- foreach(1:sim.reps, .combine = rbind, 
                                      .export = ls(), .packages = pack) %dopar% opt.ct.homo.honest.simu.func(nn = 1000)


opt.CT.hetero.nonhonest <- colMeans(opt.ct.hetero.res.nonhonest, na.rm = TRUE)
opt.CT.homo.nonhonest   <- colMeans(opt.ct.homo.res.nonhonest, na.rm = TRUE)
opt.CT.hetero.honest <- colMeans(opt.ct.hetero.res.honest, na.rm = TRUE)
opt.CT.homo.honest   <- colMeans(opt.ct.homo.res.honest, na.rm = TRUE)

stopImplicitCluster()

save.image(file = "E:/бшнд/groupTest/MethodCompare/simulation/simu.CT.RData")

save(CT.hetero.res.nonhonest, CT.homo.res.nonhonest, 
     opt.CT.hetero.res.nonhonest, opt.CT.homo.res.nonhonest, 
     CT.hetero.nonhonest, CT.homo.nonhonest, 
     opt.CT.hetero.nonhonest, opt.CT.homo.nonhonest,  
     file = "E:/project/groupTest/MethodCompare/simulation/res.CT.RData")



     # CT.hetero.res.honest, CT.homo.res.honest, CT.hetero.honest, CT.homo.honest,
     # opt.ct.hetero.res.honest, opt.ct.homo.res.honest,
     # opt.ct.hetero.honest, opt.ct.homo.honest,

