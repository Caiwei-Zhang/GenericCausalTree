setwd(".\\Functions")
functions <- list.files(getwd())
functions <- paste(paste(getwd(), "\\", sep = ""), functions, sep = "")
for (f in functions) {
  source(i)
}
# source("library.R")
# source("integr.model.fit.R")
# source("gbm.predict.R")
# source("SearchChilds.R")
# source("ulist.R")
# source("dgp.R")
# source("Treebuild.R")
# source("Treebuild.honest.R")
# source("eval.tree.R")
# source("eval.honest.tree.R")

sim.reps <- 1000
set.seed(999)
##############################################################################################
################## 1. heterogenuous: honest, modelfitBefore, est.fitBefore  ##################
##############################################################################################
hetero.fitBefore.honest.simu.func <- function(nn, dim = 5, unmeas.conf = F) {
  
  data.lst <- hetero.DataGen(N = nn, p = dim, dgp = "cov", N.test = 1000)
  if (unmeas.conf == F) {
    data <- data.lst$data
    data.test <- data.lst$data.test
  } else {
    data <- data.lst$data %>% dplyr::select(-X2)
    data.test <- data.lst$data.test %>% dplyr::select(-X2)
  }
  
  t0 <- Sys.time()
  tree.lst <- honest.Treebuild(data = data, form.outcome = "Y ~ .", form.pscore = "A ~ .", type.outcome = "continuous", 
                               train.frac = 0.8, modelFitBefore = TRUE, est.modelFitBefore = TRUE,
                               outcome.mthd = "GBM", pscore.mthd = "GLM", cp = qchisq(0.95, 1), 
                               minsplit = 30, minsize = 10)
  
  opt.tree <- tree.lst$optree 
  t1 <- Sys.time()
  time.used <- as.numeric(difftime(t1, t0, units = "secs"))
  
  # evaluate opt.tree
  tmp <- eval.honest.tree(tree = opt.tree, 
                          data.test = data.test,
                          eff.modifier = data.lst$eff.modifier, 
                          true.trt.eff = data.lst$true.trt.eff,
                          trueGroupNum = data.lst$trueGroupNum,
                          splt.point   = data.lst$splt.point,
                          noise.splt.var = data.lst$noise.splt.var,
                          type.eff.mod  = data.lst$type.eff.mod)
  
  tmp$time <- time.used
  
  tmp
}



##############################################################################################
################### 2. homogeneous: honest, modelfitBefore, est.fitBefore ####################
##############################################################################################
homo.fitBefore.honest.simu.func <- function(nn, dim = 5, unmeas.conf = F) {
  
  data.lst <- homo.DataGen(N = nn, p = dim, dgp = "cov", N.test = 1000)
  if (unmeas.conf == F) {
    data <- data.lst$data
    data.test <- data.lst$data.test
  } else {
    data <- data.lst$data %>% dplyr::select(-X2)
    data.test <- data.lst$data.test %>% dplyr::select(-X2)
  }
  
  t0 <- Sys.time()
  tree.lst <- honest.Treebuild(data = data, form.outcome = "Y ~ .", form.pscore = "A ~ .", type.outcome = "continuous", 
                               train.frac = 0.8, modelFitBefore = TRUE, est.modelFitBefore = TRUE,
                               outcome.mthd = "GBM", pscore.mthd = "GLM", cp = qchisq(0.95, 1), 
                               minsplit = 100, minsize = 50)
  
  opt.tree <- tree.lst$optree 
  t1 <- Sys.time()
  time.used <- as.numeric(difftime(t1, t0, units = "secs"))  
  
  # evaluate opt.tree
  tmp <- eval.honest.tree(tree = opt.tree, 
                          data.test = data.test,
                          eff.modifier = data.lst$eff.modifier, 
                          true.trt.eff = data.lst$true.trt.eff,
                          trueGroupNum = data.lst$trueGroupNum,
                          splt.point   = data.lst$splt.point,
                          noise.splt.var = data.lst$noise.splt.var,
                          type.eff.mod  = data.lst$type.eff.mod)
  
  tmp$time <- time.used
  
  tmp
}



##############################################################################################
################### 3. heterogenuous: honest, modelfitinparent, est.fitBefore ################
##############################################################################################
hetero.fitinparent.honest.simu.func <- function(nn, dim = 5, unmeas.conf = F) {
  
  data.lst <- hetero.DataGen(N = nn, p = dim, dgp = "cov", N.test = 1000)
  if (unmeas.conf == F) {
    data <- data.lst$data
    data.test <- data.lst$data.test
  } else {
    data <- data.lst$data %>% dplyr::select(-X2)
    data.test <- data.lst$data.test %>% dplyr::select(-X2)
  }
  
  t0 <- Sys.time()
  tree.lst <- honest.Treebuild(data = data, form.outcome = "Y ~ .", form.pscore = "A ~ .", type.outcome = "continuous", 
                               train.frac = 0.8, modelFitBefore = FALSE, est.modelFitBefore = TRUE,
                               outcome.mthd = "GBM", pscore.mthd = "GLM", cp = qchisq(0.95, 1), 
                               minsplit = 60, minsize = 30)
  
  opt.tree <- tree.lst$optree 
  t1 <- Sys.time()
  time.used <- as.numeric(difftime(t1, t0, units = "secs"))  
  
  # evaluate opt.tree
  tmp <- eval.honest.tree(tree = opt.tree, 
                          data.test = data.test,
                          eff.modifier = data.lst$eff.modifier, 
                          true.trt.eff = data.lst$true.trt.eff,
                          trueGroupNum = data.lst$trueGroupNum,
                          splt.point   = data.lst$splt.point,
                          noise.splt.var = data.lst$noise.splt.var,
                          type.eff.mod  = data.lst$type.eff.mod)
  
  tmp$time <- time.used
  
  tmp
}


##############################################################################################
##################### 4. homogenuous: honest, modelfitinparent, est.fitBefore ################
##############################################################################################
homo.fitinparent.honest.simu.func <- function(nn, dim = 5, unmeas.conf = F) {
  
  data.lst <- homo.DataGen(N = nn, p = dim, dgp = "cov", N.test = 1000)
  if (unmeas.conf == F) {
    data <- data.lst$data
    data.test <- data.lst$data.test
  } else {
    data <- data.lst$data %>% dplyr::select(-X2)
    data.test <- data.lst$data.test %>% dplyr::select(-X2)
  }
  
  t0 <- Sys.time()
  tree.lst <- honest.Treebuild(data = data, form.outcome = "Y ~ .", form.pscore = "A ~ .", type.outcome = "continuous", 
                               train.frac = 0.8, modelFitBefore = FALSE, est.modelFitBefore = TRUE,
                               outcome.mthd = "GBM", pscore.mthd = "GLM", cp = qchisq(0.95, 1), 
                               minsplit = 100, minsize = 50)
  
  opt.tree <- tree.lst$optree 
  t1 <- Sys.time()
  time.used <- as.numeric(difftime(t1, t0, units = "secs"))
  
  # evaluate opt.tree
  tmp <- eval.honest.tree(tree = opt.tree, 
                          data.test = data.test,
                          eff.modifier = data.lst$eff.modifier, 
                          true.trt.eff = data.lst$true.trt.eff,
                          trueGroupNum = data.lst$trueGroupNum,
                          splt.point   = data.lst$splt.point,
                          noise.splt.var = data.lst$noise.splt.var,
                          type.eff.mod  = data.lst$type.eff.mod)
  
  tmp$time <- time.used
  
  tmp
}

 # organize the results
Output <- function(origin.res, hetero = TRUE){
  num.corr.tree <- sum(origin.res[,"corr.tree"])
  if(hetero){
    output <- origin.res %>% as.data.frame() %>%
      summarize(corr.tree.prop = mean(corr.tree, na.rm = TRUE),
                num.leaf = mean(num.leaf, na.rm = TRUE),
                num.noise = mean(num.noise.splt, na.rm = TRUE),
                fir.splt.acc = mean(fir.splt.ind, na.rm = TRUE),
                pps = mean(pps, na.rm = TRUE),
                abs.split.err = mean(abs.err.sp, na.rm = TRUE),
                mse = mean(mse, na.rm = TRUE), 
                cvr_1 = sum(cvr1)/num.corr.tree,
                cvr_2 = sum(cvr2)/num.corr.tree,
                cvr_sim = sum(sim_cvr)/num.corr.tree,
                widthCI = mean(widthCI, na.rm = TRUE),
                time = mean(time, na.rm = TRUE))
  } else {
    output <- origin.res %>% as.data.frame() %>%
      summarize(corr.tree.prop = mean(corr.tree, na.rm = TRUE),
                num.leaf = mean(num.leaf, na.rm = TRUE),
                num.noise = mean(num.noise.splt, na.rm = TRUE),
                fir.splt.acc = mean(fir.splt.ind, na.rm = TRUE),
                pps = mean(pps, na.rm = TRUE),
                abs.split.err = mean(abs.err.sp, na.rm = TRUE),
                mse = mean(mse, na.rm = TRUE), 
                cvr = sum(cvr)/num.corr.tree,
                widthCI = mean(widthCI, na.rm = TRUE),
                time = mean(time, na.rm = TRUE))
  }
  return(output)
}

# foreach
clnum <- detectCores()
cl <- makeCluster(getOption("cl.cores", clnum - 1))
registerDoParallel(cl)
pack <- c("caret", "devtools",  "dplyr", "plm", "gbm", "MASS", "rpart", "randomForestSRC")

GCT.hetero.res.fitBefore.honest   <- foreach(1:sim.reps, .combine = rbind, .packages = pack, .export = ls()) %dopar% hetero.fitBefore.honest.simu.func(nn = 2000)
GCT.homo.res.fitinparent.honest   <- foreach(1:sim.reps, .combine = rbind, .packages = pack, .export = ls()) %dopar% homo.fitinparent.honest.simu.func(nn = 2000)
GCT.hetero.res.fitBefore.honest   <- foreach(1:sim.reps, .combine = rbind, .packages = pack, .export = ls()) %dopar% hetero.fitBefore.honest.simu.func(nn = 2000)
GCT.homo.res.fitinparent.honest   <- foreach(1:sim.reps, .combine = rbind, .packages = pack, .export = ls()) %dopar% homo.fitinparent.honest.simu.func(nn = 2000)

GCT.hetero.fitBefore.honest.2000   <- Output(GCT.hetero.res.fitBefore.honest)
GCT.hetero.fitinparent.honest.2000 <- Output(GCT.hetero.res.fitinparent.honest)
GCT.homo.fitBefore.honest.2000     <- Output(GCT.homo.res.fitBefore.honest, hetero = FALSE)
GCT.homo.fitinparent.honest.2000   <- Output(GCT.homo.res.fitinparent.honest, hetero = FALSE)

stopImplicitCluster()

####################### Comparion with non-honest version GCT ################################

##############################################################################################
############# 2-1. heterogenuous: non honest, modelfitBefore, est.fitBefore  #################
##############################################################################################
hetero.fitBefore.nonhonest.simu.func <- function(nn, dim = 5, unmeas.conf = F) {
  
  data.lst <- hetero.DataGen(N = nn, p = dim, dgp = "cov", N.test = 1000)
  if (unmeas.conf == F) {
    data <- data.lst$data
    data.test <- data.lst$data.test
  } else {
    data <- data.lst$data %>% dplyr::select(-X2)
    data.test <- data.lst$data.test %>% dplyr::select(-X2)
  }
  
  t0 <- Sys.time()
  tree.lst <- Treebuild.GCTwithVar(data = data, form.outcome = "Y ~ .", form.pscore = "A ~ .", type.outcome = "continuous", 
                                   train.frac = 0.8, modelFitBefore = TRUE, outcome.mthd = "GLM", pscore.mthd = "GLM", 
                                   cp = qchisq(0.95, 1), minsplit = 30, minsize = 10)
  
  opt.tree <- tree.lst$optree 
  t1 <- Sys.time()
  time.used <- as.numeric(difftime(t1, t0, units = "secs"))
  
  # evaluate opt.tree
  tmp <- eval.honest.tree(tree = opt.tree, 
                          data.test = data.test,
                          eff.modifier = data.lst$eff.modifier, 
                          true.trt.eff = data.lst$true.trt.eff,
                          trueGroupNum = data.lst$trueGroupNum,
                          splt.point   = data.lst$splt.point,
                          noise.splt.var = data.lst$noise.splt.var,
                          type.eff.mod  = data.lst$type.eff.mod)
  
  tmp$time <- time.used
  
  tmp
}

##############################################################################################
################ 2-2. homogeneous: nonhonest, modelfitBefore, est.fitBefore ##################
##############################################################################################
homo.fitBefore.nonhonest.simu.func <- function(nn, dim = 5, unmeas.conf = F) {
  
  data.lst <- homo.DataGen(N = nn, p = dim, dgp = "cov", N.test = 1000)
  if (unmeas.conf == F) {
    data <- data.lst$data
    data.test <- data.lst$data.test
  } else {
    data <- data.lst$data %>% dplyr::select(-X2)
    data.test <- data.lst$data.test %>% dplyr::select(-X2)
  }
  
  t0 <- Sys.time()
  tree.lst <- Treebuild.GCTwithVar(data = data, form.outcome = "Y ~ .", form.pscore = "A ~ .", type.outcome = "continuous", 
                                   train.frac = 0.8, modelFitBefore = TRUE,  outcome.mthd = "GLM", pscore.mthd = "GLM", 
                                   cp = qchisq(0.95, 1), minsplit = 100, minsize = 50)
  
  opt.tree <- tree.lst$optree 
  t1 <- Sys.time()
  time.used <- as.numeric(difftime(t1, t0, units = "secs"))  
  
  # evaluate opt.tree
  tmp <- eval.honest.tree(tree = opt.tree, 
                          data.test = data.test,
                          eff.modifier = data.lst$eff.modifier, 
                          true.trt.eff = data.lst$true.trt.eff,
                          trueGroupNum = data.lst$trueGroupNum,
                          splt.point   = data.lst$splt.point,
                          noise.splt.var = data.lst$noise.splt.var,
                          type.eff.mod  = data.lst$type.eff.mod)
  
  tmp$time <- time.used
  
  tmp
}


##############################################################################################
################ 2-3. heterogenuous: nonhonest, modelfitinparent, est.fitBefore ##############
##############################################################################################
hetero.fitinparent.nonhonest.simu.func <- function(nn, dim = 5, unmeas.conf = F) {
  
  data.lst <- hetero.DataGen(N = nn, p = dim, dgp = "cov", N.test = 1000)
  if (unmeas.conf == F) {
    data <- data.lst$data
    data.test <- data.lst$data.test
  } else {
    data <- data.lst$data %>% dplyr::select(-X2)
    data.test <- data.lst$data.test %>% dplyr::select(-X2)
  }
  
  t0 <- Sys.time()
  tree.lst <- Treebuild.GCTwithVar(data = data, form.outcome = "Y ~ .", form.pscore = "A ~ .", type.outcome = "continuous", 
                                   train.frac = 0.8, modelFitBefore = FALSE, outcome.mthd = "GBM", pscore.mthd = "GLM", 
                                   cp = qchisq(0.95, 1), minsplit = 60, minsize = 30)
  
  opt.tree <- tree.lst$optree 
  t1 <- Sys.time()
  time.used <- as.numeric(difftime(t1, t0, units = "secs"))  
  
  # evaluate opt.tree
  tmp <- eval.honest.tree(tree = opt.tree, 
                          data.test = data.test,
                          eff.modifier = data.lst$eff.modifier, 
                          true.trt.eff = data.lst$true.trt.eff,
                          trueGroupNum = data.lst$trueGroupNum,
                          splt.point   = data.lst$splt.point,
                          noise.splt.var = data.lst$noise.splt.var,
                          type.eff.mod  = data.lst$type.eff.mod)
  
  tmp$time <- time.used
  
  tmp
}


##############################################################################################
################## 2-4. homogenuous: nonhonest, modelfitinparent, est.fitBefore ##############
##############################################################################################
homo.fitinparent.nonhonest.simu.func <- function(nn, dim = 5, unmeas.conf = F) {
  
  data.lst <- homo.DataGen(N = nn, p = dim, dgp = "cov", N.test = 1000)
  if (unmeas.conf == F) {
    data <- data.lst$data
    data.test <- data.lst$data.test
  } else {
    data <- data.lst$data %>% dplyr::select(-X2)
    data.test <- data.lst$data.test %>% dplyr::select(-X2)
  }
  
  t0 <- Sys.time()
  tree.lst <- Treebuild.GCTwithVar(data = data, form.outcome = "Y ~ .", form.pscore = "A ~ .", type.outcome = "continuous", 
                                   train.frac = 0.8, modelFitBefore = FALSE, outcome.mthd = "GBM", pscore.mthd = "GLM", 
                                   cp = qchisq(0.95, 1), minsplit = 100, minsize = 50)
  
  opt.tree <- tree.lst$optree 
  t1 <- Sys.time()
  time.used <- as.numeric(difftime(t1, t0, units = "secs"))
  
  # evaluate opt.tree
  tmp <- eval.honest.tree(tree = opt.tree, 
                          data.test = data.test,
                          eff.modifier = data.lst$eff.modifier, 
                          true.trt.eff = data.lst$true.trt.eff,
                          trueGroupNum = data.lst$trueGroupNum,
                          splt.point   = data.lst$splt.point,
                          noise.splt.var = data.lst$noise.splt.var,
                          type.eff.mod  = data.lst$type.eff.mod)
  
  tmp$time <- time.used
  
  tmp
}

# implement
# GCT.hetero.res.fitBefore.nonhonest   <- apply(sapply(1:sim.reps, function(iter) hetero.fitBefore.nonhonest.simu.func(nn = 2000)), 1, as.numeric)
# GCT.hetero.res.fitinparent.nonhonest <- apply(sapply(1:sim.reps, function(iter) hetero.fitinparent.nonhonest.simu.func(nn = 2000)), 1, as.numeric)
# GCT.homo.res.fitBefore.nonhonest     <- apply(sapply(1:sim.reps, function(iter) homo.fitBefore.nonhonest.simu.func(nn = 2000)), 1, as.numeric)
# GCT.homo.res.fitinparent.nonhonest   <- apply(sapply(1:sim.reps, function(iter) homo.fitinparent.nonhonest.simu.func(nn = 2000)), 1, as.numeric)
# 
# GCT.hetero.fitBefore.nonhonest.2000   <- Output(GCT.hetero.res.fitBefore.nonhonest)
# GCT.hetero.fitinparent.nonhonest.2000 <- Output(GCT.hetero.res.fitinparent.nonhonest)
# GCT.homo.fitBefore.nonhonest.2000     <- Output(GCT.homo.res.fitBefore.nonhonest, hetero = FALSE)
# GCT.homo.fitinparent.nonhonest.2000   <- Output(GCT.homo.res.fitinparent.nonhonest, hetero = FALSE)


# foreach
clnum <- detectCores()
cl <- makeCluster(getOption("cl.cores", clnum - 1))
registerDoParallel(cl)
pack <- c("caret", "devtools",  "dplyr", "plm", "gbm", "MASS", "rpart", "randomForestSRC")

GCT.hetero.res.fitBefore.nonhonest   <- foreach(1:sim.reps, .combine = rbind, .packages = pack, .export = ls()) %dopar% hetero.fitBefore.nonhonest.simu.func(nn = 1000)
GCT.hetero.res.fitinparent.nonhonest <- foreach(1:sim.reps, .combine = rbind, .packages = pack, .export = ls()) %dopar% hetero.fitinparent.nonhonest.simu.func(nn = 1000)
GCT.homo.res.fitBefore.nonhonest     <- foreach(1:sim.reps, .combine = rbind, .packages = pack, .export = ls()) %dopar% homo.fitBefore.nonhonest.simu.func(nn = 1000)
GCT.homo.res.fitinparent.nonhonest   <- foreach(1:sim.reps, .combine = rbind, .packages = pack, .export = ls()) %dopar% homo.fitinparent.nonhonest.simu.func(nn =1000)

GCT.hetero.fitBefore.nonhonest.1000.glm   <- Output(GCT.hetero.res.fitBefore.nonhonest)
GCT.hetero.fitinparent.nonhonest.1000.glm <- Output(GCT.hetero.res.fitinparent.nonhonest)
GCT.homo.fitBefore.nonhonest.1000.glm     <- Output(GCT.homo.res.fitBefore.nonhonest, hetero = FALSE)
GCT.homo.fitinparent.nonhonest.1000.glm   <- Output(GCT.homo.res.fitinparent.nonhonest, hetero = FALSE)
stopImplicitCluster()
