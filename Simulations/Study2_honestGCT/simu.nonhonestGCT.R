setwd(".\\Functions")
functions <- list.files(getwd())
functions <- paste(paste(getwd(), "\\", sep = ""), functions, sep = "")
for (f in functions) {
  source(i)
}

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
                                   train.frac = 0.8, modelFitBefore = TRUE, outcome.mthd = "RF", pscore.mthd = "GLM", 
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
                                   train.frac = 0.8, modelFitBefore = TRUE,  outcome.mthd = "RF", pscore.mthd = "GLM", 
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



# foreach
clnum <- detectCores()
cl <- makeCluster(getOption("cl.cores", clnum - 1))
registerDoParallel(cl)
pack <- c("caret", "devtools",  "dplyr", "plm", "gbm", "MASS", "rpart", "randomForestSRC")

GCT.hetero.res.fitBefore.nonhonest   <- foreach(1:sim.reps, .combine = rbind, .packages = pack, .export = ls()) %dopar% hetero.fitBefore.nonhonest.simu.func(nn = 1000)
GCT.hetero.res.fitinparent.nonhonest <- foreach(1:sim.reps, .combine = rbind, .packages = pack, .export = ls()) %dopar% hetero.fitinparent.nonhonest.simu.func(nn = 1000)
GCT.homo.res.fitBefore.nonhonest     <- foreach(1:sim.reps, .combine = rbind, .packages = pack, .export = ls()) %dopar% homo.fitBefore.nonhonest.simu.func(nn = 1000)
GCT.homo.res.fitinparent.nonhonest   <- foreach(1:sim.reps, .combine = rbind, .packages = pack, .export = ls()) %dopar% homo.fitinparent.nonhonest.simu.func(nn =1000)

GCT.hetero.fitBefore.nonhonest.1000.rf   <- Output(GCT.hetero.res.fitBefore.nonhonest)
GCT.hetero.fitinparent.nonhonest.1000.rf <- Output(GCT.hetero.res.fitinparent.nonhonest)
GCT.homo.fitBefore.nonhonest.1000.rf     <- Output(GCT.homo.res.fitBefore.nonhonest, hetero = FALSE)
GCT.homo.fitinparent.nonhonest.1000.rf   <- Output(GCT.homo.res.fitinparent.nonhonest, hetero = FALSE)

stopImplicitCluster()



# implement
# GCT.hetero.res.fitBefore.nonhonest   <- apply(sapply(1:sim.reps, function(iter) hetero.fitBefore.nonhonest.simu.func(nn = 2000)), 1, as.numeric)
# GCT.hetero.res.fitinparent.nonhonest <- apply(sapply(1:sim.reps, function(iter) hetero.fitinparent.nonhonest.simu.func(nn = 2000)), 1, as.numeric)
# GCT.homo.res.fitBefore.nonhonest     <- apply(sapply(1:sim.reps, function(iter) homo.fitBefore.nonhonest.simu.func(nn = 2000)), 1, as.numeric)
# GCT.homo.res.fitinparent.nonhonest   <- apply(sapply(1:sim.reps, function(iter) homo.fitinparent.nonhonest.simu.func(nn = 2000)), 1, as.numeric)
