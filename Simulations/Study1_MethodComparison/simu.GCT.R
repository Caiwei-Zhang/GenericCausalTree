setwd("E:\\3_GCT\\Functions")
source("library.R")
source("integr.model.fit.R")
source("gbm.predict.R")
source("SearchChilds.R")
source("ulist.R")
source("dgp.R")
source("Treebuild.R")
source("Treebuild.honest.R")
source("eval.tree.R")
source("eval.honest.tree.R")

sim.reps <- 1000
set.seed(999)
## non-honest
################################################################################################
########################## 1. heterogenuous: nonhonest, modelfitBefore #########################
################################################################################################
hetero.fitBefore.simu.func <- function(nn, dim = 5, unmeas.conf = F) {
  
  data.lst <- hetero.DataGen(N = nn, p = dim, dgp = "cov", N.test = 1000)
  if (unmeas.conf == F) {
    data <- data.lst$data
    data.test <- data.lst$data.test
  } else {
    data <- data.lst$data %>% dplyr::select(-X2)
    data.test <- data.lst$data.test %>% dplyr::select(-X2)
  }
  
  t0 <- Sys.time()
  tree.lst <- Treebuild(data = data, form.outcome = "Y ~ .", form.pscore = "A ~ .", 
                        type.outcome = "continuous", train.frac = 0.8, 
                        modelFitBefore = TRUE, outcome.mthd = "GBM", pscore.mthd = "GLM", 
                        cp = qchisq(0.95, 1), minsplit = 30, minsize = 10)
  opt.tree <- tree.lst$optree 
  t1 <- Sys.time()
  time.used <- as.numeric(difftime(t1, t0, units = "secs"))
  
  # evaluate opt.tree
  tmp <- eval.tree(tree = opt.tree,
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



###################################################################################################
########################## 2. homogenuous: nonhonest, modelfitBefore ##############################
###################################################################################################

homo.fitBefore.simu.func <- function(nn, dim = 5, unmeas.conf = F) {
  
  data.lst <- homo.DataGen(N = nn, p = dim, dgp = "cov", N.test = 1000)
  if (unmeas.conf == F) {
    data <- data.lst$data
    data.test <- data.lst$data.test
  } else {
    data <- data.lst$data %>% dplyr::select(-X2)
    data.test <- data.lst$data.test %>% dplyr::select(-X2)
  }
  
  t0 <- Sys.time()
  tree.lst <- Treebuild(data = data, form.outcome = "Y ~ .", form.pscore = "A ~ .", type.outcome = "continuous", 
                        train.frac = 0.8, modelFitBefore = TRUE, outcome.mthd = "GBM", pscore.mthd = "GLM", 
                        cp = qchisq(0.95, 1), minsplit = 30, minsize = 10)
  
  opt.tree <- tree.lst$optree 
  t1 <- Sys.time()
  time.used <- as.numeric(difftime(t1, t0, units = "secs")) 
  
  # evaluate opt.tree
  tmp <- eval.tree(tree = opt.tree, 
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



################################################################################################
######################## 3. heterogeneous: nonhonest, modelfitinparent #########################
################################################################################################
hetero.fitinparent.simu.func <- function(nn, dim = 5, unmeas.conf = F) {
  
  data.lst <- hetero.DataGen(N = nn, p = dim, dgp = "cov", N.test = 1000)
  if (unmeas.conf == F) {
    data <- data.lst$data
    data.test <- data.lst$data.test
  } else {
    data <- data.lst$data %>% dplyr::select(-X2)
    data.test <- data.lst$data.test %>% dplyr::select(-X2)
  }
  
  t0 <- Sys.time()
  tree.lst <- Treebuild(data = data, form.outcome = "Y ~ .", form.pscore = "A ~ .", type.outcome = "continuous", 
                        train.frac = 0.8, modelFitBefore = FALSE, outcome.mthd = "GBM", pscore.mthd = "GLM", 
                        cp = qchisq(0.95, 1), minsplit = 30, minsize = 20)
  
  opt.tree <- tree.lst$optree 
  t1 <- Sys.time()
  time.used <- as.numeric(difftime(t1, t0, units = "secs"))
  
  # evaluate opt.tree
  tmp <- eval.tree(tree = opt.tree, 
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


################################################################################################
######################## 4. homogeneous: nonhonest, modelfitinparent ###########################
################################################################################################
homo.fitinparent.simu.func <- function(nn, dim = 5, unmeas.conf = F) {
  
  data.lst <- homo.DataGen(N = nn, p = dim, dgp = "cov", N.test = 1000)
  if (unmeas.conf == F) {
    data <- data.lst$data
    data.test <- data.lst$data.test
  } else {
    data <- data.lst$data %>% dplyr::select(-X2)
    data.test <- data.lst$data.test %>% dplyr::select(-X2)
  }
  
  t0 <- Sys.time()
  tree.lst <- Treebuild(data = data, form.outcome = "Y ~ .", form.pscore = "A ~ .", type.outcome = "continuous", 
                        train.frac = 0.8, modelFitBefore = FALSE, outcome.mthd = "GBM", pscore.mthd = "GLM", 
                        cp =qchisq(0.95, 1), minsplit = 30, minsize = 20)
  
  opt.tree <- tree.lst$optree 
  t1 <- Sys.time()
  time.used <- as.numeric(difftime(t1, t0, units = "secs"))
  
  # evaluate opt.tree
  tmp <- eval.tree(tree = opt.tree, 
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





# sapply
set.seed(999)
GCT.hetero.res.fitBefore   <- apply(sapply(1:sim.reps, function(iter) hetero.fitBefore.simu.func(nn = 1000)), 1, as.numeric)
GCT.homo.res.fitBefore     <- apply(sapply(1:sim.reps, function(iter) homo.fitBefore.simu.func(nn = 1000)), 1, as.numeric)
GCT.hetero.res.fitinparent <- apply(sapply(1:sim.reps, function(iter) hetero.fitinparent.simu.func(nn = 1000)), 1, as.numeric)
GCT.homo.res.fitinparent   <- apply(sapply(1:sim.reps, function(iter) homo.fitinparent.simu.func(nn = 1000)), 1, as.numeric)

GCT.hetero.fitBefore   <- colMeans(GCT.hetero.res.fitBefore, na.rm = TRUE)
GCT.homo.fitBefore     <- colMeans(GCT.homo.res.fitBefore, na.rm = TRUE)
GCT.hetero.fitinparent <- colMeans(GCT.hetero.res.fitinparent, na.rm = TRUE)
GCT.homo.fitinparent   <- colMeans(GCT.homo.res.fitinparent, na.rm = TRUE)

# foreach
clnum <- detectCores()
cl <- makeCluster(getOption("cl.cores", clnum - 1))
registerDoParallel(cl)
pack <- c("caret", "devtools", "causalTree", "causalForest", "dplyr", "grf",
          "gbm", "MASS", "rpart", "randomForestSRC", "rattle", "rpart.plot")
hetero.res.fitBefore <- foreach(1:sim.reps, .combine = rbind, .export = ls(), .packages = pack) %dopar% hetero.fitBefore.simu.func(nn=1000)
homo.res.fitBefore <- foreach(1:sim.reps, .combine = rbind, .export = ls(), .packages = pack) %dopar% homo.fitBefore.simu.func(nn=1000)
hetero.res.fitinparent <- foreach(1:sim.reps, .combine = rbind, .export = ls(), .packages = pack) %dopar% hetero.fitinparent.simu.func(nn=1000)
homo.res.fitinparent <- foreach(1:sim.reps, .combine = rbind, .export = ls(), .packages = pack) %dopar% homo.fitinparent.simu.func(nn=1000)

stopImplicitCluster()

# organize the results
GCT.hetero.fitBefore   <- colMeans(hetero.res.fitBefore, na.rm = TRUE)
GCT.homo.fitBefore     <- colMeans(homo.res.fitBefore, na.rm = TRUE)
GCT.hetero.fitinparent <- colMeans(hetero.res.fitinparent, na.rm = TRUE)
GCT.homo.fitinparent   <- colMeans(homo.res.fitinparent, na.rm = TRUE)

# save the results
save(hetero.res.fitBefore, homo.res.fitBefore, hetero.res.fitinparent, homo.res.fitinparent,
     GCT.hetero.fitBefore, GCT.homo.fitBefore, GCT.hetero.fitinparent, GCT.homo.fitinparent,
     file = "./res.GCT.non_honest.RData")

