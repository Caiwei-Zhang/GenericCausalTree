load("seed.RData")
setwd("./Functions")
file <- list.files(getwd())
filePath <- paste(paste(getwd(), "/", sep = ""), file, sep = "")
filePath <- grep("all", filePath, invert = TRUE, value = TRUE)
for (f in filePath) {
  print(f)
  source(f)
}


##############################################################################################
############# 2-1. heterogenuous: non honest, modelfitBefore, est.fitBefore  #################
##############################################################################################
hetero.fitBefore.nonhonest.simu.func <- function(nn, dim = 5, unmeas.conf = F, outcome.mthd = "GBM", seed = NULL) {
  
  if(!is.null(seed)){
    set.seed(seed)
  }
  
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
                                   train.frac = 0.8, modelFitBefore = TRUE, outcome.mthd = outcome.mthd, pscore.mthd = "GLM", 
                                   cp = qchisq(0.95, 1), minsplit = 30, minsize = 10) #minsplit = nn/100, minsize = nn/200
  
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
homo.fitBefore.nonhonest.simu.func <- function(nn, dim = 5, unmeas.conf = F, outcome.mthd = "GBM", seed = NULL) {
  
  if(!is.null(seed)){
    set.seed(seed)
  }
  
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
                                   train.frac = 0.8, modelFitBefore = TRUE,  outcome.mthd = outcome.mthd, pscore.mthd = "GLM", 
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
hetero.fitinparent.nonhonest.simu.func <- function(nn, dim = 5, unmeas.conf = F, outcome.mthd = "GBM", seed = NULL) {
  
  if(!is.null(seed)){
    set.seed(seed)
  }
  
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
                                   train.frac = 0.8, modelFitBefore = FALSE, outcome.mthd = outcome.mthd, pscore.mthd = "GLM", 
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
homo.fitinparent.nonhonest.simu.func <- function(nn, dim = 5, unmeas.conf = F, outcome.mthd = "GBM", seed = NULL) {
  
  if(!is.null(seed)){
    set.seed(seed)
  }
  
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
                                   train.frac = 0.8, modelFitBefore = FALSE, outcome.mthd = outcome.mthd, pscore.mthd = "GLM", 
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
sim.reps <- 1000
clnum <- detectCores()
cl <- makeCluster(getOption("cl.cores", clnum - 4))
registerDoParallel(cl)
pack <- c("caret", "devtools", "causalTree", "dplyr",  "gbm", "MASS", "rpart", "randomForestSRC", "xgboost")


# GBM
GCT.hetero.res.fitBefore.nonhonest   <- foreach(i = 1:sim.reps, .combine = rbind, .packages = pack, .export = ls()) %dopar% 
  hetero.fitBefore.nonhonest.simu.func(nn = 1000, seed = seed[i])
GCT.hetero.res.fitinparent.nonhonest <- foreach(i = 1:sim.reps, .combine = rbind, .packages = pack, .export = ls()) %dopar% 
  hetero.fitinparent.nonhonest.simu.func(nn = 1000, seed = seed[i])
GCT.homo.res.fitinparent.nonhonest   <- foreach(i = 1:sim.reps, .combine = rbind, .packages = pack, .export = ls()) %dopar% 
  homo.fitinparent.nonhonest.simu.func(nn = 1000, seed = seed[i + sim.reps])
GCT.homo.res.fitBefore.nonhonest     <- foreach(i = 1:sim.reps, .combine = rbind, .packages = pack, .export = ls()) %dopar% 
  homo.fitBefore.nonhonest.simu.func(nn = 1000,  seed = seed[i + sim.reps])

GCT.hetero.fitBefore.nonhonest.1000   <- Output(GCT.hetero.res.fitBefore.nonhonest)
GCT.hetero.fitinparent.nonhonest.1000 <- Output(GCT.hetero.res.fitinparent.nonhonest)
GCT.homo.fitBefore.nonhonest.1000     <- Output(GCT.homo.res.fitBefore.nonhonest, hetero = FALSE)
GCT.homo.fitinparent.nonhonest.1000   <- Output(GCT.homo.res.fitinparent.nonhonest, hetero = FALSE)


save(GCT.hetero.res.fitBefore.nonhonest,  
     GCT.hetero.res.fitinparent.nonhonest,
     GCT.homo.res.fitBefore.nonhonest,    
     GCT.homo.res.fitinparent.nonhonest, 
     GCT.hetero.fitBefore.nonhonest.1000,   
     GCT.hetero.fitinparent.nonhonest.1000, 
     GCT.homo.fitBefore.nonhonest.1000,     
     GCT.homo.fitinparent.nonhonest.1000,
     file = "./Simulations/Study1_MethodComparison/simu.GCT/seed_res/GCT.TreebuildwithVar.nonhonest.1000.RData")

# RF & GLM
GCT.hetero.res.fitBefore.nonhonest.rf   <- foreach(i = 1:sim.reps, .combine = rbind, .packages = pack, .export = ls()) %dopar% 
  hetero.fitBefore.nonhonest.simu.func(nn = 1000, outcome.mthd = "RF", seed = seed[i])
GCT.homo.res.fitBefore.nonhonest.rf   <- foreach(i = 1:sim.reps, .combine = rbind, .packages = pack, .export = ls()) %dopar% 
  homo.fitBefore.nonhonest.simu.func(nn = 1000, outcome.mthd = "RF", seed = seed[i + sim.reps])

GCT.hetero.res.fitBefore.nonhonest.glm   <- foreach(i = 1:sim.reps, .combine = rbind, .packages = pack, .export = ls()) %dopar% 
  hetero.fitBefore.nonhonest.simu.func(nn = 1000, outcome.mthd = "GLM", seed = seed[i])
GCT.homo.res.fitBefore.nonhonest.glm   <- foreach(i = 1:sim.reps, .combine = rbind, .packages = pack, .export = ls()) %dopar% 
  homo.fitBefore.nonhonest.simu.func(nn = 1000, outcome.mthd = "GLM", seed = seed[i + sim.reps])

GCT.hetero.fitBefore.nonhonest.1000.rf   <- Output(GCT.hetero.res.fitBefore.nonhonest.rf)
GCT.homo.fitBefore.nonhonest.1000.rf     <- Output(GCT.homo.res.fitBefore.nonhonest.rf, hetero = FALSE)
GCT.hetero.fitBefore.nonhonest.1000.glm   <- Output(GCT.hetero.res.fitBefore.nonhonest.glm)
GCT.homo.fitBefore.nonhonest.1000.glm     <- Output(GCT.homo.res.fitBefore.nonhonest.glm, hetero = FALSE)

load("./Simulations/Study1_MethodComparison/simu.GCT/seed_res/GCT.TreebuildwithVar.nonhonest.1000.RData")
GCT.hetero.fitBefore.nonhonest.1000.gbm <- GCT.hetero.fitBefore.nonhonest.1000
GCT.homo.fitBefore.nonhonest.1000.gbm   <- GCT.homo.fitBefore.nonhonest.1000
GCT.hetero.res.fitBefore.nonhonest.gbm  <- GCT.hetero.res.fitBefore.nonhonest
GCT.homo.res.fitBefore.nonhonest.gbm    <- GCT.homo.res.fitBefore.nonhonest

save(GCT.hetero.fitBefore.nonhonest.1000.gbm,
     GCT.homo.fitBefore.nonhonest.1000.gbm,
     GCT.hetero.fitBefore.nonhonest.1000.rf,  
     GCT.homo.fitBefore.nonhonest.1000.rf,    
     GCT.hetero.fitBefore.nonhonest.1000.glm, 
     GCT.homo.fitBefore.nonhonest.1000.glm, 
     GCT.hetero.res.fitBefore.nonhonest.gbm,
     GCT.homo.res.fitBefore.nonhonest.gbm,  
     GCT.hetero.res.fitBefore.nonhonest.rf, 
     GCT.homo.res.fitBefore.nonhonest.rf, 
     GCT.hetero.res.fitBefore.nonhonest.glm,
     GCT.homo.res.fitBefore.nonhonest.glm,
     file = "./Simulations/Appendix_ML/seed_res/GCT.nonhonest.1000.ML.RData")

stopImplicitCluster()

memory.size()
gc()



