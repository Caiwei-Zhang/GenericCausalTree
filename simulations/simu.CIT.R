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

################################################################################################
########################### 1. CIT-DR: heterogenuous, modelfitBefore ###########################
################################################################################################

CIT.DR.hetero.fitBefore.simu.func <- function(nn, dim = 5, unmeas.conf = F) {
  
  data.lst <- makeData.cont.eff.cont(N = nn, n.test = 1000, p = dim)
  
  if (unmeas.conf == F) {
    data <- data.lst$data.used
    test.data <- data.lst$test.data
  } else {
    data <- data.lst$data.used %>% dplyr::select(-X2)
    test.data <- data.lst$test.data %>% dplyr::select(-X2)
  }
  
  # val.sample: used in EstIpw.CvMethod1, the order of the columns must be A, Y, X
  # data <- data %>% dplyr::select("A", "Y", everything())
  data.train  <- data[1:800, ]
  data.valid  <- data[801:1000, ]  

  t0 <- Sys.time()
  prune.CIT.dr <- create.sequence(data.used         = data.train,
                                  est.used          = "DR",
                                  type.var          = "cont",
                                  propsc.mod.out    = T,
                                  propsc.mthd       = "GLM",
                                  propsc.form.true  = NULL,
                                  propsc.mod.insplt = F,
                                  adj.mod.out       = T, 
                                  adj.mthd          = "GBM", 
                                  adj.form.true     = NULL, 
                                  adj.mod.insplt    = F, 
                                  num.truc.obs      = 30,
                                  min.node          = 20)
  
  final.CIT.dr <- EstDr.CvMethod1(data.used         = data.train,
                                  tree.list         = prune.CIT.dr$tree.list, 
                                  lambda.used       = qchisq(0.95, 1),
                                  val.sample        = data.valid,
                                  type.var          = "cont",
                                  propsc.mod.out    = T,
                                  propsc.mthd       = "GLM",
                                  propsc.form.true  = NULL,
                                  propsc.mod.insplt = F,
                                  adj.mod.out       = T,
                                  adj.mthd          = "GBM",
                                  adj.form.true     = NULL,  
                                  adj.mod.insplt    = F,
                                  min.obs.mod       = 20)
  t1 <- Sys.time() 
  
  eval.CIT.dr <- eval.measures.eff(final.tree   = final.CIT.dr[[1]],
                                   test.data    = test.data,
                                   true.trt.eff = data.lst$true.trt.eff,
                                   noise.var    = data.lst$noise.splt.var,
                                   corr.split   = data.lst$eff.modifier,
                                   type.eff.mod  = data.lst$type.eff.mod,
                                   where.split  = data.lst$where.split,
                                   dir.split    = data.lst$dir.split)

  # tmp <- eval.tree(tree = final.CIT.dr[[1]],
  #                  data.test = data.lst$test.data,
  #                  true.trt.eff = data.lst$true.trt.eff,
  #                  noise.splt.var = data.lst$noise.splt.var,
  #                  eff.modifier = data.lst$eff.modifier,
  #                  trueGroupNum = data.lst$trueGroupNum,
  #                  splt.point   = data.lst$splt.point, 
  #                  type.eff.mod  = data.lst$type.eff.mod)

  eval.CIT.dr$time <- as.numeric(difftime(t1, t0, units = "secs"))
  
  eval.CIT.dr
  
}


################################################################################################
########################### 2. CIT-DR: homogeneous, modelfitBefore ############################
################################################################################################

CIT.DR.homo.fitBefore.simu.func <- function(nn, dim = 5, unmeas.conf = F) {
  
  data.lst <- makeData.cont.noeff.cont(N = nn, n.test = 1000, p = dim)
  
  if (unmeas.conf == F) {
    data <- data.lst$data.used
    test.data <- data.lst$test.data
  } else {
    data <- data.lst$data.used %>% dplyr::select(-X2)
    test.data <- data.lst$test.data %>% dplyr::select(-X2)
  }
  
  # val.sample: used in EstIpw.CvMethod1, the order of the columns must be A, Y, X
  # data <- data %>% dplyr::select("A", "Y", everything())
  data.train  <- data[1:800, ]
  data.valid  <- data[801:1000, ]  
  
  t0 <- Sys.time()
  prune.CIT.dr <- create.sequence(data.used         = data.train,
                                  est.used          = "DR",
                                  type.var          = "cont",
                                  propsc.mod.out    = T,
                                  propsc.mthd       = "GLM",
                                  propsc.form.true  = NULL,
                                  propsc.mod.insplt = F,
                                  adj.mod.out       = T, 
                                  adj.mthd          = "GBM", 
                                  adj.form.true     = NULL, 
                                  adj.mod.insplt    = F, 
                                  num.truc.obs      = 30,
                                  min.node          = 20)
  
  final.CIT.dr <- EstDr.CvMethod1(data.used         = data.train,
                                  tree.list         = prune.CIT.dr$tree.list, 
                                  lambda.used       = qchisq(0.95, 1),
                                  val.sample        = data.valid,
                                  type.var          = "cont",
                                  propsc.mod.out    = T,
                                  propsc.mthd       = "GLM",
                                  propsc.form.true  = NULL,
                                  propsc.mod.insplt = F,
                                  adj.mod.out       = T,
                                  adj.mthd          = "GBM",
                                  adj.form.true     = NULL, 
                                  adj.mod.insplt    = F,
                                  min.obs.mod       = 20)
  t1 <- Sys.time() 
  
  eval.CIT.dr <- eval.measures.eff(final.tree   = final.CIT.dr[[1]],
                                   test.data    = test.data,
                                   true.trt.eff = data.lst$true.trt.eff,
                                   noise.var    = data.lst$noise.splt.var,
                                   corr.split   = data.lst$eff.modifier,
                                   type.eff.mod  = data.lst$type.eff.mod,
                                   where.split  = data.lst$where.split,
                                   dir.split    = data.lst$dir.split)
  
  eval.CIT.dr$time <- as.numeric(difftime(t1, t0, units = "secs"))
  
  eval.CIT.dr
  
}



################################################################################################
########################### 3. CIT-DR: heterogenuous, modelfitinparent #########################
################################################################################################

CIT.DR.hetero.fitinparent.simu.func <- function(nn, dim = 5, unmeas.conf = F) {
  
  data.lst <- makeData.cont.eff.cont(N = nn, n.test = 1000, p = dim)
  
  if (unmeas.conf == F) {
    data <- data.lst$data.used
    test.data <- data.lst$test.data
  } else {
    data <- data.lst$data.used %>% dplyr::select(-X2)
    test.data <- data.lst$test.data %>% dplyr::select(-X2)
  }
  # val.sample: used in EstIpw.CvMethod1, the order of the columns must be A, Y, X
  # data <- data %>% dplyr::select("A", "Y", everything())
  
  data.train  <- data[1:800, ]
  data.valid  <- data[801:1000, ]  
  
  
  t0 <- Sys.time()
  prune.CIT.dr <- create.sequence(data.used         = data.train,
                                  est.used          = "DR",
                                  type.var          = "cont",
                                  propsc.mod.out    = F,
                                  propsc.mthd       = "GLM",
                                  propsc.form.true  = NULL,
                                  propsc.mod.insplt = F,
                                  adj.mod.out       = F, 
                                  adj.mthd          = "GBM", 
                                  adj.form.true     = NULL, 
                                  adj.mod.insplt    = F, 
                                  num.truc.obs      = 30,
                                  min.node          = 30)
  
  final.CIT.dr <- EstDr.CvMethod1(data.used         = data.train,
                                  tree.list         = prune.CIT.dr$tree.list, 
                                  lambda.used       = qchisq(0.95, 1),
                                  val.sample        = data.valid,
                                  type.var          = "cont",
                                  propsc.mod.out    = F,
                                  propsc.mthd       = "GLM",
                                  propsc.form.true  = NULL,
                                  propsc.mod.insplt = F,
                                  adj.mod.out       = F,
                                  adj.mthd          = "GBM",
                                  adj.form.true     = NULL,  
                                  adj.mod.insplt    = F,
                                  min.obs.mod       = 30)
  t1 <- Sys.time() 
  
  eval.CIT.dr <- eval.measures.eff(final.tree   = final.CIT.dr[[1]],
                                   test.data    = test.data,
                                   true.trt.eff = data.lst$true.trt.eff,
                                   noise.var    = data.lst$noise.splt.var,
                                   corr.split   = data.lst$eff.modifier,
                                   type.eff.mod  = data.lst$type.eff.mod,
                                   where.split  = data.lst$where.split,
                                   dir.split    = data.lst$dir.split)
  eval.CIT.dr$time <- as.numeric(difftime(t1, t0, units = "secs"))
  
  eval.CIT.dr
  
}


################################################################################################
########################## 4. CIT-DR: homogeneous, modelfitinparent ############################
################################################################################################

CIT.DR.homo.fitinparent.simu.func <- function(nn, dim = 5, unmeas.conf = F) {
  
  data.lst <- makeData.cont.noeff.cont(N = nn, n.test = 1000, p = dim)
  
  if (unmeas.conf == F) {
    data <- data.lst$data.used
    test.data <- data.lst$test.data
  } else {
    data <- data.lst$data.used %>% dplyr::select(-X2)
    test.data <- data.lst$test.data %>% dplyr::select(-X2)
  }
  
  # val.sample: used in EstIpw.CvMethod1, the order of the columns must be A, Y, X
  # data <- data %>% dplyr::select("A", "Y", everything())
  
  data.train  <- data[1:800, ]
  data.valid  <- data[801:1000, ]  
  
  
  t0 <- Sys.time()
  prune.CIT.dr <- create.sequence(data.used         = data.train,
                                  est.used          = "DR",
                                  type.var          = "cont",
                                  propsc.mod.out    = F,
                                  propsc.mthd       = "GLM",
                                  propsc.form.true  = NULL,
                                  propsc.mod.insplt = F,
                                  adj.mod.out       = F, 
                                  adj.mthd          = "GBM", 
                                  adj.form.true     = NULL, 
                                  adj.mod.insplt    = F, 
                                  num.truc.obs      = 30,
                                  min.node          = 30)
  
  final.CIT.dr <- EstDr.CvMethod1(data.used         = data.train,
                                  tree.list         = prune.CIT.dr$tree.list, 
                                  lambda.used       = qchisq(0.95, 1),
                                  val.sample        = data.valid,
                                  type.var          = "cont",
                                  propsc.mod.out    = F,
                                  propsc.mthd       = "GLM",
                                  propsc.form.true  = NULL,
                                  propsc.mod.insplt = F,
                                  adj.mod.out       = F,
                                  adj.mthd          = "GBM",
                                  adj.form.true     = NULL,  
                                  adj.mod.insplt    = F,
                                  min.obs.mod       = 30)
  t1 <- Sys.time() 
  
  eval.CIT.dr <- eval.measures.eff(final.tree   = final.CIT.dr[[1]],
                                   test.data    = test.data,
                                   true.trt.eff = data.lst$true.trt.eff,
                                   noise.var    = data.lst$noise.splt.var,
                                   corr.split   = data.lst$eff.modifier,
                                   type.eff.mod = data.lst$type.eff.mod,
                                   where.split  = data.lst$where.split,
                                   dir.split    = data.lst$dir.split)
  
  eval.CIT.dr$time <- as.numeric(difftime(t1, t0, units = "secs"))
  
  eval.CIT.dr
  
}


set.seed(999)
CIT.DR.hetero.res.fitBefore   <- apply(sapply(1:sim.reps, function(iter) CIT.DR.hetero.fitBefore.simu.func(nn = 1000)), 1, as.numeric)
CIT.DR.homo.res.fitBefore     <- apply(sapply(1:sim.reps, function(iter) CIT.DR.homo.fitBefore.simu.func(nn = 1000)), 1, as.numeric)
CIT.DR.hetero.res.fitinparent <- apply(sapply(1:sim.reps, function(iter) CIT.DR.hetero.fitinparent.simu.func(nn = 1000)), 1, as.numeric)
CIT.DR.homo.res.fitinparent   <- apply(sapply(1:sim.reps, function(iter) CIT.DR.homo.fitinparent.simu.func(nn = 1000)), 1, as.numeric)

CIT.DR.hetero.fitBefore   <- colMeans(CIT.DR.hetero.res.fitBefore, na.rm = TRUE)
CIT.DR.homo.fitBefore     <- colMeans(CIT.DR.homo.res.fitBefore, na.rm = TRUE)
CIT.DR.hetero.fitinparent <- colMeans(CIT.DR.hetero.res.fitinparent, na.rm = TRUE)
CIT.DR.homo.fitinparent   <- colMeans(CIT.DR.homo.res.fitinparent, na.rm = TRUE)



# foreach (seed can not be consistent)
# clnum <- detectCores()
# cl <- makeCluster(getOption("cl.cores", clnum - 1))
# registerDoParallel(cl)
# pack <- c("caret", "devtools", "causalTree", "causalForest", "dplyr", "grf",
#           "gbm", "MASS", "rpart", "randomForestSRC", "rattle", "rpart.plot")
# CIT.DR.hetero.res.fitBefore <- foreach(1:sim.reps, .combine = rbind, 
#                                 .export = ls(), .packages = pack) %dopar% CIT.DR.hetero.fitBefore.simu.func(nn = 1000)
# 
# 
# CIT.DR.homo.res.fitBefore <- foreach(1:sim.reps, .combine = rbind, 
#                               .export = ls(), .packages = pack) %dopar% CIT.DR.homo.fitBefore.simu.func(nn = 1000)
# 
# 
# CIT.DR.hetero.res.fitinparent <- foreach(1:sim.reps, .combine = rbind, 
#                                   .export = ls(), .packages = pack) %dopar% CIT.DR.hetero.fitinparent.simu.func(nn = 1000)
# 
# 
# CIT.DR.homo.res.fitinparent <- foreach(1:sim.reps, .combine = rbind, 
#                                 .export = ls(), .packages = pack) %dopar% CIT.DR.homo.fitinparent.simu.func(nn = 1000)
# 
# 
# stopImplicitCluster()
# 
# CIT.DR.hetero.fitBefore <- colMeans(CIT.DR.hetero.res.fitBefore, na.rm = TRUE)
# CIT.DR.homo.fitBefore   <- colMeans(CIT.DR.homo.res.fitBefore, na.rm = TRUE)
# CIT.DR.hetero.fitinparent <- colMeans(CIT.DR.hetero.res.fitinparent, na.rm = TRUE)
# CIT.DR.homo.fitinparent   <- colMeans(CIT.DR.homo.res.fitinparent, na.rm = TRUE)


save.image(file = "E:/Project/groupTest/MethodCompare/simulation/simu.CIT.RData")

save(CIT.DR.hetero.res.fitBefore, CIT.DR.homo.res.fitBefore, 
     CIT.DR.hetero.res.fitinparent, CIT.DR.homo.res.fitinparent,
     CIT.DR.hetero.fitBefore, CIT.DR.homo.fitBefore, 
     CIT.DR.hetero.fitinparent, CIT.DR.homo.fitinparent,
     file = "E:/Project/groupTest/MethodCompare/simulation/res.CIT.RData")

