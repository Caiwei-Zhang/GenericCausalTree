setwd("./CIT-master/Functions/")  #### setwd("..")
functions <- list.files(getwd())
functions <- grep("main.R", functions, invert = T, value = T)
functions <- paste(getwd(), functions, sep = "/")
for (i in functions){
  source(i)
}


################################################################################################
########################### 1. CIT-DR: heterogenuous, modelfitBefore ###########################
################################################################################################
CIT.DR.hetero.fitBefore.simu.func <- function(nn, dim = 5, unmeas.conf = F) {
  
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
  

  t0 <- Sys.time()
  # val.sample: used in EstIpw.CvMethod1, the order of the columns must be A, Y, X
  # data <- data %>% dplyr::select("A", "Y", everything())
  data.train  <- data[1:800, ]
  data.valid  <- data[801:1000, ]
  
  prune.CIT.dr <- create.sequence(data.used         = data.train,
                                  est.used          = "DR",
                                  type.var          = "cont",
                                  propsc.mod.loc    = "out",
                                  propsc.mthd       = "GLM",
                                  propsc.form.true  = NULL,
                                  adj.mod.loc       = "out", 
                                  adj.mthd          = "GBM", 
                                  adj.form.true     = NULL, 
                                  num.truc.obs      = 30,
                                  min.node          = 10)
  
  final.CIT.dr <- EstDr.CvMethod1(data.used         = data.train,
                                  tree.list         = prune.CIT.dr$tree.list, 
                                  lambda.used       = qchisq(0.95, 1),
                                  val.sample        = data.valid,
                                  type.var          = "cont",
                                  propsc.mod.loc    = "out",
                                  propsc.mthd       = "GLM",
                                  propsc.form.true  = NULL,
                                  adj.mod.loc       = "out",
                                  adj.mthd          = "GBM",
                                  adj.form.true     = NULL)
  t1 <- Sys.time() 
  
  eval.CIT.dr <- eval.measures.eff(final.tree   = final.CIT.dr[[1]],
                                   test.data    = data.test,
                                   true.trt.eff = data.lst$true.trt.eff,
                                   noise.var    = data.lst$noise.splt.var,
                                   corr.split   = data.lst$eff.modifier,
                                   #type.eff.mod  = data.lst$type.eff.mod,
                                   where.split  = data.lst$where.split,
                                   dir.split    = data.lst$dir.split)


  eval.CIT.dr$time <- as.numeric(difftime(t1, t0, units = "secs"))
  
  return(unlist(eval.CIT.dr))
  
}


################################################################################################
########################### 2. CIT-DR: homogeneous, modelfitBefore ############################
################################################################################################

CIT.DR.homo.fitBefore.simu.func <- function(nn, dim = 5, unmeas.conf = F) {
  
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
  
  
  # val.sample: used in EstIpw.CvMethod1, the order of the columns must be A, Y, X
  # data <- data %>% dplyr::select("A", "Y", everything())
  data.train  <- data[1:800, ]
  data.valid  <- data[801:1000, ]  
  
  t0 <- Sys.time()
  prune.CIT.dr <- create.sequence(data.used         = data.train,
                                  est.used          = "DR",
                                  type.var          = "cont",
                                  propsc.mod.loc    = "out",
                                  propsc.mthd       = "GLM",
                                  propsc.form.true  = NULL,
                                  adj.mod.loc       = "out", 
                                  adj.mthd          = "GBM", 
                                  adj.form.true     = NULL, 
                                  num.truc.obs      = 100,
                                  min.node          = 50)
  
  final.CIT.dr <- EstDr.CvMethod1(data.used         = data.train,
                                  tree.list         = prune.CIT.dr$tree.list, 
                                  lambda.used       = qchisq(0.95, 1),
                                  val.sample        = data.valid,
                                  type.var          = "cont",
                                  propsc.mod.loc    = "out",
                                  propsc.mthd       = "GLM",
                                  propsc.form.true  = NULL,
                                  adj.mod.loc       = "out",
                                  adj.mthd          = "GBM",
                                  adj.form.true     = NULL)
  t1 <- Sys.time() 
  
  eval.CIT.dr <- eval.measures.eff(final.tree   = final.CIT.dr[[1]],
                                   test.data    = data.test,
                                   true.trt.eff = data.lst$true.trt.eff,
                                   noise.var    = data.lst$noise.splt.var,
                                   corr.split   = data.lst$eff.modifier,
                                   #type.eff.mod  = data.lst$type.eff.mod,
                                   where.split  = data.lst$where.split,
                                   dir.split    = data.lst$dir.split)
  
  eval.CIT.dr$time <- as.numeric(difftime(t1, t0, units = "secs"))
  
  return(unlist(eval.CIT.dr))
  
}



################################################################################################
########################### 3. CIT-DR: heterogenuous, modelfitinparent #########################
################################################################################################

CIT.DR.hetero.fitinparent.simu.func <- function(nn, dim = 5, unmeas.conf = F) {
  
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
  
  # val.sample: used in EstIpw.CvMethod1, the order of the columns must be A, Y, X
  # data <- data %>% dplyr::select("A", "Y", everything())
  
  data.train  <- data[1:800, ]
  data.valid  <- data[801:1000, ]  
  
  
  t0 <- Sys.time()
  prune.CIT.dr <- create.sequence(data.used         = data.train,
                                  est.used          = "DR",
                                  type.var          = "cont",
                                  propsc.mod.loc    = "node",
                                  propsc.mthd       = "GLM",
                                  propsc.form.true  = NULL,
                                  adj.mod.loc       = "node", 
                                  adj.mthd          = "GBM", 
                                  adj.form.true     = NULL, 
                                  num.truc.obs      = 60,
                                  min.node          = 30)
  
  final.CIT.dr <- EstDr.CvMethod1(data.used         = data.train,
                                  tree.list         = prune.CIT.dr$tree.list, 
                                  lambda.used       = qchisq(0.95, 1),
                                  val.sample        = data.valid,
                                  type.var          = "cont",
                                  propsc.mod.loc    = "node",
                                  propsc.mthd       = "GLM",
                                  propsc.form.true  = NULL,
                                  adj.mod.loc       = "node",
                                  adj.mthd          = "GBM",
                                  adj.form.true     = NULL)
  t1 <- Sys.time() 
  
  eval.CIT.dr <- eval.measures.eff(final.tree   = final.CIT.dr[[1]],
                                   test.data    = data.test,
                                   true.trt.eff = data.lst$true.trt.eff,
                                   noise.var    = data.lst$noise.splt.var,
                                   corr.split   = data.lst$eff.modifier,
                                   #type.eff.mod  = data.lst$type.eff.mod,
                                   where.split  = data.lst$where.split,
                                   dir.split    = data.lst$dir.split)
  eval.CIT.dr$time <- as.numeric(difftime(t1, t0, units = "secs"))
  
  return(unlist(eval.CIT.dr))
  
}


################################################################################################
########################## 4. CIT-DR: homogeneous, modelfitinparent ############################
################################################################################################

CIT.DR.homo.fitinparent.simu.func <- function(nn, dim = 5, unmeas.conf = F) {
  
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
  
  # val.sample: used in EstIpw.CvMethod1, the order of the columns must be A, Y, X
  # data <- data %>% dplyr::select("A", "Y", everything())
  
  data.train  <- data[1:800, ]
  data.valid  <- data[801:1000, ]  
  
  
  t0 <- Sys.time()
  prune.CIT.dr <- create.sequence(data.used         = data.train,
                                  est.used          = "DR",
                                  type.var          = "cont",
                                  propsc.mod.loc    = "node",
                                  propsc.mthd       = "GLM",
                                  propsc.form.true  = NULL,
                                  adj.mod.loc       = "node", 
                                  adj.mthd          = "GBM", 
                                  adj.form.true     = NULL, 
                                  num.truc.obs      = 100,
                                  min.node          = 50)
  
  final.CIT.dr <- EstDr.CvMethod1(data.used         = data.train,
                                  tree.list         = prune.CIT.dr$tree.list, 
                                  lambda.used       = qchisq(0.95, 1),
                                  val.sample        = data.valid,
                                  type.var          = "cont",
                                  propsc.mod.loc    = "node",
                                  propsc.mthd       = "GLM",
                                  propsc.form.true  = NULL,
                                  adj.mod.loc       = "node",
                                  adj.mthd          = "GBM",
                                  adj.form.true     = NULL)
  t1 <- Sys.time() 
  
  eval.CIT.dr <- eval.measures.eff(final.tree   = final.CIT.dr[[1]],
                                   test.data    = data.test,
                                   true.trt.eff = data.lst$true.trt.eff,
                                   noise.var    = data.lst$noise.splt.var,
                                   corr.split   = data.lst$eff.modifier,
                                   #type.eff.mod = data.lst$type.eff.mod,
                                   where.split  = data.lst$where.split,
                                   dir.split    = data.lst$dir.split)
  
  eval.CIT.dr$time <- as.numeric(difftime(t1, t0, units = "secs"))
  
  return(unlist(eval.CIT.dr))
  
}

sim.reps <- 1000
set.seed(999)

# foreach
clnum <- detectCores()
cl <- makeCluster(getOption("cl.cores", clnum - 4))
registerDoParallel()
pack <- c("caret", "devtools", "gbm", "dplyr", "rpart", "MASS", "randomForestSRC", "plm")

CIT.DR.hetero.res.fitBefore   <- foreach(1:sim.reps, .combine = rbind, .export = ls(), .packages = pack) %dopar% CIT.DR.hetero.fitBefore.simu.func(nn = 1000)
CIT.DR.homo.res.fitBefore     <- foreach(1:sim.reps, .combine = rbind, .export = ls(), .packages = pack) %dopar% CIT.DR.homo.fitBefore.simu.func(nn = 1000)
CIT.DR.hetero.res.fitinparent <- foreach(1:sim.reps, .combine = rbind, .export = ls(), .packages = pack) %dopar% CIT.DR.hetero.fitinparent.simu.func(nn = 1000)
CIT.DR.homo.res.fitinparent   <- foreach(1:sim.reps, .combine = rbind, .export = ls(), .packages = pack) %dopar% CIT.DR.homo.fitinparent.simu.func(nn = 1000)

CIT.DR.hetero.fitBefore.1000   <- colMeans(CIT.DR.hetero.res.fitBefore, na.rm = TRUE)
CIT.DR.homo.fitBefore.1000     <- colMeans(CIT.DR.homo.res.fitBefore, na.rm = TRUE)
CIT.DR.hetero.fitinparent.1000   <- colMeans(CIT.DR.hetero.res.fitinparent, na.rm = TRUE)
CIT.DR.homo.fitinparent.1000     <- colMeans(CIT.DR.homo.res.fitinparent, na.rm = TRUE)

# CIT.DR.hetero.fitBefore.1000.rf   <- colMeans(CIT.DR.hetero.res.fitBefore, na.rm = TRUE)
# CIT.DR.homo.fitBefore.1000.rf     <- colMeans(CIT.DR.homo.res.fitBefore, na.rm = TRUE)
# CIT.DR.hetero.fitBefore.1000.glm   <- colMeans(CIT.DR.hetero.res.fitBefore, na.rm = TRUE)
# CIT.DR.homo.fitBefore.1000.glm     <- colMeans(CIT.DR.homo.res.fitBefore, na.rm = TRUE)

stopImplicitCluster()



# apply
# CIT.DR.hetero.res.fitBefore   <- apply(sapply(1:sim.reps, function(iter) CIT.DR.hetero.fitBefore.simu.func(nn = 1000)), 1, as.numeric)
# CIT.DR.homo.res.fitBefore     <- apply(sapply(1:sim.reps, function(iter) CIT.DR.homo.fitBefore.simu.func(nn = 1000)), 1, as.numeric)
# CIT.DR.hetero.res.fitinparent <- apply(sapply(1:sim.reps, function(iter) CIT.DR.hetero.fitinparent.simu.func(nn = 1000)), 1, as.numeric)
# CIT.DR.homo.res.fitinparent   <- apply(sapply(1:sim.reps, function(iter) CIT.DR.homo.fitinparent.simu.func(nn = 1000)), 1, as.numeric)


