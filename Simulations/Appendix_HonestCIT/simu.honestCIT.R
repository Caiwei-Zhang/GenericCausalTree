# #!/usr/bin/env Rscript
# setwd("./")
# folder <- paste(getwd(), "/Functions/", sep="")
# functions <- list.files(folder)
# functions <- grep("main.R", functions, invert = T, value = T)
# functions <- paste(folder, functions, sep = "")
# for (i in functions){
#   source(i)
# }

setwd("..")
load("./seed.RData")

setwd("./Simulations/Study1_MethodComparison/CIT-master/Functions/")  
functions <- list.files(getwd())
functions <- grep("main.R", functions, invert = T, value = T)
functions <- paste(getwd(), functions, sep = "/")
for (i in functions){
  source(i)
}

source("./Functions/dgp.R")




sim.reps <- 1000
boot.R <- 1000
propsc.mthd <- "GLM"
adj.mthd    <- "RF"
type.var    <- "cont"

clnum <- detectCores()
cl <- makeCluster(getOption("cl.cores", clnum - 4))
registerDoParallel(cl)
pack <- c("caret", "devtools", "gbm", "dplyr", "rpart", "MASS", "randomForestSRC", "plm", "xgboost")

res.all <-
  foreach(i = 1:sim.reps, .combine  = "rbind", .export = ls(), .packages = pack) %dopar% {


# sim.CIT.ci <- function() {
# 
#   tab <- NULL
  
  # for(i in 1:sim.reps){
    set.seed(seed[i + 2000]) 
    
    # data generation
    data.lst <- hetero.DataGen(N = 2000, p = 5, dgp = "cov", N.test = 1000)
    data.lst$where.split <- list(c(1))   
    data.lst$dir.split  <- list(c(NULL))
    
    data.used.full.mixed  <- data.lst$data[1:1000, ]
    data.test <- data.lst$data.test
  
    
    data.used.mixed       <- data.used.full.mixed[1:800, ]
    
    # val.sample: used in EstIpw.CvMethod1, the order of the columns must be A, Y, X
    data.validation.mixed <- data.used.full.mixed[801:1000, ]  
    split.point <- numeric(2) 
    ########################################################################################################
    ####################################### 1 hetero performance ###########################################
    ########################################################################################################
    
    ## 1.1 hetero fitBefore performance
    t0 <- Sys.time()
    seq.created.estdr.fitBefore <- create.sequence(data.used         = data.used.mixed,
                                                   est.used          = "DR",
                                                   type.var          = "cont",
                                                   propsc.mod.loc    = "out",
                                                   propsc.mthd       = "GLM",
                                                   propsc.form.true  = NULL,
                                                   adj.mod.loc       = "out", 
                                                   adj.mthd          = adj.mthd, 
                                                   adj.form.true     = NULL, 
                                                   num.truc.obs      = 30,
                                                   min.node          = 10)
  
    final.tree.estdr.fitBefore <- EstDr.CvMethod1(data.used         = data.used.mixed,
                                                  tree.list         = seq.created.estdr.fitBefore$tree.list, 
                                                  lambda.used       = qchisq(0.95, 1),
                                                  val.sample        = data.validation.mixed,
                                                  type.var          = "cont",
                                                  propsc.mod.loc    = "out",
                                                  propsc.mthd       = "GLM",
                                                  propsc.form.true  = NULL,
                                                  adj.mod.loc       = "out",
                                                  adj.mthd          = adj.mthd,
                                                  adj.form.true     = NULL)
    t1 <- Sys.time()
    
    eval.final.estdr.fitBefore <- eval.measures.eff(final.tree   = final.tree.estdr.fitBefore[[1]],
                                                    test.data    = data.test,
                                                    true.trt.eff = data.lst$true.trt.eff,
                                                    noise.var    = data.lst$noise.splt.var,
                                                    corr.split   = data.lst$eff.modifier,
                                                    #type.eff.mod  = data.lst$type.eff.mod,
                                                    where.split  = data.lst$where.split,
                                                    dir.split    = data.lst$dir.split)
    
  
    eval.final.estdr.fitBefore$time <- as.numeric(difftime(t1, t0, units = "secs"))
    
    if (eval.final.estdr.fitBefore$exact.corr == 1){
      split.point[1] <- final.tree.estdr.fitBefore[[1]]$splits[4] #[,4] split point
    } else {
      split.point[1] <- 0
    }
    
    
    ## 1.2 hetero fitinparent performance 
    t0 <- Sys.time()
    seq.created.estdr.fitinparent <- create.sequence(data.used       = data.used.mixed,
                                                   est.used          = "DR",
                                                   type.var          = "cont",
                                                   propsc.mod.loc    = "node",
                                                   propsc.mthd       = "GLM",
                                                   propsc.form.true  = NULL,
                                                   adj.mod.loc       = "node", 
                                                   adj.mthd          = adj.mthd, 
                                                   adj.form.true     = NULL, 
                                                   num.truc.obs      = 60,
                                                   min.node          = 30)
    
    final.tree.estdr.fitinparent <- EstDr.CvMethod1(data.used       = data.used.mixed,
                                                  tree.list         = seq.created.estdr.fitinparent$tree.list, 
                                                  lambda.used       = qchisq(0.95, 1),
                                                  val.sample        = data.validation.mixed,
                                                  type.var          = "cont",
                                                  propsc.mod.loc    = "node",
                                                  propsc.mthd       = "GLM",
                                                  propsc.form.true  = NULL,
                                                  adj.mod.loc       = "node",
                                                  adj.mthd          = adj.mthd,
                                                  adj.form.true     = NULL)
    t1 <- Sys.time()
    
    eval.final.estdr.fitinparent <- eval.measures.eff(final.tree   = final.tree.estdr.fitinparent[[1]],
                                                    test.data    = data.test,
                                                    true.trt.eff = data.lst$true.trt.eff,
                                                    noise.var    = data.lst$noise.splt.var,
                                                    corr.split   = data.lst$eff.modifier,
                                                    #type.eff.mod  = data.lst$type.eff.mod,
                                                    where.split  = data.lst$where.split,
                                                    dir.split    = data.lst$dir.split)
    
    
    eval.final.estdr.fitinparent$time <- as.numeric(difftime(t1, t0, units = "secs"))
    
    
    if (eval.final.estdr.fitinparent$exact.corr == 1){
      split.point[2] <- final.tree.estdr.fitinparent[[1]]$splits[4] #[,4] split point
    } else {
      split.point[2] <- 0
    }
    
    
    performance.hetero.drInnd <- unlist(c(eval.final.estdr.fitBefore, 
                                   eval.final.estdr.fitinparent))
    
    
    ########################################################################################################
    ######################################### 2 hetero coverage ############################################
    ########################################################################################################
    data.used.full.est <- (data.lst$data)[1001:2000, ]
    hetero.trt.eff <- matrix(NA, nrow = boot.R, ncol = 2*2)
        
    # Calculate the treatment effect and CI later anyway since might be used for a correct tree in either final tree selection method
    # will only calculate coverage using those iterations that get correct trees
    for (boot.i in 1:boot.R) {
      
      # seet seed for each bootstrap sample so can be replicated for each parameter setting
      set.seed(seed[boot.i])
      data.boot <- data.used.full.est[sample.int(1000, size = 1000, replace = T), ]
      
      ################################################################################
      ######################### fit models on the whole dataset ######################
      ################################################################################
      data.l.out <- data.boot %>% filter(X4 <= split.point[1]) # true.trt.eff = 2
      data.r.out <- data.boot %>% filter(X4 > split.point[1])  # true.trt.eff = 5
      
      data.boot.noy     <- data.boot[, !colnames(data.boot) %in% c("Y")]
      
      whl.propsc <- gen.fullrank.ipw(df.noy = data.boot.noy, propsc.form.true = NULL)
      whl.propsc <- withWarnings(est.prop.sc(df.noy    = whl.propsc$df.noy.fullrank,
                                             method    = propsc.mthd,
                                             form.true = NULL))
      
      
      whl.g <- gen.fullrank.g(df = data.boot, adj.form.true = NULL)
      whl.g <- withWarnings(est.cond.eff(df        = whl.g$df.fullrank,
                                         method    = adj.mthd,
                                         form.true = NULL,
                                         type.var  = type.var))
      
      propsc.l <- whl.propsc$value$prop.sc[data.boot$X4 <= split.point[1]]
      propsc.r <- whl.propsc$value$prop.sc[data.boot$X4 > split.point[1]]
      
      est.cond.eff.0.l <- whl.g$value$pred.A.0[data.boot$X4 <= split.point[1]]
      est.cond.eff.1.l <- whl.g$value$pred.A.1[data.boot$X4 <= split.point[1]]
      
      est.cond.eff.0.r <- whl.g$value$pred.A.0[data.boot$X4 > split.point[1]]
      est.cond.eff.1.r <- whl.g$value$pred.A.1[data.boot$X4 > split.point[1]]
      
      
      hetero.trt.eff[boot.i, 1] <- mean(data.l.out$A * (data.l.out$Y - est.cond.eff.1.l) / propsc.l + est.cond.eff.1.l) - 
        mean((1 - data.l.out$A) * (data.l.out$Y - est.cond.eff.0.l) / (1 - propsc.l) + est.cond.eff.0.l)
      hetero.trt.eff[boot.i, 2] <- mean(data.r.out$A * (data.r.out$Y - est.cond.eff.1.r) / propsc.r + est.cond.eff.1.r) - 
        mean((1 - data.r.out$A) * (data.r.out$Y - est.cond.eff.0.r) / (1 - propsc.r) + est.cond.eff.0.r)
    
        
      ################################################################################
      ####################### fit models in terminal nodes ###########################
      ################################################################################
      data.l <- data.boot %>% filter(X4 <= split.point[2])  # true.trt.eff = 2
      data.r <- data.boot %>% filter(X4 > split.point[2])   # true.trt.eff = 5
      
      data.l.noy  <- data.l[, !colnames(data.l) %in% c("Y")]
      
      left.propsc <- gen.fullrank.ipw(df.noy           = data.l.noy, 
                                      propsc.form.true = NULL)
      left.propsc <- withWarnings(est.prop.sc(df.noy    = left.propsc$df.noy.fullrank,
                                              method    = propsc.mthd,
                                              form.true = NULL))
      
      left.g <- gen.fullrank.g(df            = data.l,
                               adj.form.true = NULL)
      left.g <- withWarnings(est.cond.eff(df        = left.g$df.fullrank,
                                          method    = adj.mthd,
                                          form.true = NULL,
                                          type.var  = type.var))
      
      
      data.r.noy     <- data.r[, !colnames(data.r) %in% c("Y")]
      
      right.propsc <- gen.fullrank.ipw(df.noy           = data.r.noy, 
                                       propsc.form.true = NULL)
      right.propsc <- withWarnings(est.prop.sc(df.noy    = right.propsc$df.noy.fullrank,
                                               method    = propsc.mthd,
                                               form.true = NULL))
      
      
      right.g <- gen.fullrank.g(df            = data.r,
                                adj.form.true = NULL)
      right.g <- withWarnings(est.cond.eff(df        = right.g$df.fullrank,
                                           method    = adj.mthd,
                                           form.true = NULL,
                                           type.var  = type.var))
      
      
      ################################################################################
      ##################### find propensity scores in terminal nodes #################
      ################################################################################
      # if 1) propensity score model fitting produces warning or 2) only produce thresholded propensity scores
      # use outside model fitting results
      if (!is.null(left.propsc$warnings)) {
        
        cond.warnings <- T
        for (warnings.i in 1:length(left.propsc$warnings)) {
          cond.warnings <- cond.warnings & (left.propsc$warnings[[warnings.i]]$message == "prediction from a rank-deficient fit may be misleading")
        }
        
        # when there is only rank-deficient warning, use the local prediction
        if (cond.warnings) {
          propsc.l <- left.propsc$value$prop.sc
        } else {
          propsc.l <- whl.propsc$value$prop.sc[data.boot$X4 <= split.point[2]]
        }
        
      } else if (identical(sort(unique(left.propsc$value$prop.sc)), c(0.1, 0.9))) {
        propsc.l <- whl.propsc$value$prop.sc[data.boot$X4 <= split.point[2]]
        
      } else {
        propsc.l <- left.propsc$value$prop.sc
      }
      
      
      # when the outcome model fitting produces warning
      # use outside model fitting results for conditional means
      if (!is.null(left.g$warnings)) {
        
        # when there is only rank-deficient warning, use the local prediction
        cond.warnings <- T
        for (warnings.i in 1:length(left.g$warnings)) {
          cond.warnings <- cond.warnings & (left.g$warnings[[warnings.i]]$message == "prediction from a rank-deficient fit may be misleading")
        }
        
        if (cond.warnings) {
          est.cond.eff.0.l <- left.g$value$pred.A.0
          est.cond.eff.1.l <- left.g$value$pred.A.1
        } else {
          est.cond.eff.0.l <- whl.g$value$pred.A.0[data.boot$X4 <= split.point[2]]
          est.cond.eff.1.l <- whl.g$value$pred.A.1[data.boot$X4 <= split.point[2]]
        }
        
      } else { # if (!is.null(left.g$warnings)) {
        est.cond.eff.0.l <- left.g$value$pred.A.0
        est.cond.eff.1.l <- left.g$value$pred.A.1
      }
      
      
      ################################################################################
      ################# find propensity scores in right terminal node ################
      ################################################################################
      if (!is.null(right.propsc$warnings)) {
        
        cond.warnings <- T
        for (warnings.i in 1:length(right.propsc$warnings)) {
          cond.warnings <- cond.warnings & (right.propsc$warnings[[warnings.i]]$message == "prediction from a rank-deficient fit may be misleading")
        }
        
        # when there is only rank-deficient warning, use the local prediction
        if (cond.warnings) {
          propsc.r <- right.propsc$value$prop.sc
        } else {
          propsc.r <- whl.propsc$value$prop.sc[data.boot$X4 > split.point[2]]
        }
        
      } else if (identical(sort(unique(right.propsc$value$prop.sc)), c(0.1, 0.9))) {
        propsc.r <- whl.propsc$value$prop.sc[data.boot$X4 > split.point[2]]
        
      } else {
        propsc.r <- right.propsc$value$prop.sc
      }
      
      
      ################################################################################
      ################ find potential outcomes in right terminal node ################
      ################################################################################
      # when the outcome model fitting produces warning
      # use outside model fitting results for conditional means
      if (!is.null(right.g$warnings)) {
        
        # when there is only rank-deficient warning, use the local prediction
        cond.warnings <- T
        for (warnings.i in 1:length(right.g$warnings)) {
          cond.warnings <- cond.warnings & (right.g$warnings[[warnings.i]]$message == "prediction from a rank-deficient fit may be misleading")
        }
        
        if (cond.warnings) {
          est.cond.eff.0.r <- right.g$value$pred.A.0
          est.cond.eff.1.r <- right.g$value$pred.A.1
        } else {
          est.cond.eff.0.r <- whl.g$value$pred.A.0[data.boot$X4 > split.point[2]]
          est.cond.eff.1.r <- whl.g$value$pred.A.1[data.boot$X4 > split.point[2]]
        }
        
      } else { # if (!is.null(right.g$warnings)) {
        est.cond.eff.0.r <- right.g$value$pred.A.0
        est.cond.eff.1.r <- right.g$value$pred.A.1
      }
      
      
      ################################################################################
      ################## Find treatment effects in two terminal nodes ################
      ################################################################################
      
      tmp.l <- mean(data.l$A * (data.l$Y - est.cond.eff.1.l) / propsc.l + est.cond.eff.1.l) - 
        mean((1 - data.l$A) * (data.l$Y - est.cond.eff.0.l) / (1 - propsc.l) + est.cond.eff.0.l)
      tmp.r <- mean(data.r$A * (data.r$Y - est.cond.eff.1.r) / propsc.r + est.cond.eff.1.r) - 
        mean((1 - data.r$A) * (data.r$Y - est.cond.eff.0.r) / (1 - propsc.r) + est.cond.eff.0.r)
      
      if(abs(tmp.l-2) < abs(tmp.l-5) | abs(tmp.r-2) > abs(tmp.r-5)) { # hetero.trt.eff[, 3] true trt.eff = 2; hetero.trt.eff[, 4] true trt.eff = 5
        
        hetero.trt.eff[boot.i, 3] <- tmp.l
        hetero.trt.eff[boot.i, 4] <- tmp.r
        
      } else {
        
        hetero.trt.eff[boot.i, 3] <- tmp.r
        hetero.trt.eff[boot.i, 4] <- tmp.l
        
      }
      
    } # for (boot.i in 1:boot.R) {
    
    res <- as.vector(apply(hetero.trt.eff, 2, quantile, probs = c(0.025, 0.5, 0.975)))
    
    # hetero.trt.eff[boot.i, 3] <- mean(data.l$A * (data.l$Y - est.cond.eff.1.l) / propsc.l + est.cond.eff.1.l) - 
    #   mean((1 - data.l$A) * (data.l$Y - est.cond.eff.0.l) / (1 - propsc.l) + est.cond.eff.0.l)
    # hetero.trt.eff[boot.i, 4] <- mean(data.r$A * (data.r$Y - est.cond.eff.1.r) / propsc.r + est.cond.eff.1.r) - 
    #   mean((1 - data.r$A) * (data.r$Y - est.cond.eff.0.r) / (1 - propsc.r) + est.cond.eff.0.r)
    
    
    
    ##########################################################################################################
    ####################################### 3 homo performance ###############################################
    ##########################################################################################################
    set.seed(seed[i+3000]) 
    
    
    # data generation
    data.lst <- homo.DataGen(N = 2000, p = 5, dgp = "cov", N.test = 1000)
    data.lst$where.split <- list(c(1))   
    data.lst$dir.split  <- list(c(NULL))
    
    data.used.full.homo  <- data.lst$data[1:1000, ]
    data.test.homo <- data.lst$data.test
    
    
    data.used.homo      <- data.used.full.homo[1:800, ]
    
    # val.sample: used in EstIpw.CvMethod1, the order of the columns must be A, Y, X
    data.validation.homo <- data.used.full.homo[801:1000, ]  
    
    
    ## homo fitBefore
    t0 <- Sys.time()
    seq.created.estdr.fitBefore <- create.sequence(data.used         = data.used.homo,
                                                   est.used          = "DR",
                                                   type.var          = "cont",
                                                   propsc.mod.loc    = "out",
                                                   propsc.mthd       = "GLM",
                                                   propsc.form.true  = NULL,
                                                   adj.mod.loc       = "out", 
                                                   adj.mthd          = adj.mthd, 
                                                   adj.form.true     = NULL, 
                                                   num.truc.obs      = 100,
                                                   min.node          = 50)
    
    final.tree.estdr.fitBefore <- EstDr.CvMethod1(data.used         = data.used.homo,
                                                  tree.list         = seq.created.estdr.fitBefore$tree.list, 
                                                  lambda.used       = qchisq(0.95, 1),
                                                  val.sample        = data.validation.homo,
                                                  type.var          = "cont",
                                                  propsc.mod.loc    = "out",
                                                  propsc.mthd       = "GLM",
                                                  propsc.form.true  = NULL,
                                                  adj.mod.loc       = "out",
                                                  adj.mthd          = adj.mthd,
                                                  adj.form.true     = NULL)
    t1 <- Sys.time()
    
    eval.final.estdr.fitBefore <- eval.measures.eff(final.tree   = final.tree.estdr.fitBefore[[1]],
                                                    test.data    = data.test.homo,
                                                    true.trt.eff = data.lst$true.trt.eff,
                                                    noise.var    = data.lst$noise.splt.var,
                                                    corr.split   = data.lst$eff.modifier,
                                                    #type.eff.mod  = data.lst$type.eff.mod,
                                                    where.split  = data.lst$where.split,
                                                    dir.split    = data.lst$dir.split)
    
    
    eval.final.estdr.fitBefore$time <- as.numeric(difftime(t1, t0, units = "secs"))
    
    
    
    ## homo fitinparent
    t0 <- Sys.time()
    seq.created.estdr.fitinparent <- create.sequence(data.used         = data.used.homo,
                                                     est.used          = "DR",
                                                     type.var          = "cont",
                                                     propsc.mod.loc    = "node",
                                                     propsc.mthd       = "GLM",
                                                     propsc.form.true  = NULL,
                                                     adj.mod.loc       = "node", 
                                                     adj.mthd          = adj.mthd, 
                                                     adj.form.true     = NULL, 
                                                     num.truc.obs      = 100,
                                                     min.node          = 50)
    
    final.tree.estdr.fitinparent <- EstDr.CvMethod1(data.used         = data.used.homo,
                                                    tree.list         = seq.created.estdr.fitinparent$tree.list, 
                                                    lambda.used       = qchisq(0.95, 1),
                                                    val.sample        = data.validation.homo,
                                                    type.var          = "cont",
                                                    propsc.mod.loc    = "node",
                                                    propsc.mthd       = "GLM",
                                                    propsc.form.true  = NULL,
                                                    adj.mod.loc       = "node",
                                                    adj.mthd          = adj.mthd,
                                                    adj.form.true     = NULL)
    t1 <- Sys.time()
    
    eval.final.estdr.fitinparent <- eval.measures.eff(final.tree   = final.tree.estdr.fitinparent[[1]],
                                                      test.data    = data.test.homo,
                                                      true.trt.eff = data.lst$true.trt.eff,
                                                      noise.var    = data.lst$noise.splt.var,
                                                      corr.split   = data.lst$eff.modifier,
                                                      #type.eff.mod  = data.lst$type.eff.mod,
                                                      where.split  = data.lst$where.split,
                                                      dir.split    = data.lst$dir.split)
    
    
    eval.final.estdr.fitinparent$time <- as.numeric(difftime(t1, t0, units = "secs"))
    
    performance.homo.drInnd <- unlist(c(eval.final.estdr.fitBefore, eval.final.estdr.fitinparent))
    
    
    
    performance <- c(performance.hetero.drInnd, performance.homo.drInnd) # 16+16, dim = 32
    
    
    ##########################################################################################################
    ######################################### 3 homo coverage ################################################
    ##########################################################################################################
    homo.trt.eff <- matrix(NA, nrow = boot.R, ncol = 1)
    data.used.full.est <- data.lst$data[1001:2000, ]
    
    # Calculate the treatment effect and CI later anyway since might be used for a correct tree in either final tree selection method
    # will only calculate coverage using those iterations that get correct trees
    for (boot.i in 1:boot.R) {
      
      # set seed for each bootstrap sample so can be replicated for each parameter setting
      set.seed(seed[boot.i])
      
      data.boot <- data.used.full.est[sample.int(1000, size = 1000, replace = T), ]
      
      ################################################################################
      ############## fit the propensity score model on the whole dataset #############
      ################################################################################
      data.boot.noy     <- data.boot[, !colnames(data.boot) %in% c("Y")]
      
      whl.propsc <- gen.fullrank.ipw(df.noy           = data.boot.noy, 
                                     propsc.form.true = NULL)
      whl.propsc <- withWarnings(est.prop.sc(df.noy    = whl.propsc$df.noy.fullrank,
                                             method    = propsc.mthd,
                                             form.true =NULL))
      
      
      ################################################################################
      ################## fit the outcome model on the whole dataset ##################
      ################################################################################
      whl.g <- gen.fullrank.g(df            = data.boot,
                              adj.form.true = NULL)
      whl.g <- withWarnings(est.cond.eff(df        = whl.g$df.fullrank,
                                         method    = adj.mthd,
                                         form.true = NULL,
                                         type.var  = type.var))
      
      ################################################################################
      ################## Find treatment effects in two terminal nodes ################
      ################################################################################
      homo.trt.eff[boot.i, 1] <- mean(data.boot$A * (data.boot$Y - whl.g$value$pred.A.1) / whl.propsc$value$prop.sc + whl.g$value$pred.A.1) - 
        mean((1 - data.boot$A) * (data.boot$Y - whl.g$value$pred.A.0) / (1 - whl.propsc$value$prop.sc) + whl.g$value$pred.A.0)
      
    } # for (boot.i in 1:boot.R) {
    
    
    res <- c(res, as.vector(apply(homo.trt.eff, 2, quantile, probs = c(0.025, 0.5, 0.975)))) # 12+3, dim = 15
    
    # print must be put before output, otherwise the output will be print output
    if (i%%10 == 0) {print(i)}
  
    #tab <- rbind(tab, c(performance, res))
  
  # }
  # tab
  c(performance, res) # dim = 32+15 = 47
  
 
}


# performance <- res.all[, 1:32]
# CIT.DR.1000 <- rbind(performance[8 * 0 + 1:8],
#                      performance[8 * 1 + 1:8],
#                      performance[8 * 2 + 1:8],
#                      performance[8 * 3 + 1:8])
# rownames(CIT.DR.1000) <- c("CIT.DR.hetero.fitBefore.1000", 
#                            "CIT.DR.hetero.fitBeinparent.1000",
#                            "CIT.DR.homo.fitBefore.1000", 
#                            "CIT.DR.homo.fitinparent.1000")
# 
# save(performance.drInnd, CIT.DR.1000, file = "E:/GCT/Study1_MethodComparison/simu.CIT/res/gbm_res/CIT.DR.1000.RData")

  
performance.drInnd <- res.all[, 1:32]
ci.dr <- res.all[, 33:ncol(res.all)]

save(performance.drInnd, file = "./Simulations/Study2_honestGCT/res/CIT.DR.ci.performance.RData")
save(ci.dr, file = "./Simulations/Study2_honestGCT/res/CIT.DR.ci.RData")

gc()


dim(performance.drInnd)
performance.all <- rbind(performance.drInnd[, 8 * 0 + 1:8],
                         performance.drInnd[, 8 * 1 + 1:8],
                         performance.drInnd[, 8 * 2 + 1:8],
                         performance.drInnd[, 8 * 3 + 1:8])



mode <- c("fitBe", "fitIn")
setting   <- c("Heterogeneous", "Homogeneous")
Method   <- c("Both unknown DR-CIT")
scnrs.cit <- expand.grid(Method, setting, mode)
scnrs.cit <- scnrs.cit[c(1,3,2,4),]
colnames(scnrs.cit) <- c("Method", "Setting", "Mode")

performance.all <- cbind(scnrs.cit[rep(1:nrow(scnrs.cit), each = sim.reps), ],
                         performance.all)
dim(performance.all)


# confidence intervals
ci.dr <- rbind(ci.dr[, 6 * 0 + 1:6], #fitBe
               ci.dr[, 6 * 1 + 1:6], #fitIn
               cbind(ci.dr[, 6 * 2 + 3 * 0 + 1:3], matrix(NA, nrow = sim.reps, ncol = 3)), # fitBe or fitIn
               cbind(ci.dr[, 6 * 2 + 3 * 0 + 1:3], matrix(NA, nrow = sim.reps, ncol = 3))) # fitBe or fitIn
ci.all <- ci.dr
dim(ci.all)

performance.all <- cbind(performance.all, ci.all)
colnames(performance.all)[(ncol(performance.all) - 6 + 1):ncol(performance.all)] <- c(paste0("l", c("025", "50", "975")),
                                                                                      paste0("r", c("025", "50", "975")))
#colnames(performance.all) <- gsub("hetero.dr.adjTGlmInnd.propscTGlmInnd.cv1.", "", colnames(performance.all))

performance.all <- performance.all %>%
  mutate(l025 = ifelse(exact.corr == 0, NA, l025),
         l50  = ifelse(exact.corr == 0, NA, l50),
         l975 = ifelse(exact.corr == 0, NA, l975),
         r025 = ifelse(exact.corr == 0, NA, r025),
         r50  = ifelse(exact.corr == 0, NA, r50),
         r975 = ifelse(exact.corr == 0, NA, r975))

performance.all <- performance.all %>%
  mutate(cover_l = ((l025 <= 2) & (l975 >= 2)),
         cover_r = ((r025 <= 5) & (r975 >= 5)))


summ_perf <- performance.all %>% group_by(Setting, Mode) %>% summarise(mse        = mean(mse, na.rm = TRUE),
                                                               exact.corr = mean(exact.corr, na.rm = TRUE),
                                                               size.tree  = mean(size.tree, na.rm = TRUE),
                                                               num.splt   = mean(num.splt, na.rm = TRUE),
                                                               fir.splt   = mean(fir.splt, na.rm = TRUE),
                                                               numb.noise = mean(numb.noise, na.rm = TRUE),
                                                               pps        = mean(pps, na.rm = TRUE),
                                                               time       = mean(time, na.rm = TRUE),
                                                               cov.l      = mean(cover_l, na.rm = TRUE), 
                                                               cov.r      = mean(cover_r, na.rm = TRUE), 
                                                               .group     = FALSE)

save(performance.all, summ_perf, file = "./Simulations/Study2_honestGCT/res/CIT.DR.performance+ci.summ.RData")










# registerDoParallel(cores=(Sys.getenv("SLURM_NTASKS_PER_NODE")))
# 
# # Read in the arguments from command line
# option_list = list(
#   make_option(c("-s", "--start"), type="integer"),
#   make_option(c("-e", "--end"), type = "integer"))
# 
# opt_parser <- OptionParser(option_list=option_list)
# opt <- parse_args(opt_parser)
# start <- opt$start
# end   <- opt$end


# ci.dr <-
#   foreach(i = 1:sim.reps, 
#           .combine  = "rbind", .export = ls(), .packages = pack) %dopar% {
#   
#   #tab <- NULL
#   
#   #for (i in 1:sim.reps) {
# 
#             
#       # newly generated data should have seed different from before
#       set.seed(seed[i])  
#       
#       data.mixed            <-  hetero.DataGen(N = 1000, p = 5, dgp = "cov", N.test = 1000)
#       data.used.full.mixed  <- data.mixed$data
#       data.used.mixed       <- data.used.full.mixed[1:800, ]
#       # val.sample: used in EstIpw.CvMethod1, the order of the columns must be A, Y, X
#       data.validation.mixed <- data.used.full.mixed[801:1000, ]  
#       
#       # / 4 because only hetero and both final tree selection method should have the same treatment effect estimate
#       # * 2 because two terminal nodes
#       hetero.trt.eff <- matrix(NA, nrow = boot.R, ncol = 2 * 2)
#       
#       # Calculate the treatment effect and CI later anyway since might be used for a correct tree in either final tree selection method
#       # will only calculate coverage using those iterations that get correct trees
#       for (boot.i in 1:boot.R) {
#       
#         # seet seed for each bootstrap sample so can be replicated for each parameter setting
#         set.seed(seed[boot.i])
#         data.boot <- data.used.full.mixed[sample.int(1000, size = 1000, replace = T), ]
#         
#         data.l <- data.boot %>% filter(X4 <= 0)
#         data.r <- data.boot %>% filter(X4 > 0)
#         
#         ################################################################################
#         ######################### fit models on the whole dataset ######################
#         ################################################################################
#         data.boot.noy     <- data.boot[, !colnames(data.boot) %in% c("Y")]
#         
#         whl.propsc <- gen.fullrank.ipw(df.noy = data.boot.noy,propsc.form.true = NULL)
#         whl.propsc <- withWarnings(est.prop.sc(df.noy    = whl.propsc$df.noy.fullrank,
#                                                  method    = propsc.mthd,
#                                                  form.true = NULL))
# 
#         
#         whl.g <- gen.fullrank.g(df = data.boot, adj.form.true = NULL)
#         whl.g <- withWarnings(est.cond.eff(df        = whl.g$df.fullrank,
#                                            method    = adj.mthd,
#                                            form.true = NULL,
#                                            type.var  = type.var))
#           
#         propsc.l <- whl.propsc$value$prop.sc[data.boot$X4 <= 0]
#         propsc.r <- whl.propsc$value$prop.sc[data.boot$X4 > 0]
#         
#         est.cond.eff.0.l <- whl.g$value$pred.A.0[data.boot$X4 <= 0]
#         est.cond.eff.1.l <- whl.g$value$pred.A.1[data.boot$X4 <= 0]
#         
#         est.cond.eff.0.r <- whl.g$value$pred.A.0[data.boot$X4 > 0]
#         est.cond.eff.1.r <- whl.g$value$pred.A.1[data.boot$X4 > 0]
#         
#         
#         hetero.trt.eff[boot.i, 1] <- mean(data.l$A * (data.l$Y - est.cond.eff.1.l) / propsc.l + est.cond.eff.1.l) - 
#           mean((1 - data.l$A) * (data.l$Y - est.cond.eff.0.l) / (1 - propsc.l) + est.cond.eff.0.l)
#         hetero.trt.eff[boot.i, 2] <- mean(data.r$A * (data.r$Y - est.cond.eff.1.r) / propsc.r + est.cond.eff.1.r) - 
#           mean((1 - data.r$A) * (data.r$Y - est.cond.eff.0.r) / (1 - propsc.r) + est.cond.eff.0.r)
#           
# 
#         ################################################################################
#         ####################### fit models in terminal nodes ###########################
#         ################################################################################
#         data.l.noy     <- data.l[, !colnames(data.l) %in% c("Y")]
#         
#         left.propsc <- gen.fullrank.ipw(df.noy           = data.l.noy, 
#                                         propsc.form.true = NULL)
#         left.propsc <- withWarnings(est.prop.sc(df.noy    = left.propsc$df.noy.fullrank,
#                                                   method    = propsc.mthd,
#                                                   form.true = NULL))
# 
#         left.g <- gen.fullrank.g(df            = data.l,
#                                  adj.form.true = NULL)
#         left.g <- withWarnings(est.cond.eff(df        = left.g$df.fullrank,
#                                             method    = adj.mthd,
#                                             form.true = NULL,
#                                             type.var  = type.var))
# 
#         
#         data.r.noy     <- data.r[, !colnames(data.r) %in% c("Y")]
#         
#         right.propsc <- gen.fullrank.ipw(df.noy           = data.r.noy, 
#                                          propsc.form.true = NULL)
#         right.propsc <- withWarnings(est.prop.sc(df.noy    = right.propsc$df.noy.fullrank,
#                                                    method    = propsc.mthd,
#                                                    form.true = NULL))
#         
# 
#         right.g <- gen.fullrank.g(df            = data.r,
#                                   adj.form.true = NULL)
#         right.g <- withWarnings(est.cond.eff(df        = right.g$df.fullrank,
#                                              method    = adj.mthd,
#                                              form.true = NULL,
#                                              type.var  = type.var))
#         
#         
#         ################################################################################
#         ##################### find propensity scores in terminal nodes #################
#         ################################################################################
#         # if 1) propensity score model fitting produces warning or 2) only produce thresholded propensity scores
#         # use outside model fitting results
#         if (!is.null(left.propsc$warnings)) {
#           
#           cond.warnings <- T
#           for (warnings.i in 1:length(left.propsc$warnings)) {
#             cond.warnings <- cond.warnings & (left.propsc$warnings[[warnings.i]]$message == "prediction from a rank-deficient fit may be misleading")
#           }
#           
#           # when there is only rank-deficient warning, use the local prediction
#           if (cond.warnings) {
#             propsc.l <- left.propsc$value$prop.sc
#           } else {
#             propsc.l <- whl.propsc$value$prop.sc[data.boot$X4 <= 0]
#           }
#           
#         } else if (identical(sort(unique(left.propsc$value$prop.sc)), c(0.1, 0.9))) {
#           propsc.l <- whl.propsc$value$prop.sc[data.boot$X4 <= 0]
#           
#         } else {
#           propsc.l <- left.propsc$value$prop.sc
#         }
#         
#         
#         # when the outcome model fitting produces warning
#         # use outside model fitting results for conditional means
#         if (!is.null(left.g$warnings)) {
#           
#           # when there is only rank-deficient warning, use the local prediction
#           cond.warnings <- T
#           for (warnings.i in 1:length(left.g$warnings)) {
#             cond.warnings <- cond.warnings & (left.g$warnings[[warnings.i]]$message == "prediction from a rank-deficient fit may be misleading")
#           }
#           
#           if (cond.warnings) {
#             est.cond.eff.0.l <- left.g$value$pred.A.0
#             est.cond.eff.1.l <- left.g$value$pred.A.1
#           } else {
#             est.cond.eff.0.l <- whl.g$value$pred.A.0[data.boot$X4 <= 0]
#             est.cond.eff.1.l <- whl.g$value$pred.A.1[data.boot$X4 <= 0]
#           }
#           
#         } else { # if (!is.null(left.g$warnings)) {
#           est.cond.eff.0.l <- left.g$value$pred.A.0
#           est.cond.eff.1.l <- left.g$value$pred.A.1
#         }
#         
# 
#         ################################################################################
#         ################# find propensity scores in right terminal node ################
#         ################################################################################
#         if (!is.null(right.propsc$warnings)) {
#           
#           cond.warnings <- T
#           for (warnings.i in 1:length(right.propsc$warnings)) {
#             cond.warnings <- cond.warnings & (right.propsc$warnings[[warnings.i]]$message == "prediction from a rank-deficient fit may be misleading")
#           }
#           
#           # when there is only rank-deficient warning, use the local prediction
#           if (cond.warnings) {
#             propsc.r <- right.propsc$value$prop.sc
#           } else {
#             propsc.r <- whl.propsc$value$prop.sc[data.boot$X4 > 0]
#           }
#           
#         } else if (identical(sort(unique(right.propsc$value$prop.sc)), c(0.1, 0.9))) {
#           propsc.r <- whl.propsc$value$prop.sc[data.boot$X4 > 0]
#           
#         } else {
#           propsc.r <- right.propsc$value$prop.sc
#         }
#         
#         
#         ################################################################################
#         ################ find potential outcomes in right terminal node ################
#         ################################################################################
#         # when the outcome model fitting produces warning
#         # use outside model fitting results for conditional means
#         if (!is.null(right.g$warnings)) {
#           
#           # when there is only rank-deficient warning, use the local prediction
#           cond.warnings <- T
#           for (warnings.i in 1:length(right.g$warnings)) {
#             cond.warnings <- cond.warnings & (right.g$warnings[[warnings.i]]$message == "prediction from a rank-deficient fit may be misleading")
#           }
#           
#           if (cond.warnings) {
#             est.cond.eff.0.r <- right.g$value$pred.A.0
#             est.cond.eff.1.r <- right.g$value$pred.A.1
#           } else {
#             est.cond.eff.0.r <- whl.g$value$pred.A.0[data.boot$X4 > 0]
#             est.cond.eff.1.r <- whl.g$value$pred.A.1[data.boot$X4 > 0]
#           }
#           
#         } else { # if (!is.null(right.g$warnings)) {
#           est.cond.eff.0.r <- right.g$value$pred.A.0
#           est.cond.eff.1.r <- right.g$value$pred.A.1
#         }
#         
# 
#         ################################################################################
#         ################## Find treatment effects in two terminal nodes ################
#         ################################################################################
# 
#         tmp.l <- mean(data.l$A * (data.l$Y - est.cond.eff.1.l) / propsc.l + est.cond.eff.1.l) - 
#           mean((1 - data.l$A) * (data.l$Y - est.cond.eff.0.l) / (1 - propsc.l) + est.cond.eff.0.l)
#         tmp.r <- mean(data.r$A * (data.r$Y - est.cond.eff.1.r) / propsc.r + est.cond.eff.1.r) - 
#           mean((1 - data.r$A) * (data.r$Y - est.cond.eff.0.r) / (1 - propsc.r) + est.cond.eff.0.r)
#         
#         if(abs(tmp.l-2) < abs(tmp.l-5) | abs(tmp.r-2) > abs(tmp.r-5)) { # hetero.trt.eff[, 3] true trt.eff = 2
#           
#           hetero.trt.eff[boot.i, 3] <- tmp.l
#           hetero.trt.eff[boot.i, 4] <- tmp.r
#           
#         } else {
#            
#           hetero.trt.eff[boot.i, 3] <- tmp.r
#           hetero.trt.eff[boot.i, 4] <- tmp.l
#           
#         }
#         
#         # hetero.trt.eff[boot.i, 3] <- mean(data.l$A * (data.l$Y - est.cond.eff.1.l) / propsc.l + est.cond.eff.1.l) - 
#         #   mean((1 - data.l$A) * (data.l$Y - est.cond.eff.0.l) / (1 - propsc.l) + est.cond.eff.0.l)
#         # hetero.trt.eff[boot.i, 4] <- mean(data.r$A * (data.r$Y - est.cond.eff.1.r) / propsc.r + est.cond.eff.1.r) - 
#         #   mean((1 - data.r$A) * (data.r$Y - est.cond.eff.0.r) / (1 - propsc.r) + est.cond.eff.0.r)
#         
#        
#       
#       } # for (boot.i in 1:boot.R) {
#             
#             
#       res <- as.vector(apply(hetero.trt.eff, 2, quantile, probs = c(0.025, 0.5, 0.975)))
#             
#       #####################################################################################################################
#       ################################################### Homogeneous #####################################################
#       #####################################################################################################################
#       # newly generated data should have seed different from before
#       set.seed(seed[i])  
#       
#       data.mixed            <- homo.DataGen(N = 1000, p = 5, dgp = "cov", N.test = 1000)
#       data.used.full.mixed  <- data.mixed$data
#       data.used.mixed       <- data.used.full.mixed[1:800, ]
#       # val.sample: used in EstIpw.CvMethod1, the order of the columns must be A, Y, X
#       data.validation.mixed <- data.used.full.mixed[801:1000, ]  
#       
#       # / 4 because only hetero and both final tree selection method should have the same treatment effect estimate
#       homo.trt.eff <- matrix(NA, nrow = boot.R, ncol = 1)
#       
#       # Calculate the treatment effect and CI later anyway since might be used for a correct tree in either final tree selection method
#       # will only calculate coverage using those iterations that get correct trees
#       for (boot.i in 1:boot.R) {
#         
#         # seet seed for each bootstrap sample so can be replicated for each parameter setting
#         set.seed(seed[boot.i])
#         
#         data.boot <- data.used.full.mixed[sample.int(1000, size = 1000, replace = T), ]
#         
#         ################################################################################
#         ############## fit the propensity score model on the whole dataset #############
#         ################################################################################
#         data.boot.noy     <- data.boot[, !colnames(data.boot) %in% c("Y")]
#         
#         whl.propsc <- gen.fullrank.ipw(df.noy           = data.boot.noy, 
#                                        propsc.form.true = NULL)
#         whl.propsc <- withWarnings(est.prop.sc(df.noy    = whl.propsc$df.noy.fullrank,
#                                                method    = propsc.mthd,
#                                                form.true =NULL))
# 
#         
#         ################################################################################
#         ################## fit the outcome model on the whole dataset ##################
#         ################################################################################
#         whl.g <- gen.fullrank.g(df            = data.boot,
#                                 adj.form.true = NULL)
#         whl.g <- withWarnings(est.cond.eff(df        = whl.g$df.fullrank,
#                                            method    = adj.mthd,
#                                            form.true = NULL,
#                                            type.var  = type.var))
#         
#         ################################################################################
#         ################## Find treatment effects in two terminal nodes ################
#         ################################################################################
#         # Both True
#         homo.trt.eff[boot.i, 1] <- mean(data.boot$A * (data.boot$Y - whl.g$value$pred.A.1) / whl.propsc$value$prop.sc + whl.g$value$pred.A.1) - 
#           mean((1 - data.boot$A) * (data.boot$Y - whl.g$value$pred.A.0) / (1 - whl.propsc$value$prop.sc) + whl.g$value$pred.A.0)
#         
#       } # for (boot.i in 1:boot.R) {
#         
#         
#       res <- c(res, as.vector(apply(homo.trt.eff, 2, quantile, probs = c(0.025, 0.5, 0.975))))
#         
#       # print must be put before output, otherwise the output will be print output
#       if (i%%30 == 0) {print(i)}
#       
#       res
#         
#     }
# 
# gc()
# memory.size()


# # output the result
# load("../BinMixedDr.RData")
# dim(performance.drInnd)
# performance.dr <- performance.drInnd
# 
# performance.all <- rbind(performance.dr[, 8 * 0 + 1:8],
#                          performance.dr[, 8 * 1 + 1:8],
#                          performance.dr[, 8 * 2 + 1:8],
#                          performance.dr[, 8 * 3 + 1:8])
# 
# rm(performance.dr,performance.drInnd)
# 
# 
# mode <- c("fitBe", "fitIn")
# setting   <- c("Heterogeneous", "Homogeneous")
# Method   <- c("Both unknown DR-CIT")
# scnrs.cit <- expand.grid(Method, setting, mode)
# scnrs.cit <- scnrs.cit[c(1,3,2,4),]
# colnames(scnrs.cit) <- c("Method", "Setting", "Mode")
# 
# performance.all <- cbind(scnrs.cit[rep(1:nrow(scnrs.cit), each = sim.reps), ],
#                          performance.all)
# 
# # confidence intervals
# load("../BinMixedDrCover.RData")
# 
# 
# ci.dr <- rbind(ci.dr[, 6 * 0 + 1:6], #fitBew
#                ci.dr[, 6 * 1 + 1:6], #fitIn
#                cbind(ci.dr[, 6 * 2 + 3 * 0 + 1:3], matrix(NA, nrow = sim.reps, ncol = 3)), # fitBe or fitIn
#                cbind(ci.dr[, 6 * 2 + 3 * 0 + 1:3], matrix(NA, nrow = sim.reps, ncol = 3))) # fitBe or fitIn
# ci.all <- ci.dr
# dim(ci.all)
# 
# performance.all <- cbind(performance.all, ci.all)
# colnames(performance.all)[(ncol(performance.all) - 6 + 1):ncol(performance.all)] <- c(paste0("l", c("025", "50", "975")),
#                                                                                       paste0("r", c("025", "50", "975")))
# #colnames(performance.all) <- gsub("hetero.dr.adjTGlmInnd.propscTGlmInnd.cv1.", "", colnames(performance.all))
# 
# performance.all <- performance.all %>%
#   mutate(l025 = ifelse(exact.corr == 0, NA, l025),
#          l50  = ifelse(exact.corr == 0, NA, l50),
#          l975  = ifelse(exact.corr == 0, NA, l975),
#          r025 = ifelse(exact.corr == 0, NA, r025),
#          r50  = ifelse(exact.corr == 0, NA, r50),
#          r975  = ifelse(exact.corr == 0, NA, r975))
# 
# performance.all <- performance.all %>%
#   mutate(cover_l = (l025 <= 2) & (l975 >= 2),
#          cover_r = (r025 <= 5) & (r975 >= 5))
# 
# performance.all %>% group_by(Setting, Mode) %>% summarise(cov.l = mean(cover_l, na.rm = TRUE), 
#                                                           cov.r = mean(cover_r, na.rm = TRUE), .group = FALSE)


