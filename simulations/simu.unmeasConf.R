set.seed(999)

# foreach
clnum <- detectCores()
cl <- makeCluster(getOption("cl.cores", clnum - 1))
registerDoParallel(cl)
pack <- c("caret", "devtools", "causalTree", "causalForest", "dplyr", "grf",
          "gbm", "MASS", "rpart", "randomForestSRC", "rattle", "rpart.plot")
sim.reps <- 1000

##################################################################################################################
######################################### 1. GCT: Unmeasured confounding #########################################
##################################################################################################################

GCT.unmeasConf.hetero.res.fitBefore <- foreach(1:sim.reps, .combine = rbind, .export = ls(), .packages = pack) %dopar% hetero.fitBefore.simu.func(nn=1000, unmeas.conf = T)


GCT.unmeasConf.homo.res.fitBefore <- foreach(1:sim.reps, .combine = rbind, .export = ls(), .packages = pack) %dopar% homo.fitBefore.simu.func(nn=1000, unmeas.conf = T)


GCT.unmeasConf.hetero.res.fitinparent <- foreach(1:sim.reps, .combine = rbind, .export = ls(), .packages = pack) %dopar% hetero.fitinparent.simu.func(nn=1000, unmeas.conf = T)


GCT.unmeasConf.homo.res.fitinparent <- foreach(1:sim.reps, .combine = rbind, .export = ls(), .packages = pack) %dopar% homo.fitinparent.simu.func(nn=1000, unmeas.conf = T)


# hetero.res.fitBefore.honest <- foreach(1:sim.reps, .combine = rbind, 
#                                   .export = ls(), .packages = pack) %dopar% hetero.fitBefore.honest.simu.func()
# print("5.honest")
# 
# 
# homo.res.fitBefore.honest <- foreach(1:sim.reps, .combine = rbind, 
#                                   .export = ls(), .packages = pack) %dopar% homo.fitBefore.honest.simu.func()
# print("6.honest")
# 
# 
# hetero.res.fitinparent.honest <- foreach(1:sim.reps, .combine = rbind, 
#                                   .export = ls(), .packages = pack) %dopar% hetero.fitinparent.honest.simu.func()
# print("7.honest")
# 
# 
# homo.res.fitinparent.honest <- foreach(1:sim.reps, .combine = rbind, 
#                                   .export = ls(), .packages = pack) %dopar% homo.fitinparent.honest.simu.func()
# print("8.honest")

GCT.unmeasConf.nonhonest <- list(hetero.fitBefore   = colMeans(GCT.unmeasConf.hetero.res.fitBefore, na.rm = TRUE),
                                 homo.fitBefore     = colMeans(GCT.unmeasConf.homo.res.fitBefore, na.rm = TRUE),
                                 hetero.fitinparent = colMeans(GCT.unmeasConf.hetero.res.fitinparent, na.rm = TRUE),
                                 homo.fitinparent   = colMeans(GCT.unmeasConf.homo.res.fitinparent, na.rm = TRUE))

##################################################################################################################
######################################### 2. CIT: Unmeasured confounding #########################################
##################################################################################################################

CIT.DR.unmeasConf.hetero.res.fitBefore <- foreach(1:sim.reps, .combine = rbind, .export = ls(), .packages = pack) %dopar% CIT.DR.hetero.fitBefore.simu.func(nn = 1000, unmeas.conf = T)


CIT.DR.unmeasConf.homo.res.fitBefore <- foreach(1:sim.reps, .combine = rbind, .export = ls(), .packages = pack) %dopar% CIT.DR.homo.fitBefore.simu.func(nn = 1000, unmeas.conf = T)


CIT.DR.unmeasConf.hetero.res.fitinparent <- foreach(1:sim.reps, .combine = rbind, .export = ls(), .packages = pack) %dopar% CIT.DR.hetero.fitinparent.simu.func(nn = 1000, unmeas.conf = T)


CIT.DR.unmeasConf.homo.res.fitinparent <- foreach(1:sim.reps, .combine = rbind, .export = ls(), .packages = pack) %dopar% CIT.DR.homo.fitinparent.simu.func(nn = 1000, unmeas.conf = T)

CIT.DR.unmeasConf <- list(hetero.fitBefore   = colMeans(CIT.DR.unmeasConf.hetero.res.fitBefore, na.rm = TRUE),
                          homo.fitBefore     = colMeans(CIT.DR.unmeasConf.homo.res.fitBefore, na.rm = TRUE),
                          hetero.fitinparent = colMeans(CIT.DR.unmeasConf.hetero.res.fitinparent, na.rm = TRUE),
                          homo.fitinparent   = colMeans(CIT.DR.unmeasConf.homo.res.fitinparent, na.rm = TRUE))


 ##################################################################################################################
######################################### 3. CT: Unmeasured confounding ##########################################
##################################################################################################################

CT.unmeasConf.hetero.res.nonhonest <- foreach(1:sim.reps, .combine = rbind, .export = ls(), .packages = pack) %dopar% CT.hetero.nonhonest.simu.func(nn = 1000, unmeas.conf = T)


CT.unmeasConf.homo.res.nonhonest <- foreach(1:sim.reps, .combine = rbind, .export = ls(), .packages = pack) %dopar% CT.homo.nonhonest.simu.func(nn = 1000, unmeas.conf = T)


# CT.unmeasConf.hetero.res.honest <- foreach(1:sim.reps, .combine = rbind, 
#                                            .export = ls(), .packages = pack) %dopar% CT.hetero.honest.simu.func(nn = 1000, unmeas.conf = T)
# 
# 
# CT.unmeasConf.homo.res.honest <- foreach(1:sim.reps, .combine = rbind, 
#                                          .export = ls(), .packages = pack) %dopar% 

CT.unmeasConf.nonhonest <- list(hetero.fitBefore   = colMeans(CT.unmeasConf.hetero.res.nonhonest, na.rm = TRUE),
                                homo.fitBefore     = colMeans(CT.unmeasConf.homo.res.nonhonest, na.rm = TRUE))


##################################################################################################################
###################################### 4. opt.ct: Unmeasured confounding #########################################
##################################################################################################################
opt.ct.unmeasConf.hetero.res.nonhonest <- foreach(1:sim.reps, .combine = rbind, .export = ls(), .packages = pack) %dopar% opt.ct.hetero.nonhonest.simu.func(nn = 1000, unmeas.conf = T)


opt.ct.unmeasConf.homo.res.nonhonest <- foreach(1:sim.reps, .combine = rbind, .export = ls(), .packages = pack) %dopar% opt.ct.homo.nonhonest.simu.func(nn = 1000, unmeas.conf = T)


 # CT.unmeasConf.hetero.res.honest <- foreach(1:sim.reps, .combine = rbind, 
#                                            .export = ls(), .packages = pack) %dopar% CT.hetero.honest.simu.func(nn = 1000, unmeas.conf = T)
# 
# 
# CT.unmeasConf.homo.res.honest <- foreach(1:sim.reps, .combine = rbind, 
#                                          .export = ls(), .packages = pack) %dopar% CT.homo.honest.simu.func(nn = 1000, unmeas.conf = T)

opt.CT.unmeasConf.nonhonest <- list(hetero.nonhonest = colMeans(opt.ct.unmeasConf.hetero.res.nonhonest, na.rm = TRUE),
                                    homo.fitBefore   = colMeans(opt.ct.unmeasConf.homo.res.nonhonest, na.rm = TRUE))



save(GCT.unmeasConf.nonhonest, CIT.DR.unmeasConf, CT.unmeasConf.nonhonest, opt.CT.unmeasConf.nonhonest,
     file = "E:/project/groupTest/MethodCompare/Simulation/unmeasConf.RData")



stopImplicitCluster()
