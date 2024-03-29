# seed: from 5000
# hetero: {i in 1:sim.reps  seed = 5000 + i}
# homo:   {i in 1:sim.reps  seed = 6000 + i}

setwd("..")
load("./seed.RData")

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

GCT.unmeasConf.hetero.res.fitBefore <- foreach(i = 1:sim.reps, .combine = rbind, .export = ls(), .packages = pack) %dopar% 
  hetero.fitBefore.nonhonest.simu.func(nn=1000, unmeas.conf = T, seed = seed[5000+i])
GCT.unmeasConf.homo.res.fitBefore   <- foreach(i = 1:sim.reps, .combine = rbind, .export = ls(), .packages = pack) %dopar% 
  homo.fitBefore.nonhonest.simu.func(nn=1000, unmeas.conf = T, seed = seed[5000+sim.reps+i])
GCT.unmeasConf.hetero.res.fitinparent <- foreach(i = 1:sim.reps, .combine = rbind, .export = ls(), .packages = pack) %dopar% 
  hetero.fitinparent.nonhonest.simu.func(nn=1000, unmeas.conf = T, seed = seed[5000+i])
GCT.unmeasConf.homo.res.fitinparent   <- foreach(i = 1:sim.reps, .combine = rbind, .export = ls(), .packages = pack) %dopar% 
  homo.fitinparent.nonhonest.simu.func(nn=1000, unmeas.conf = T, seed = seed[5000+sim.reps+i])

GCT.hetero.fitBefore.nonhonest.unmsConf   <- Output(GCT.unmeasConf.hetero.res.fitBefore)
GCT.hetero.fitinparent.nonhonest.unmsConf <- Output(GCT.unmeasConf.hetero.res.fitinparent)
GCT.homo.fitBefore.nonhonest.unmsConf     <- Output(GCT.unmeasConf.homo.res.fitBefore, hetero = FALSE)
GCT.homo.fitinparent.nonhonest.unmsConf   <- Output(GCT.unmeasConf.homo.res.fitinparent, hetero = FALSE)

GCT.unmeasConf <- list(hetero.fitBefore   = GCT.hetero.fitBefore.nonhonest.unmsConf,
                       homo.fitBefore     = GCT.homo.fitBefore.nonhonest.unmsConf,
                       hetero.fitinparent = GCT.hetero.fitinparent.nonhonest.unmsConf,
                       homo.fitinparent   = GCT.homo.fitinparent.nonhonest.unmsConf)
save(GCT.unmeasConf, file = "./Simulations/Appendix_UnmeasConf/res/GCT.unmeasConf.nonhonest.1000.RData" )





##################################################################################################################
######################################### 2. CIT: Unmeasured confounding #########################################
##################################################################################################################

CIT.DR.unmeasConf.hetero.res.fitBefore <- foreach(i = 1:sim.reps, .combine = rbind, .export = ls(), .packages = pack) %dopar% CIT.DR.hetero.fitBefore.simu.func(nn = 1000, unmeas.conf = T, seed = seed[5000+i])
CIT.DR.unmeasConf.homo.res.fitBefore   <- foreach(i = 1:sim.reps, .combine = rbind, .export = ls(), .packages = pack) %dopar% CIT.DR.homo.fitBefore.simu.func(nn = 1000, unmeas.conf = T, seed = seed[5000+sim.reps+i])
CIT.DR.unmeasConf.hetero.res.fitinparent <- foreach(i = 1:sim.reps, .combine = rbind, .export = ls(), .packages = pack) %dopar% CIT.DR.hetero.fitinparent.simu.func(nn = 1000, unmeas.conf = T, seed = seed[5000+i])
CIT.DR.unmeasConf.homo.res.fitinparent   <- foreach(i = 1:sim.reps, .combine = rbind, .export = ls(), .packages = pack) %dopar% CIT.DR.homo.fitinparent.simu.func(nn = 1000, unmeas.conf = T, seed = seed[5000+sim.reps+i])

CIT.DR.unmeasConf <- list(hetero.fitBefore   = colMeans(CIT.DR.unmeasConf.hetero.res.fitBefore, na.rm = TRUE),
                          homo.fitBefore     = colMeans(CIT.DR.unmeasConf.homo.res.fitBefore, na.rm = TRUE),
                          hetero.fitinparent = colMeans(CIT.DR.unmeasConf.hetero.res.fitinparent, na.rm = TRUE),
                          homo.fitinparent   = colMeans(CIT.DR.unmeasConf.homo.res.fitinparent, na.rm = TRUE))
save(CIT.DR.unmeasConf, file = "./Simulations/Appendix_UnmeasConf/res/CIT.DR.unmeasConf.1000.RData" )





##################################################################################################################
######################################### 3. CT: Unmeasured confounding ##########################################
##################################################################################################################

CT.unmeasConf.hetero.res.nonhonest <- foreach(i = 1:sim.reps, .combine = rbind, .export = ls(), .packages = pack) %dopar% CT.hetero.nonhonest.simu.func(nn = 1000, unmeas.conf = T, seed = seed[5000+i])
CT.unmeasConf.homo.res.nonhonest   <- foreach(i = 1:sim.reps, .combine = rbind, .export = ls(), .packages = pack) %dopar% CT.homo.nonhonest.simu.func(nn = 1000, unmeas.conf = T, seed = seed[5000+sim.reps+i])


CT.unmeasConf.nonhonest <- list(hetero.fitBefore   = colMeans(CT.unmeasConf.hetero.res.nonhonest, na.rm = TRUE),
                                homo.fitBefore     = colMeans(CT.unmeasConf.homo.res.nonhonest, na.rm = TRUE))

save(CT.unmeasConf.nonhonest, file = "./Simulations/Appendix_UnmeasConf/res/CT.unmeasConf.nonhonest.1000.RData" )






##################################################################################################################
###################################### 4. opt.ct: Unmeasured confounding #########################################
##################################################################################################################
opt.ct.unmeasConf.hetero.res.nonhonest <- foreach(i = 1:sim.reps, .combine = rbind, .export = ls(), .packages = pack) %dopar% opt.ct.hetero.nonhonest.simu.func(nn = 1000, unmeas.conf = T, seed = seed[5000+i])
opt.ct.unmeasConf.homo.res.nonhonest   <- foreach(i = 1:sim.reps, .combine = rbind, .export = ls(), .packages = pack) %dopar% opt.ct.homo.nonhonest.simu.func(nn = 1000, unmeas.conf = T, seed = seed[5000+sim.reps+i])


opt.CT.unmeasConf.nonhonest <- list(hetero.nonhonest = colMeans(opt.ct.unmeasConf.hetero.res.nonhonest, na.rm = TRUE),
                                    homo.fitBefore   = colMeans(opt.ct.unmeasConf.homo.res.nonhonest, na.rm = TRUE))

save(opt.CT.unmeasConf.nonhonest, file = "./Simulations/Appendix_UnmeasConf/res/Opt.CT.unmeasConf.nonhonest.1000.RData" )


gc()
stopImplicitCluster()
