n.seq <- c(3000, 5000, 10000, 15000, 20000)
sim.r <- 100

load("./seed.RData")
# increase minsplit and minsize for growing sample size, for example:
# hetero.fitBefore: minsplit = nn/100, minsize = nn/200
# hetero.fittInparent:minsplit = nn/50, minsize = nn/100

# seed: from 4000
# hetero setting: i in 1:length(n.seq) {s in 1:sim.r  seed = 4000 + sim.r * (i-1) + s}
# homo   setting: i in 1:length(n.seq) {s in 1:sim.r  seed = 4000 + sim.r * length(n.seq) + sim.r * (i-1) + s}


# Start Cluster 
clnum <- detectCores()
cl <- makeCluster(getOption("cl.cores", clnum - 4))
registerDoParallel(cl)
pack <- c("caret", "devtools", "causalTree", "dplyr", "gbm", "MASS", "rpart", "randomForestSRC")


## 1. Heterogeneous setting
############################################################################################################
##################################### 1.1. GCT: hetero, nonhonest ##########################################
############################################################################################################
resList.GCT.1 <- list(length(n.seq))
for(i in 1:length(n.seq)) {
  tmp <- foreach(s = 1:sim.r, .combine = rbind, .export = ls(), .packages = pack) %dopar% hetero.fitBefore.nonhonest.simu.func(nn = n.seq[i], seed = seed[4000+sim.r*(i-1)+s])
  resList.GCT.1[[i]] <- Output(tmp)
}

sampleSize.GCT.hetero.fitBefore <- data.frame(Reduce(rbind, resList.GCT.1))


############################################################################################################
##################################### 1.2. GCT: hetero, nonhonest ##########################################
############################################################################################################
resList.GCT.2 <- list(length(n.seq))
for(i in 1:length(n.seq)) {
  tmp <- foreach(s = 1:sim.r, .combine = rbind, .export = ls(), .packages = pack) %dopar% hetero.fitinparent.nonhonest.simu.func(nn = n.seq[i], seed = seed[4000+sim.r*(i-1)+s])
  resList.GCT.2[[i]] <- Output(tmp)
}

sampleSize.GCT.hetero.fitinparent <- data.frame(Reduce(rbind, resList.GCT.2))
rownames(sampleSize.GCT.hetero.fitBefore) <- rownames(sampleSize.GCT.hetero.fitinparent) <- n.seq

save(sampleSize.GCT.hetero.fitBefore, sampleSize.GCT.hetero.fitinparent, file = "E:/GCT/Appendix_SampleSize/seed_res/GCT.hetero.3000_20000.RData")


############################################################################################################
####################################### 1.3. CT: hetero, nonhonest #########################################
############################################################################################################
resList.CT.1 <- list(length(n.seq))
for (i in 1:length(n.seq)) {
  tmp <- foreach(s = 1:sim.r, .combine = rbind, .export = ls(), .packages = pack) %dopar%  CT.hetero.nonhonest.simu.func(nn = n.seq[i], seed = seed[4000+sim.r*(i-1)+s])
  resList.CT.1[[i]] <- colMeans(tmp, na.rm = TRUE)
}

sampleSize.CT.hetero.nonhonest <- data.frame(Reduce(rbind, resList.CT.1))
rownames(sampleSize.CT.hetero.nonhonest) <- n.seq

save(sampleSize.CT.hetero.nonhonest, file = "E:/GCT/Appendix_SampleSize/seed_res/CT.hetero.3000_20000.RData")



############################################################################################################
##################################### 1.4. opt.CT: hetero, nonhonest #######################################
############################################################################################################
resList.optCT.1 <- list(length(n.seq))
for (i in 1:length(n.seq)) {
  tmp <- foreach(s = 1:sim.r, .combine = rbind, .export = ls(), .packages = pack) %dopar% opt.ct.hetero.nonhonest.simu.func(nn = n.seq[i], seed = seed[4000+sim.r*(i-1)+s])
  resList.optCT.1[[i]] <- colMeans(tmp, na.rm = TRUE)
}

sampleSize.opt.ct.hetero.nonhonest <- data.frame(Reduce(rbind, resList.optCT.1))
rownames(sampleSize.opt.ct.hetero.nonhonest) <- n.seq

save(sampleSize.opt.ct.hetero.nonhonest, file = "E:/GCT/Appendix_SampleSize/seed_res/Opt.CT.hetero.3000_20000.RData")




############################################################################################################
###################################### 1.5. CIT: hetero, fitBefore #########################################
############################################################################################################
resList.CIT.1 <- list(length(n.seq))
for (i in 1:length(n.seq)) {
  tmp <- foreach(s = 1:sim.r, .combine = rbind, .export = ls(), .packages = pack) %dopar% CIT.DR.hetero.fitBefore.simu.func(nn = n.seq[i], adj.mthd = "GBM", seed = seed[4000+sim.r*(i-1)+s])
  resList.CIT.1[[i]] <- colMeans(tmp, na.rm = TRUE)
}

sampleSize.CIT.DR.hetero.fitBefore <- data.frame(Reduce(rbind, resList.CIT.1))


############################################################################################################
###################################### 1.6. CIT: hetero, fitinparent #######################################
############################################################################################################
resList.CIT.2 <- list(length(n.seq))
for (i in 1:length(n.seq)) {
  tmp <- foreach(s = 1:sim.r, .combine = rbind, .export = ls(), .packages = pack) %dopar% CIT.DR.hetero.fitinparent.simu.func(nn = n.seq[i], adj.mthd = "GBM", seed = seed[4000+sim.r*(i-1)+s])
  resList.CIT.2[[i]] <- colMeans(tmp, na.rm = TRUE)
}

sampleSize.CIT.DR.hetero.fitinparent <- data.frame(Reduce(rbind, resList.CIT.2))
rownames(sampleSize.CIT.DR.hetero.fitBefore) <- rownames(sampleSize.CIT.DR.hetero.fitinparent) <- n.seq

save(sampleSize.CIT.DR.hetero.fitBefore, sampleSize.CIT.DR.hetero.fitinparent, file = "E:/GCT/Appendix_SampleSize/seed_res/CIT.hetero.3000_20000.RData")

# save(sampleSize.CIT.DR.hetero.fitBefore.GLM, sampleSize.CIT.DR.hetero.fitinparent.GLM, file = "E:/GCT/Appendix_SampleSize/seed_res/CIT.hetero.GLM.3000_20000.RData" )
# save(sampleSize.CIT.DR.hetero.fitBefore.RF, sampleSize.CIT.DR.hetero.fitinparent.RF, file = "E:/GCT/Appendix_SampleSize/seed_res/CIT.hetero.RF.3000_20000.RData" )





## 2. Homogenous Setting
############################################################################################################
##################################### 2.1. GCT: homo, nonhonest ############################################
############################################################################################################
resList.GCT.3 <- list(length(n.seq))
for (i in 1:length(n.seq)) {
  tmp <- foreach(s = 1:sim.r, .combine = rbind, .export = ls(), .packages = pack) %dopar% homo.fitBefore.nonhonest.simu.func(nn = n.seq[i], seed = seed[4000+sim.r*length(n.seq)+sim.r*(i-1)+s])
  resList.GCT.3[[i]] <- Output(tmp, hetero = FALSE)
}

sampleSize.GCT.homo.fitBefore <- data.frame(Reduce(rbind, resList.GCT.3))
rownames(sampleSize.GCT.homo.fitBefore) <- n.seq

############################################################################################################
##################################### 2.2. GCT: homo, nonhonest ############################################
############################################################################################################
resList.GCT.4 <- list(length(n.seq))
for (i in 1:length(n.seq)) {
  tmp <- foreach(s = 1:sim.r, .combine = rbind, .export = ls(), .packages = pack) %dopar% homo.fitinparent.nonhonest.simu.func(nn = n.seq[i], seed = seed[4000+sim.r*length(n.seq)+sim.r*(i-1)+s])
  resList.GCT.4[[i]] <- Output(tmp, hetero = FALSE)
}

sampleSize.GCT.homo.fitinparent <- data.frame(Reduce(rbind, resList.GCT.4))
rownames(sampleSize.GCT.homo.fitinparent) <- n.seq

save(sampleSize.GCT.homo.fitBefore, sampleSize.GCT.homo.fitinparent, file = "E:/GCT/Appendix_SampleSize/seed_res/GCT.homo.3000_20000.RData")


############################################################################################################
######################################## 2.3. CT: homo, nonhonest ##########################################
############################################################################################################
resList.CT.2 <- list(length(n.seq))
for (i in 1:length(n.seq)) {
  tmp <- foreach(s = 1:sim.r, .combine = rbind, .export = ls(), .packages = pack) %dopar% CT.homo.nonhonest.simu.func(nn = n.seq[i], seed = seed[4000+sim.r*length(n.seq)+sim.r*(i-1)+s])
  resList.CT.2[[i]] <- colMeans(tmp, na.rm = TRUE)
}

sampleSize.CT.homo.nonhonest <- data.frame(Reduce(rbind, resList.CT.2))
colnames(sampleSize.CT.homo.nonhonest) <- colnames(tmp)
rownames(sampleSize.CT.homo.nonhonest) <- n.seq
save(sampleSize.CT.homo.nonhonest, file = "E:/GCT/Appendix_SampleSize/seed_res/CT.homo.3000_20000.RData")


############################################################################################################
###################################### 2.4. opt.CT: homo, nonhonest ########################################
############################################################################################################
resList.optCT.2 <- list(length(n.seq))
for (i in 1:length(n.seq)) {
  tmp <- foreach(s = 1:sim.r, .combine = rbind, .export = ls(), .packages = pack) %dopar% opt.ct.homo.nonhonest.simu.func(nn = n.seq[i], seed = seed[4000+sim.r*length(n.seq)+sim.r*(i-1)+s])
  resList.optCT.2[[i]] <- colMeans(tmp, na.rm = TRUE)
}

sampleSize.opt.ct.homo.nonhonest <- data.frame(Reduce(rbind, resList.optCT.2))
colnames(sampleSize.opt.ct.homo.nonhonest) <- colnames(tmp)
rownames(sampleSize.opt.ct.homo.nonhonest) <- n.seq
save(sampleSize.opt.ct.homo.nonhonest, file = "E:/GCT/Appendix_SampleSize/seed_res/Opt.CT.homo.3000_20000.RData")



############################################################################################################
####################################### 2.5. CIT: homo, fitBefore ##########################################
############################################################################################################
resList.CIT.3 <- list(length(n.seq))
for (i in 1:length(n.seq)) {
  tmp <- foreach(s = 1:sim.r, .combine = rbind, .export = ls(), .packages = pack) %dopar% CIT.DR.homo.fitBefore.simu.func(nn = n.seq[i], seed = seed[4000+sim.r*length(n.seq)+sim.r*(i-1)+s])
  resList.CIT.3[[i]] <- colMeans(tmp, na.rm = TRUE)
}

sampleSize.CIT.DR.homo.fitBefore <- data.frame(Reduce(rbind, resList.CIT.3))
colnames(sampleSize.CIT.DR.homo.fitBefore) <- colnames(tmp)  
rownames(sampleSize.CIT.DR.homo.fitBefore) <- n.seq


############################################################################################################
####################################### 2.6. CIT: homo, fitinparent ########################################
############################################################################################################
resList.CIT.4 <- list(length(n.seq))
for (i in 1:length(n.seq)) {
  tmp <- foreach(s = 1:sim.r, .combine = rbind, .export = ls(), .packages = pack) %dopar% CIT.DR.homo.fitinparent.simu.func(nn = n.seq[i], seed = seed[4000+sim.r*length(n.seq)+sim.r*(i-1)+s])
  resList.CIT.4[[i]] <- colMeans(tmp, na.rm = TRUE)
}

sampleSize.CIT.DR.homo.fitinparent <- data.frame(Reduce(rbind, resList.CIT.4))
colnames(sampleSize.CIT.DR.homo.fitinparent) <- colnames(tmp) 
rownames(sampleSize.CIT.DR.homo.fitinparent) <- n.seq

save(sampleSize.CIT.DR.homo.fitBefore, sampleSize.CIT.DR.homo.fitinparent, file = "E:/GCT/Appendix_SampleSize/seed_res/CIT.homo.3000_20000.RData")


# close the clusters
stopImplicitCluster()


# extract the results of three algorithm when sample size = 1000
load("E:/GCT/Study1_MethodComparison/simu.GCT/seed_res/GCT.TreebuildwithVar.nonhonest.1000.RData")
load("E:/GCT/Study1_MethodComparison/simu.CT/seed_res/CT.1000.RData")
load("E:/GCT/Study1_MethodComparison/simu.CT/seed_res/opt.CT.1000.RData")
load("E:/GCT/Study1_MethodComparison/simu.CIT/seed_res/CIT.nonhonest.1000.RData")


sampleSize.GCT.hetero.fitBefore.1000   <- Output(GCT.hetero.res.fitBefore.nonhonest[1:100, ])
sampleSize.GCT.hetero.fitinparent.1000 <- Output(GCT.hetero.res.fitinparent.nonhonest[1:100, ])
sampleSize.GCT.homo.fitBefore.1000     <- Output(GCT.homo.res.fitBefore.nonhonest[1:100, ], hetero = F)
sampleSize.GCT.homo.fitinparent.1000   <- Output(GCT.homo.res.fitinparent.nonhonest[1:100, ], hetero = F)


sampleSize.CIT.DR.hetero.fitBefore.1000   <- colMeans(CIT.DR.hetero.res.fitBefore[1:100, ], na.rm = T)
sampleSize.CIT.DR.hetero.fitinparent.1000 <- colMeans(CIT.DR.hetero.res.fitinparent[1:100, ], na.rm = T)
sampleSize.CIT.DR.homo.fitBefore.1000     <- colMeans(CIT.DR.homo.res.fitBefore[1:100, ], na.rm = T)
sampleSize.CIT.DR.homo.fitinparent.1000   <- colMeans(CIT.DR.homo.res.fitinparent[1:100, ], na.rm = T)


sampleSize.CT.hetero.nonhonest.1000      <- colMeans(CT.hetero.res.nonhonest[1:100, ], na.rm = T)
sampleSize.CT.homo.nonhonest.1000        <- colMeans(CT.homo.res.nonhonest[1:100, ], na.rm = T)
sampleSize.opt.ct.hetero.nonhonest.1000  <- colMeans(opt.ct.hetero.res.nonhonest[1:100, ], na.rm = T)
sampleSize.opt.ct.homo.nonhonest.1000    <- colMeans(opt.ct.homo.res.nonhonest[1:100, ], na.rm = T)


sampleSize.GCT.hetero.fitBefore   <- rbind(sampleSize.GCT.hetero.fitBefore.1000, sampleSize.GCT.hetero.fitBefore)
sampleSize.GCT.homo.fitBefore     <- rbind(sampleSize.GCT.homo.fitBefore.1000, sampleSize.GCT.homo.fitBefore)
sampleSize.GCT.hetero.fitinparent <- rbind(sampleSize.GCT.hetero.fitinparent.1000, sampleSize.GCT.hetero.fitinparent)
sampleSize.GCT.homo.fitinparent   <- rbind(sampleSize.GCT.homo.fitinparent.1000, sampleSize.GCT.homo.fitinparent)


sampleSize.CIT.DR.hetero.fitBefore   <- rbind(sampleSize.CIT.DR.hetero.fitBefore.1000, sampleSize.CIT.DR.hetero.fitBefore)
sampleSize.CIT.DR.homo.fitBefore     <- rbind(sampleSize.CIT.DR.homo.fitBefore.1000, sampleSize.CIT.DR.homo.fitBefore)
sampleSize.CIT.DR.hetero.fitinparent <- rbind(sampleSize.CIT.DR.hetero.fitinparent.1000, sampleSize.CIT.DR.hetero.fitinparent)
sampleSize.CIT.DR.homo.fitinparent   <- rbind(sampleSize.CIT.DR.homo.fitinparent.1000, sampleSize.CIT.DR.homo.fitinparent)


sampleSize.CT.hetero.nonhonest     <- rbind(sampleSize.CT.hetero.nonhonest.1000, sampleSize.CT.hetero.nonhonest)
sampleSize.CT.homo.nonhonest       <- rbind(sampleSize.CT.homo.nonhonest.1000, sampleSize.CT.homo.nonhonest)
sampleSize.opt.ct.hetero.nonhonest <- rbind(sampleSize.opt.ct.hetero.nonhonest.1000, sampleSize.opt.ct.hetero.nonhonest)
sampleSize.opt.ct.homo.nonhonest   <- rbind(sampleSize.opt.ct.homo.nonhonest.1000, sampleSize.opt.ct.homo.nonhonest)

rownames(sampleSize.GCT.hetero.fitBefore) <- rownames(sampleSize.GCT.homo.fitBefore) <- rownames(sampleSize.GCT.hetero.fitinparent) <-
  rownames(sampleSize.GCT.homo.fitinparent) <- rownames(sampleSize.CIT.DR.hetero.fitBefore) <- rownames(sampleSize.CIT.DR.homo.fitBefore) <- 
  rownames(sampleSize.CIT.DR.hetero.fitinparent) <- rownames(sampleSize.CIT.DR.homo.fitinparent) <- rownames(sampleSize.CT.hetero.nonhonest) <-
  rownames(sampleSize.CT.homo.nonhonest) <- rownames(sampleSize.opt.ct.hetero.nonhonest) <- rownames(sampleSize.opt.ct.homo.nonhonest) <- n.seq   

save(sampleSize.GCT.hetero.fitBefore, sampleSize.CT.hetero.nonhonest, 
     sampleSize.opt.ct.hetero.nonhonest, sampleSize.CIT.DR.hetero.fitBefore,
     sampleSize.GCT.homo.fitBefore, sampleSize.CT.homo.nonhonest, 
     sampleSize.opt.ct.homo.nonhonest, sampleSize.CIT.DR.homo.fitBefore,
     sampleSize.GCT.hetero.fitinparent, sampleSize.GCT.homo.fitinparent,
     sampleSize.CIT.DR.hetero.fitinparent, sampleSize.CIT.DR.homo.fitinparent,
     file = "E:/GCT/Appendix_SampleSize/seed_res/sampleSize.RData")
