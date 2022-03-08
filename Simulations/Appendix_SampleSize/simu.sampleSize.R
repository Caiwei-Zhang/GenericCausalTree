# increase minsplit and minsize for growing sample size
n.seq <- c(3000, 5000, 10000, 15000, 20000)
sim.r <- 100

# Start Cluster 
clnum <- detectCores()
cl <- makeCluster(getOption("cl.cores", clnum - 4))
registerDoParallel(cl)
pack <- c("caret", "devtools", "causalTree", "dplyr", "gbm", "MASS", "rpart", "randomForestSRC")


## 1. Heterogeneous setting
############################################################################################################
##################################### 1.1. GCT: hetero, nonhonest ##########################################
############################################################################################################
resList <- list(length(n.seq))
for(i in 1:length(n.seq)) {
  tmp <- foreach(1:sim.r, .combine = rbind, .export = ls(), .packages = pack) %dopar% hetero.fitBefore.nonhonest.simu.func(nn = n.seq[i])
  resList[[i]] <- Output(tmp)
}

sampleSize.GCT.hetero.fitBefore <- data.frame(Reduce(rbind, resList))


############################################################################################################
##################################### 1.2. GCT: hetero, nonhonest ##########################################
############################################################################################################
resList <- list(length(n.seq))
for(i in 1:length(n.seq)) {
  tmp <- foreach(1:sim.r, .combine = rbind, .export = ls(), .packages = pack) %dopar% hetero.fitinparent.nonhonest.simu.func(nn = n.seq[i])
  resList[[i]] <- Output(tmp)
}

sampleSize.GCT.hetero.fitinparent <- data.frame(Reduce(rbind, resList))


save(sampleSize.GCT.hetero.fitBefore, sampleSize.GCT.hetero.fitinparent, file = "E:/GCT/Appendix_SampleSize/res/GCT.hetero.15000_20000.RData")


############################################################################################################
####################################### 1.3. CT: hetero, nonhonest #########################################
############################################################################################################
resList <- list(length(n.seq))
for (i in 1:length(n.seq)) {
  tmp <- foreach(1:sim.r, .combine = rbind, .export = ls(), .packages = pack) %dopar%  CT.hetero.nonhonest.simu.func(nn = n.seq[i])
  resList[[i]] <- colMeans(tmp, na.rm = TRUE)
}

sampleSize.CT.hetero.nonhonest <- data.frame(Reduce(rbind, resList))

save(sampleSize.CT.hetero.nonhonest, file = "E:/GCT/Appendix_SampleSize/res/CT.hetero.15000_20000.RData")



############################################################################################################
##################################### 1.4. opt.CT: hetero, nonhonest #######################################
############################################################################################################
resList <- list(length(n.seq))
for (i in 1:length(n.seq)) {
  tmp <- foreach(1:sim.r, .combine = rbind, .export = ls(), .packages = pack) %dopar% opt.ct.hetero.nonhonest.simu.func(nn = n.seq[i])
  resList[[i]] <- colMeans(tmp, na.rm = TRUE)
}

sampleSize.opt.ct.hetero.nonhonest <- data.frame(Reduce(rbind, resList))

save(sampleSize.opt.ct.hetero.nonhonest, file = "E:/GCT/Appendix_SampleSize/res/Opt.CT.hetero.15000_20000.RData")




############################################################################################################
###################################### 1.5. CIT: hetero, fitBefore #########################################
############################################################################################################
resList <- list(length(n.seq))
for (i in 1:length(n.seq)) {
  tmp <- foreach(1:sim.r, .combine = rbind, .export = ls(), .packages = pack) %dopar% CIT.DR.hetero.fitBefore.simu.func(nn = n.seq[i])
  resList[[i]] <- colMeans(tmp, na.rm = TRUE)
}

sampleSize.CIT.DR.hetero.fitBefore <- data.frame(Reduce(rbind, resList))



############################################################################################################
###################################### 1.6. CIT: hetero, fitinparent #######################################
############################################################################################################
resList <- list(length(n.seq))
for (i in 1:length(n.seq)) {
  tmp <- foreach(1:sim.r, .combine = rbind, .export = ls(), .packages = pack) %dopar% CIT.DR.hetero.fitinparent.simu.func(nn = n.seq[i])
  resList[[i]] <- colMeans(tmp, na.rm = TRUE)
}

sampleSize.CIT.DR.hetero.fitinparent <- data.frame(Reduce(rbind, resList))

save(sampleSize.CIT.DR.hetero.fitBefore, sampleSize.CIT.DR.hetero.fitinparent, file = "E:/GCT/Appendix_SampleSize/res/CIT.hetero.15000_20000.RData")






## 2. Homogenous Setting
############################################################################################################
##################################### 2.1. GCT: homo, nonhonest ############################################
############################################################################################################
resList <- list(length(n.seq))
for (i in 1:length(n.seq)) {
  tmp <- foreach(1:sim.r, .combine = rbind, .export = ls(), .packages = pack) %dopar% homo.fitBefore.nonhonest.simu.func(nn = n.seq[i])
  resList[[i]] <- Output(tmp, hetero = FALSE)
}

sampleSize.GCT.homo.fitBefore <- data.frame(Reduce(rbind, resList))



############################################################################################################
##################################### 2.2. GCT: homo, nonhonest ############################################
############################################################################################################
resList <- list(length(n.seq))
for (i in 1:length(n.seq)) {
  tmp <- foreach(1:sim.r, .combine = rbind, .export = ls(), .packages = pack) %dopar% homo.fitinparent.nonhonest.simu.func(nn = n.seq[i])
  resList[[i]] <- Output(tmp, hetero = FALSE)
}

sampleSize.GCT.homo.fitinparent <- data.frame(Reduce(rbind, resList))

save(sampleSize.GCT.homo.fitBefore, sampleSize.GCT.homo.fitinparent, file = "E:/GCT/Appendix_SampleSize/res/GCT.homo.3000_10000.RData")


############################################################################################################
######################################## 2.3. CT: homo, nonhonest ##########################################
############################################################################################################
resList <- list(length(n.seq))
for (i in 1:length(n.seq)) {
  tmp <- foreach(1:sim.r, .combine = rbind, .export = ls(), .packages = pack) %dopar%  CT.homo.nonhonest.simu.func(nn = n.seq[i])
  resList[[i]] <- colMeans(tmp, na.rm = TRUE)
}

sampleSize.CT.homo.nonhonest <- data.frame(Reduce(rbind, resList))
colnames(sampleSize.CT.homo.nonhonest) <- colnames(tmp)
save(sampleSize.CT.homo.nonhonest, file = "E:/GCT/Appendix_SampleSize/res/CT.homo.3000_10000.RData")


############################################################################################################
###################################### 2.4. opt.CT: homo, nonhonest ########################################
############################################################################################################
resList <- list(length(n.seq))
for (i in 1:length(n.seq)) {
  tmp <- foreach(1:sim.r, .combine = rbind, .export = ls(), .packages = pack) %dopar%  opt.ct.homo.nonhonest.simu.func(nn = n.seq[i])
  resList[[i]] <- colMeans(tmp, na.rm = TRUE)
}

sampleSize.opt.ct.homo.nonhonest <- data.frame(Reduce(rbind, resList))
colnames(sampleSize.opt.ct.homo.nonhonest) <- colnames(tmp)
save(sampleSize.opt.ct.homo.nonhonest, file = "E:/GCT/Appendix_SampleSize/res/Opt.CT.homo.3000_10000.RData")



############################################################################################################
####################################### 2.5. CIT: homo, fitBefore ##########################################
############################################################################################################
resList <- list(length(n.seq))
for (i in 1:length(n.seq)) {
  tmp <- foreach(1:sim.r, .combine = rbind, .export = ls(), .packages = pack) %dopar% CIT.DR.homo.fitBefore.simu.func(nn = n.seq[i])
  resList[[i]] <- colMeans(tmp, na.rm = TRUE)
}

sampleSize.CIT.DR.homo.fitBefore <- data.frame(Reduce(rbind, resList))
colnames(sampleSize.CIT.DR.homo.fitBefore) <- colnames(tmp)  



############################################################################################################
####################################### 2.6. CIT: homo, fitinparent ########################################
############################################################################################################
resList <- list(length(n.seq))
for (i in 1:length(n.seq)) {
  tmp <- foreach(1:sim.r, .combine = rbind, .export = ls(), .packages = pack) %dopar% CIT.DR.homo.fitinparent.simu.func(nn = n.seq[i])
  resList[[i]] <- colMeans(tmp, na.rm = TRUE)
}

sampleSize.CIT.DR.homo.fitinparent <- data.frame(Reduce(rbind, resList))
colnames(sampleSize.CIT.DR.homo.fitinparent) <- colnames(tmp) 

save(sampleSize.CIT.DR.homo.fitBefore, sampleSize.CIT.DR.homo.fitinparent, file = "E:/GCT/Appendix_SampleSize/res/CIT.homo.3000_10000.RData")


# close the clusters
stopImplicitCluster()




save(sampleSize.GCT.hetero.fitBefore, sampleSize.CT.hetero.nonhonest, 
     sampleSize.opt.ct.hetero.nonhonest, sampleSize.CIT.DR.hetero.fitBefore,
     sampleSize.GCT.homo.fitBefore, sampleSize.CT.homo.nonhonest, 
     sampleSize.opt.ct.homo.nonhonest, sampleSize.CIT.DR.homo.fitBefore,
     sampleSize.GCT.hetero.fitinparent, sampleSize.GCT.homo.fitinparent,
     sampleSize.CIT.DR.hetero.fitinparent, sampleSize.CIT.DR.homo.fitinparent,
     file = "E:/Project/groupTest/MethodCompare/Simulation/sampleSize.RData")




### integrate the results of all methods in different sample size
subset.GCT.hetero.fitBefore    <- sampleSize.GCT.hetero.fitBefore["mse"] %>% mutate(Method = "GCT.fitB", n = n.seq, honest = "nonhonest", hetero = "heterogeneous") 
subset.CIT.DR.hetero.fitBefore <- sampleSize.CIT.DR.hetero.fitBefore["mse"] %>% mutate(Method = "CIT.DR.fitB", n = n.seq, honest = "nonhonest", hetero = "heterogeneous") 
subset.CT.hetero.nonhonest     <- sampleSize.CT.hetero.nonhonest["mse"] %>% mutate(Method = "CT", n = n.seq, honest = "nonhonest", hetero = "heterogeneous") 
subset.opt.ct.hetero.nonhonest <- sampleSize.opt.ct.hetero.nonhonest["mse"] %>% mutate(Method = "opt.ct", n = n.seq, honest = "nonhonest", hetero = "heterogeneous") 
subset.GCT.homo.fitBefore      <- sampleSize.GCT.homo.fitBefore["mse"] %>% mutate(Method = "GCT.fitB", n = n.seq, honest = "nonhonest", hetero = "homogeneous") 
subset.CIT.DR.homo.fitBefore   <- sampleSize.CIT.DR.homo.fitBefore["mse"] %>% mutate(Method = "CIT.DR.fitB", n = n.seq, honest = "nonhonest", hetero = "homogeneous") 
subset.CT.homo.nonhonest       <- sampleSize.CT.homo.nonhonest["mse"] %>% mutate(Method = "CT", n = n.seq, honest = "nonhonest", hetero = "homogeneous") 
subset.opt.ct.homo.nonhonest   <- sampleSize.opt.ct.homo.nonhonest["mse"] %>% mutate(Method = "opt.ct", n = n.seq, honest = "nonhonest", hetero = "homogeneous") 

subset.GCT.hetero.fitinparent <- sampleSize.GCT.hetero.fitinparent["mse"] %>% mutate(Method = "GCT.fitB", n = n.seq, honest = "nonhonest", hetero = "heterogeneous") 
subset.CIT.DR.hetero.fitinparent <- sampleSize.CIT.DR.hetero.fitin+parent["mse"] %>% mutate(Method = "GCT.fitB", n = n.seq, honest = "nonhonest", hetero = "heterogeneous") 
subset.GCT.homo.fitinparent <- sampleSize.GCT.homo.fitinparent["mse"] %>% mutate(Method = "GCT.fitin", n = n.seq, honest = "nonhonest", hetero = "homogeneous") 
subset.CIT.DR.homo.fitinparent <- sampleSize.CIT.DR.homo.fitinparent["mse"] %>% mutate(Method = "GCT.fitin", n = n.seq, honest = "nonhonest", hetero = "homogeneous")  


res.sampleSize <- rbind(subset.GCT.hetero.fitBefore,   
                        subset.CIT.DR.hetero.fitBefore,
                        subset.CT.hetero.nonhonest,    
                        subset.opt.ct.hetero.nonhonest,
                        subset.GCT.homo.fitBefore,     
                        subset.CIT.DR.homo.fitBefore, 
                        subset.CT.homo.nonhonest,      
                        subset.opt.ct.homo.nonhonest) 


res.sampleSize$hetero <- factor(res.sampleSize$hetero, levels = c(0, 1), labels=c("heterogeneous", "homogeneous"))
res.sampleSize$honest <- factor(res.sampleSize$honest, levels = c(0, 1), labels=c("nonhonest", "honest"))
res.sampleSize$method <- factor(res.sampleSize$method)

library(ggplot2)
ggplot(data = res.sampleSize, aes(x = n, y = mse, shape = Method, color = Method)) +
  geom_line(size = 3) +
  facet_grid(hetero~honest) +
  labs(title = "MSE by Different Methods", x = "sample size", y = "MSE")


subset.honest.GCT.hetero.fitBefore <- sampleSize.GCT.hetero.fitBefore["mse"] %>% mutate(Method = "GCT", n = n.seq, honest = "honest", hetero = "heterogeneous") 
subset.honest.CT.hetero.honest     <- sampleSize.CT.hetero.honest["mse"] %>% mutate(Method = "CT", n = n.seq, honest = "honest", hetero = "heterogeneous") 
subset.honest.opt.ct.hetero.honest <- sampleSize.opt.ct.hetero.honest["mse"] %>% mutate(Method = "opt.ct", n = n.seq, honest = "honest", hetero = "heterogeneous") 
subset.honest.GCT.homo.fitBefore   <- sampleSize.GCT.homo.fitBefore["mse"] %>% mutate(Method = "GCT", n = n.seq, honest = "honest", hetero = "homogeneous") 
subset.honest.CT.homo.honest       <- sampleSize.CT.homo.honest["mse"] %>% mutate(Method = "CT", n = n.seq, honest = "honest", hetero = "homogeneous") 
subset.honest.opt.ct.homo.honest   <- sampleSize.opt.ct.homo.honest["mse"] %>% mutate(Method = "opt.ct", n = n.seq, honest = "honest", hetero = "homogeneous") 

ggplot(data= sampleSize.GCT.homo.fitBefore, aes(x = c(10, 20, 30, 40, 50, 80, 100, 125), y = mse)) +
  geom_line(size = 0.5) +
  geom_point(size = 1.5) +
  scale_y_log10()
