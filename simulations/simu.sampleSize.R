n.seq <- c(1000, 3000, 5000)
# n.seq <- c(1000, 3000, 5000)
# increase minsplit and minsize for growing sample size (minsize * 2 in CIT)
## n = 1000, 3000, 5000, minsplit = 30, minsize = 10
## n = 8000, 10000, 15000, 20000, minsplit = 60, minsize = 30
sim.r <- 100

# Start Cluster 
clnum <- detectCores()
cl <- makeCluster(getOption("cl.cores", clnum - 1))
registerDoParallel(cl)
pack <- c("caret", "devtools", "causalTree", "causalForest", "dplyr", "grf",
          "gbm", "MASS", "rpart", "randomForestSRC", "rattle", "rpart.plot")


############################################################################################################
##################################### 1.1. GCT: hetero, nonhonest ##########################################
############################################################################################################
tmp.GCT.hetero.fitBefore <- data.frame(matrix(0, length(n.seq), 10))
for (iter in 1:sim.r) {
  tmp <- foreach(i = 1:length(n.seq), .combine = rbind, .export = ls(), .packages = pack) %dopar% hetero.fitBefore.simu.func(nn = n.seq[i])
  tmp.GCT.hetero.fitBefore <- tmp.GCT.hetero.fitBefore + tmp
}

sampleSize.GCT.hetero.fitBefore <- tmp.GCT.hetero.fitBefore/sim.r
colnames(sampleSize.GCT.hetero.fitBefore) <- colnames(tmp)



############################################################################################################
##################################### 1.2. GCT: hetero, nonhonest ##########################################
############################################################################################################
tmp.GCT.hetero.fitinparent <- data.frame(matrix(0, length(n.seq), 10))
for (iter in 1:sim.r) {
  tmp <- foreach(i = 1:length(n.seq), .combine = rbind, .export = ls(), .packages = pack) %dopar% hetero.fitinparent.simu.func(nn = n.seq[i])
  tmp.GCT.hetero.fitinparent <- tmp.GCT.hetero.fitinparent + tmp
}

sampleSize.GCT.hetero.fitinparent <- tmp.GCT.hetero.fitinparent/sim.r
colnames(sampleSize.GCT.hetero.fitinparent) <- colnames(tmp)



############################################################################################################
##################################### 2.1. GCT: homo, nonhonest ############################################
############################################################################################################
tmp.GCT.homo.fitBefore <- data.frame(matrix(0, length(n.seq), 10))
for (iter in 1:sim.r) {
  tmp <- foreach(i = 1:length(n.seq), .combine = rbind, .export = ls(), .packages = pack) %dopar% homo.fitBefore.simu.func(nn = n.seq[i])
  tmp.GCT.homo.fitBefore <- tmp.GCT.homo.fitBefore + tmp
}

sampleSize.GCT.homo.fitBefore <- tmp.GCT.homo.fitBefore/sim.r
colnames(sampleSize.GCT.homo.fitBefore) <- colnames(tmp)



############################################################################################################
##################################### 2.2. GCT: homo, nonhonest ############################################
############################################################################################################
tmp.GCT.homo.fitinparent <- data.frame(matrix(0, length(n.seq), 10))
for (iter in 1:sim.r) {
  tmp <- foreach(i = 1:length(n.seq), .combine = rbind, .export = ls(), .packages = pack) %dopar% homo.fitinparent.simu.func(nn = n.seq[i])
  tmp.GCT.homo.fitinparent <- tmp.GCT.homo.fitinparent + tmp
}

sampleSize.GCT.homo.fitinparent <- tmp.GCT.homo.fitinparent/sim.r
colnames(sampleSize.GCT.homo.fitinparent) <- colnames(tmp)



############################################################################################################
######################################## 3. CT: hetero, nonhonest ##########################################
############################################################################################################
tmp.CT.hetero.nonhonest <- data.frame(matrix(0, length(n.seq), 13))
for (iter in 1:sim.r) {
  tmp <- foreach(i = 1:length(n.seq), .combine = rbind, .export = ls(), .packages = pack) %dopar%  CT.hetero.nonhonest.simu.func(nn = n.seq[i])
  tmp.CT.hetero.nonhonest <- tmp.CT.hetero.nonhonest + tmp
}

sampleSize.CT.hetero.nonhonest <- tmp.CT.hetero.nonhonest / sim.r
colnames(sampleSize.CT.hetero.nonhonest) <- colnames(tmp)


############################################################################################################
######################################### 4. CT: homo, nonhonest ###########################################
############################################################################################################
tmp.CT.homo.nonhonest <- data.frame(matrix(0, length(n.seq), 13))
for (iter in 1:sim.r) {
  tmp <- foreach(i = 1:length(n.seq), .combine = rbind, .export = ls(), .packages = pack) %dopar%  CT.homo.nonhonest.simu.func(nn = n.seq[i])
  tmp.CT.homo.nonhonest <- tmp.CT.homo.nonhonest + tmp
}

sampleSize.CT.homo.nonhonest <- tmp.CT.homo.nonhonest / sim.r
colnames(sampleSize.CT.homo.nonhonest) <- colnames(tmp)



############################################################################################################
###################################### 5. opt.CT: hetero, nonhonest ########################################
############################################################################################################
tmp.opt.ct.hetero.nonhonest <- data.frame(matrix(0, length(n.seq), 13))
for (iter in 1:sim.r) {
  tmp <- foreach(i = 1:length(n.seq), .combine = rbind, .export = ls(), .packages = pack) %dopar% opt.ct.hetero.nonhonest.simu.func(nn = n.seq[i])
  tmp.opt.ct.hetero.nonhonest <- tmp.opt.ct.hetero.nonhonest + tmp
}

sampleSize.opt.ct.hetero.nonhonest <- tmp.opt.ct.hetero.nonhonest / sim.r
colnames(sampleSize.opt.ct.hetero.nonhonest) <- colnames(tmp)


# ²¹³äÊµÑé
# tmp.opt.ct.add <- data.frame(matrix(0, 2, 13))
# for (iter in 1:sim.r) {
#   tmp <- foreach(i = 7:8, .combine = rbind, .export = ls(), .packages = pack) %dopar% opt.ct.hetero.nonhonest.simu.func(nn = n.seq[i])
#   tmp.opt.ct.add <- tmp.opt.ct.add + tmp
# }
# 
# opt.ct.add <- tmp.opt.ct.add / sim.r
# colnames(opt.ct.add) <- colnames(tmp)



############################################################################################################
######################################## 6. opt.CT: homo, nonhonest ########################################
############################################################################################################
tmp.opt.ct.homo.nonhonest <- data.frame(matrix(0, length(n.seq), 13))
for (iter in 1:sim.r) {
  tmp <- foreach(i = 1:length(n.seq), .combine = rbind, .export = ls(), .packages = pack) %dopar%  opt.ct.homo.nonhonest.simu.func(nn = n.seq[i])
  tmp.opt.ct.homo.nonhonest <- tmp.opt.ct.homo.nonhonest + tmp
}

sampleSize.opt.ct.homo.nonhonest <- tmp.opt.ct.homo.nonhonest / sim.r
colnames(sampleSize.opt.ct.homo.nonhonest) <- colnames(tmp)



############################################################################################################
###################################### 7.1. CIT: hetero, fitBefore #########################################
############################################################################################################
tmp.CIT.DR.hetero.fitBefore <- data.frame(matrix(0, length(n.seq), 13))
for (iter in 1:sim.r) {
  tmp <- foreach(i = 1:length(n.seq), .combine = rbind, .export = ls(), .packages = pack) %dopar% CIT.DR.hetero.fitBefore.simu.func(nn = n.seq[i])
  tmp.CIT.DR.hetero.fitBefore <- tmp.CIT.DR.hetero.fitBefore + tmp
}

sampleSize.CIT.DR.hetero.fitBefore <- tmp.CIT.DR.hetero.fitBefore/sim.r
colnames(sampleSize.CIT.DR.hetero.fitBefore) <- colnames(tmp)



############################################################################################################
###################################### 7.2. CIT: hetero, fitinparent #######################################
############################################################################################################
tmp.CIT.DR.hetero.fitinparent <- data.frame(matrix(0, length(n.seq), 13))
for (iter in 1:sim.r) {
  tmp <- foreach(i = 1:length(n.seq), .combine = rbind, .export = ls(), .packages = pack) %dopar% CIT.DR.hetero.fitinparent.simu.func(nn = n.seq[i])
  tmp.CIT.DR.hetero.fitinparent <- tmp.CIT.DR.hetero.fitinparent + tmp
}

sampleSize.CIT.DR.hetero.fitinparent <- tmp.CIT.DR.hetero.fitinparent/sim.r
colnames(sampleSize.CIT.DR.hetero.fitinparent) <- colnames(tmp)



############################################################################################################
####################################### 8.1. CIT: homo, fitBefore ##########################################
############################################################################################################
tmp.CIT.DR.homo.fitBefore <- data.frame(matrix(0, length(n.seq), 13))
for (iter in 1:sim.r) {
  tmp <- foreach(i = 1:length(n.seq), .combine = rbind, .export = ls(), .packages = pack) %dopar% CIT.DR.homo.fitBefore.simu.func(nn = n.seq[i])
  tmp.CIT.DR.homo.fitBefore <- tmp.CIT.DR.homo.fitBefore + tmp
}

sampleSize.CIT.DR.homo.fitBefore <- tmp.CIT.DR.homo.fitBefore/sim.r
colnames(sampleSize.CIT.DR.homo.fitBefore) <- colnames(tmp)  



############################################################################################################
####################################### 8.2. CIT: homo, fitinparent ########################################
############################################################################################################
tmp.CIT.DR.homo.fitinparent <- data.frame(matrix(0, length(n.seq), 13))
for (iter in 1:sim.r) {
  tmp <- foreach(i = 1:length(n.seq), .combine = rbind, .export = ls(), .packages = pack) %dopar% CIT.DR.homo.fitinparent.simu.func(nn = n.seq[i])
  tmp.CIT.DR.homo.fitinparent <- tmp.CIT.DR.homo.fitinparent + tmp
}

sampleSize.CIT.DR.homo.fitinparent <- tmp.CIT.DR.homo.fitinparent/sim.r
colnames(sampleSize.CIT.DR.homo.fitinparent) <- colnames(tmp) 


# close the clusters
stopImplicitCluster()

save(sampleSize.GCT.hetero.fitBefore, sampleSize.CT.hetero.nonhonest, 
     sampleSize.opt.ct.hetero.nonhonest, sampleSize.CIT.DR.hetero.fitBefore,
     sampleSize.GCT.homo.fitBefore, sampleSize.CT.homo.nonhonest, 
     sampleSize.opt.ct.homo.nonhonest, sampleSize.CIT.DR.homo.fitBefore,
     sampleSize.GCT.hetero.fitinparent, sampleSize.GCT.homo.fitinparent,
     sampleSize.CIT.DR.hetero.fitinparent, sampleSize.CIT.DR.homo.fitinparent,
     file = "E:/Project/groupTest/MethodCompare/Simulation/sampleSize.RData")


save(b.sampleSize.GCT.hetero.fitBefore, 
     b.sampleSize.GCT.hetero.fitinparent, 
     b.sampleSize.GCT.homo.fitBefore, 
     b.sampleSize.GCT.homo.fitinparent,
     b.sampleSize.CIT.DR.hetero.fitBefore, 
     b.sampleSize.CIT.DR.hetero.fitinparent,
     b.sampleSize.CIT.DR.homo.fitBefore, 
     b.sampleSize.CIT.DR.homo.fitinparent,
     b.sampleSize.CT.hetero.nonhonest, 
     b.sampleSize.CT.homo.nonhonest,
     b.sampleSize.opt.ct.hetero.nonhonest, 
     b.sampleSize.opt.ct.homo.nonhonest,
     file = "E://Adsimu//ad_sampleSize_5.RData")

# save(add.sampleSize.opt.ct.homo.nonhonest, 
#      add.sampleSize.opt.ct.hetero.nonhonest,
#      add.sampleSize.CT.hetero.nonhonest,
#      add.sampleSize.CT.homo.nonhonest,
#      file = "E://Adsimu//ad_sampleSize_2.RData")

b.sampleSize.GCT.hetero.fitBefore   <- sampleSize.GCT.hetero.fitBefore
b.sampleSize.GCT.hetero.fitinparent <- sampleSize.GCT.hetero.fitinparent
b.sampleSize.GCT.homo.fitBefore      <- sampleSize.GCT.homo.fitBefore
b.sampleSize.GCT.homo.fitinparent    <- sampleSize.GCT.homo.fitinparent
b.sampleSize.CIT.DR.hetero.fitBefore    <- sampleSize.CIT.DR.hetero.fitBefore
b.sampleSize.CIT.DR.hetero.fitinparent  <- sampleSize.CIT.DR.hetero.fitinparent
b.sampleSize.CIT.DR.homo.fitBefore     <- sampleSize.CIT.DR.homo.fitBefore
b.sampleSize.CIT.DR.homo.fitinparent   <- sampleSize.CIT.DR.homo.fitinparent
b.sampleSize.CT.hetero.nonhonest   <- sampleSize.CT.hetero.nonhonest
b.sampleSize.CT.homo.nonhonest     <- sampleSize.CT.homo.nonhonest
b.sampleSize.opt.ct.hetero.nonhonest  <- sampleSize.opt.ct.hetero.nonhonest
b.sampleSize.opt.ct.homo.nonhonest    <- sampleSize.opt.ct.homo.nonhonest
# 
# save(sampleSize.GCT.hetero.fitBefore.2,
#      sampleSize.GCT.hetero.fitinparent.2, 
#      sampleSize.GCT.homo.fitBefore.2, 
#      sampleSize.GCT.homo.fitinparent.2, 
#      sampleSize.CIT.DR.hetero.fitBefore.2,   
#      sampleSize.CIT.DR.hetero.fitinparent.2, 
#      sampleSize.CIT.DR.homo.fitBefore.2,
#      sampleSize.CIT.DR.homo.fitinparent.2,
#      sampleSize.CT.hetero.nonhonest.2,
#      sampleSize.CT.homo.nonhonest.2,
#      sampleSize.opt.ct.hetero.nonhonest.2,
#      sampleSize.opt.ct.homo.nonhonest.2, 
#      file = "E://Adsimu//ad_sampleSize3.RData")


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
