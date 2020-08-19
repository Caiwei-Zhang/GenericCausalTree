set.seed(999)
dim.seq <- c(10, 15, 20)
sim.r <- 1000

############################################################################################################
####################################### 1.1. GCT: hetero, nonhonest ########################################
############################################################################################################
tmp.noise.GCT.hetero.fitBefore <- data.frame(matrix(0, length(dim.seq), 10))
for (iter in 1:sim.r) {
  tmp <- foreach(i = 1:length(dim.seq), .combine = rbind, .export = ls(), .packages = pack) %dopar% 
    hetero.fitBefore.simu.func(nn = 1000, dim = dim.seq[i])
  tmp.noise.GCT.hetero.fitBefore <- tmp.noise.GCT.hetero.fitBefore + tmp
}

noise.GCT.hetero.fitBefore <- tmp.noise.GCT.hetero.fitBefore/sim.r
colnames(noise.GCT.hetero.fitBefore) <- colnames(tmp)


############################################################################################################
##################################### 1.2. GCT: hetero, nonhonest ##########################################
############################################################################################################
tmp.noise.GCT.hetero.fitinparent <- data.frame(matrix(0, length(dim.seq), 10))
for (iter in 1:sim.r) {
  tmp <- foreach(i = 1:length(dim.seq), .combine = rbind, .export = ls(), .packages = pack) %dopar% 
    hetero.fitinparent.simu.func(nn = 1000, dim = dim.seq[i])
  tmp.noise.GCT.hetero.fitinparent <- tmp.noise.GCT.hetero.fitinparent + tmp
}

noise.GCT.hetero.fitinparent <- tmp.noise.GCT.hetero.fitinparent/sim.r
colnames(noise.GCT.hetero.fitinparent) <- colnames(tmp)


############################################################################################################
####################################### 2. GCT: homo, nonhonest ############################################
############################################################################################################
tmp.noise.GCT.homo.fitBefore <- data.frame(matrix(0, length(dim.seq), 10))
for (iter in 1:sim.r) {
  tmp <- foreach(i = 1:length(dim.seq), .combine = rbind, .export = ls(), .packages = pack) %dopar% 
    homo.fitBefore.simu.func(nn = 1000, dim = dim.seq[i])
  tmp.noise.GCT.homo.fitBefore <- tmp.noise.GCT.homo.fitBefore + tmp
}

noise.GCT.homo.fitBefore <- tmp.noise.GCT.homo.fitBefore/sim.r
colnames(noise.GCT.homo.fitBefore) <- colnames(tmp)


############################################################################################################
##################################### 2.2. GCT: homo, nonhonest ##########################################
############################################################################################################
tmp.noise.GCT.homo.fitinparent <- data.frame(matrix(0, length(dim.seq), 10))
for (iter in 1:sim.r) {
  tmp <- foreach(i = 1:length(dim.seq), .combine = rbind, .export = ls(), .packages = pack) %dopar% 
    homo.fitinparent.simu.func(nn = 1000, dim = dim.seq[i])
  tmp.noise.GCT.homo.fitinparent <- tmp.noise.GCT.homo.fitinparent + tmp
}

noise.GCT.homo.fitinparent <- tmp.noise.GCT.homo.fitinparent/sim.r
colnames(noise.GCT.homo.fitinparent) <- colnames(tmp)


############################################################################################################
######################################## 3. CT: hetero, nonhonest ##########################################
############################################################################################################
tmp.noise.CT.hetero.nonhonest <- data.frame(matrix(0, length(dim.seq), 13))
for (iter in 1:sim.r) {
  tmp <- foreach(i = 1:length(dim.seq), .combine = rbind, .export = ls(), .packages = pack) %dopar%  
    CT.hetero.nonhonest.simu.func(nn = 1000, dim = dim.seq[i])
  tmp.noise.CT.hetero.nonhonest <- tmp.noise.CT.hetero.nonhonest + tmp
}

noise.CT.hetero.nonhonest <- tmp.noise.CT.hetero.nonhonest / sim.r
colnames(noise.CT.hetero.nonhonest) <- colnames(tmp)


############################################################################################################
######################################### 4. CT: homo, nonhonest ###########################################
############################################################################################################
tmp.noise.CT.homo.nonhonest <- data.frame(matrix(0, length(dim.seq), 13))
for (iter in 1:sim.r) {
  tmp <- foreach(i = 1:length(dim.seq), .combine = rbind, .export = ls(), .packages = pack) %dopar%  
    CT.hetero.nonhonest.simu.func(nn = 1000, dim = dim.seq[i])
  tmp.noise.CT.homo.nonhonest <- tmp.noise.CT.homo.nonhonest + tmp
}

noise.CT.homo.nonhonest <- tmp.noise.CT.homo.nonhonest / sim.r
colnames(noise.CT.homo.nonhonest) <- colnames(tmp)



############################################################################################################
###################################### 5. opt.CT: hetero, nonhonest ########################################
############################################################################################################
tmp.noise.opt.ct.hetero.nonhonest <- data.frame(matrix(0, length(dim.seq), 13))
for (iter in 1:sim.r) {
  tmp <- foreach(i = 1:length(dim.seq), .combine = rbind, .export = ls(), .packages = pack) %dopar% 
    opt.ct.hetero.nonhonest.simu.func(nn = 1000, dim = dim.seq[i])
  tmp.noise.opt.ct.hetero.nonhonest <- tmp.noise.opt.ct.hetero.nonhonest + tmp
}

noise.opt.ct.hetero.nonhonest <- tmp.noise.opt.ct.hetero.nonhonest / sim.r
colnames(noise.opt.ct.hetero.nonhonest) <- colnames(tmp)


############################################################################################################
######################################## 6. opt.CT: homo, nonhonest ########################################
############################################################################################################
tmp.noise.opt.ct.homo.nonhonest <- data.frame(matrix(0, length(dim.seq), 13))
for (iter in 1:sim.r) {
  tmp <- foreach(i = 1:length(dim.seq), .combine = rbind, .export = ls(), .packages = pack) %dopar%  
    opt.ct.homo.nonhonest.simu.func(nn = 1000, dim = dim.seq[i])
  tmp.noise.opt.ct.homo.nonhonest <- tmp.noise.opt.ct.homo.nonhonest + tmp
}

noise.opt.ct.homo.nonhonest <- tmp.noise.opt.ct.homo.nonhonest / sim.r
colnames(noise.opt.ct.homo.nonhonest) <- colnames(tmp)


############################################################################################################
######################################## 7.1. CIT: hetero, nonhonest #######################################
############################################################################################################
tmp.noise.CIT.DR.hetero.fitBefore <- data.frame(matrix(0, length(dim.seq), 13))
for (iter in 1:sim.r) {
  tmp <- foreach(i = 1:length(dim.seq), .combine = rbind, .export = ls(), .packages = pack) %dopar% 
    CIT.DR.hetero.fitBefore.simu.func(nn = 1000, dim = dim.seq[i])
  tmp.noise.CIT.DR.hetero.fitBefore <- tmp.noise.CIT.DR.hetero.fitBefore + tmp
}

noise.CIT.DR.hetero.fitBefore <- tmp.noise.CIT.DR.hetero.fitBefore/sim.r
colnames(noise.CIT.DR.hetero.fitBefore) <- colnames(tmp)


############################################################################################################
######################################## 7.2. CIT: hetero, nonhonest #######################################
############################################################################################################
tmp.noise.CIT.DR.hetero.fitinparent <- data.frame(matrix(0, length(dim.seq), 13))
for (iter in 1:sim.r) {
  tmp <- foreach(i = 1:length(dim.seq), .combine = rbind, .export = ls(), .packages = pack) %dopar% 
    CIT.DR.hetero.fitinparent.simu.func(nn = 1000, dim = dim.seq[i])
  tmp.noise.CIT.DR.hetero.fitinparent <- tmp.noise.CIT.DR.hetero.fitinparent + tmp
}

noise.CIT.DR.hetero.fitinparent <- tmp.noise.CIT.DR.hetero.fitinparent/sim.r
colnames(noise.CIT.DR.hetero.fitinparent) <- colnames(tmp)


############################################################################################################
######################################### 8.1 CIT: homo, nonhonest #########################################
############################################################################################################
tmp.noise.CIT.DR.homo.fitBefore <- data.frame(matrix(0, length(dim.seq), 13))
for (iter in 1:sim.r) {
  tmp <- foreach(i = 1:length(dim.seq), .combine = rbind, .export = ls(), .packages = pack) %dopar% 
    CIT.DR.hetero.fitBefore.simu.func(nn = 1000, dim = dim.seq[i])
  tmp.noise.CIT.DR.homo.fitBefore <- tmp.noise.CIT.DR.homo.fitBefore + tmp
}

noise.CIT.DR.homo.fitBefore <- tmp.noise.CIT.DR.homo.fitBefore/sim.r
colnames(noise.CIT.DR.homo.fitBefore) <- colnames(tmp)  

############################################################################################################
######################################### 8.2 CIT: homo, nonhonest #########################################
############################################################################################################
tmp.noise.CIT.DR.homo.fitinparent <- data.frame(matrix(0, length(dim.seq), 13))
for (iter in 1:sim.r) {
  tmp <- foreach(i = 1:length(dim.seq), .combine = rbind, .export = ls(), .packages = pack) %dopar% 
    CIT.DR.homo.fitinparent.simu.func(nn = 1000, dim = dim.seq[i])
  tmp.noise.CIT.DR.homo.fitinparent <- tmp.noise.CIT.DR.homo.fitinparent + tmp
}

noise.CIT.DR.homo.fitinparent <- tmp.noise.CIT.DR.homo.fitinparent/sim.r
colnames(noise.CIT.DR.homo.fitinparent) <- colnames(tmp) 

# noise.GCT <- list(hetero = noise.GCT.hetero.fitBefore, homo = noise.GCT.homo.fitBefore)
# noise.CT <- list(hetero = noise.CT.hetero.nonhonest, homo = noise.CT.homo.nonhonest )
# noise.opt.CT <- list(hetero = noise.opt.ct.hetero.nonhonest, homo = noise.opt.ct.homo.nonhonest)
# noise.CIT.DR <- list(hetero = noise.CIT.DR.hetero.fitBefore, homo = noise.CIT.DR.homo.fitBefore)

save(noise.GCT.hetero.fitBefore, noise.GCT.homo.fitBefore, 
     noise.CT.hetero.nonhonest, noise.CT.homo.nonhonest,
     noise.opt.ct.hetero.nonhonest, noise.opt.ct.homo.nonhonest, 
     noise.CIT.DR.hetero.fitBefore, noise.CIT.DR.homo.fitBefore,
     file = "E:/project/groupTest/MethodCompare/Simulation/noiseCov.RData")



