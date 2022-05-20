# organize the results
Output <- function(origin.res, hetero = TRUE){
  num.corr.tree <- sum(origin.res[,"corr.tree"])
  if(hetero){
    output <- origin.res %>% as.data.frame() %>%
      summarize(corr.tree.prop = mean(corr.tree, na.rm = TRUE),
                num.leaf = mean(num.leaf, na.rm = TRUE),
                num.noise = mean(num.noise.splt, na.rm = TRUE),
                fir.splt.acc = mean(fir.splt.ind, na.rm = TRUE),
                pps = mean(pps, na.rm = TRUE),
                abs.split.err = mean(abs.err.sp, na.rm = TRUE),
                mse = mean(mse, na.rm = TRUE), 
                cvr_1 = mean(ifelse(corr.tree == 1, cvr1, NA), na.rm = TRUE),#sum(cvr1)/num.corr.tree,
                cvr_2 = mean(ifelse(corr.tree == 1, cvr2, NA), na.rm = TRUE),#sum(cvr2)/num.corr.tree,
                cvr_sim = mean(ifelse(corr.tree == 1, sim_cvr, NA), na.rm = TRUE),#sum(sim_cvr)/num.corr.tree,
                cov_ratio = mean(cov_ratio, na.rm = TRUE), 
                widthCI1 = mean(ifelse(corr.tree == 1, widthCI1, NA), na.rm = TRUE),
                widthCI2 = mean(ifelse(corr.tree == 1, widthCI2, NA), na.rm = TRUE),
                time = mean(time, na.rm = TRUE))
  } else {
    output <- origin.res %>% as.data.frame() %>%
      summarize(corr.tree.prop = mean(corr.tree, na.rm = TRUE),
                num.leaf = mean(num.leaf, na.rm = TRUE),
                num.noise = mean(num.noise.splt, na.rm = TRUE),
                fir.splt.acc = mean(fir.splt.ind, na.rm = TRUE),
                pps = mean(pps, na.rm = TRUE),
                abs.split.err = mean(abs.err.sp, na.rm = TRUE),
                mse = mean(mse, na.rm = TRUE), 
                cvr = mean(ifelse(corr.tree == 1, cvr, NA), na.rm = TRUE), #sum(cvr)/num.corr.tree,
                cov_ratio = mean(cov_ratio, na.rm = TRUE), 
                widthCI = mean(ifelse(corr.tree == 1, widthCI, NA), na.rm = TRUE),
                time = mean(time, na.rm = TRUE))
  }
  return(output)
}
