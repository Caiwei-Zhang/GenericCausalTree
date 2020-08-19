if (!require(forcats)) {install.packages("forcats"); library(forcats)}
if (!require(tidyverse)) {install.packages("tidyverse"); library(tidyverse)}

load("E:\\GenericCausalTree\\Data_analysis\\TEDSA_2015_PUF.Rda") 
# load("C:/Users/ZCW/AppData/Local/Temp/Temp1_TEDS-A-2015-DS0001-bndl-data-r.zip/TEDSA_2015_PUF.Rda")
# View(tedsa_2015)
Maryland <- tedsa_2015 %>% dplyr::filter(STFIPS %in% c(24))
# colnames(Maryland)
Maryland_select <- Maryland %>% 
  dplyr::select(AGE, GENDER, RACE, MARSTAT, EDUC, EMPLOY, PREG, VET, DETNLF,
                LIVARAG, PRIMINC, ARRESTS, SERVSETA, METHUSE, DAYWAIT, 
                PSOURCE, DETCRIM, NOPRIOR, SUB1, ROUTE1, FREQ1,
                FRSTUSE1, NUMSUBS, PSYPROB, HLTHINS, PRIMPAY) %>%
  rename(., age = AGE, gender = GENDER, race = RACE, marstat = MARSTAT, 
         educ = EDUC, employ = EMPLOY,  preg = PREG, vet = VET, detnlf = DETNLF,
         livarag = LIVARAG, priminc = PRIMINC, arrests = ARRESTS, 
         servseta = SERVSETA, methuse = METHUSE, daywait = DAYWAIT, 
         psource = PSOURCE, detcrim = DETCRIM, noprior = NOPRIOR,
         sub1 = SUB1, route1 = ROUTE1, freq1 = FREQ1, 
         frstuse1 = FRSTUSE1, numsubs = NUMSUBS, psyprob = PSYPROB, 
         hlthins = HLTHINS, primpay = PRIMPAY) 

Maryland_filter <- Maryland_select %>%
  filter(gender %in% c(1,2)) %>%
  filter(race %in% c(4, 5)) %>% 
  filter(marstat %in% c(1, 2, 3, 4)) %>%
  filter(employ %in% c(1, 2, 3, 4)) %>%
  # filter(detnlf %in% c(1, 2, 3, 5, 6, -9)) %>%
  filter(vet %in% c(1, 2)) %>%
  filter(livarag %in% c(1, 2, 3)) %>%
  # filter(priminc %in% c(1, 2, 3, 20, -9)) %>%
  filter(arrests %in% c(0, 1, 2)) %>%
  # filter(servseta %in% c(4, 6, 7)) %>%
  # filter(detcrim %in% c(1, 3, 7, 8, -9)) %>% 
  filter(methuse %in% c(1, 2)) %>%
  filter(psource %in% c(1, 2, 3, 4, 5, 6, 7)) %>%
  filter(noprior %in% c(0, 1, 2, 3, 4, 5)) %>%
  filter(sub1 %in% c(5, 7)) %>%
  filter(route1 %in% c(1, 2, 3, 4, 20)) %>%
  filter(freq1 %in% c(1, 2, 3, 4, 5)) %>% 
  filter(frstuse1 %in% c(1:12)) %>%
  filter(numsubs %in% c(1, 2, 3)) %>% 
  filter(psyprob %in% c(1, 2)) %>%
  filter(hlthins %in% c(1, 2, 3, 4)) %>%
  filter(primpay %in% c(1, 2, 3, 4, 5, 9)) %>%
  filter(!daywait %in% c(-9)) %>% 
  dplyr::select(daywait, race, everything())

# remove extra levels of factor covariates and transform categorical covariates as fators
X.tmp <- as_tibble(lapply(Maryland_filter %>% dplyr::select(age:primpay), factor))

# merge levels of covariate age and frstuse1
Maryland_merge <- tibble(X.tmp, race = Maryland_filter$race, daywait = Maryland_filter$daywait) %>% 
  dplyr::mutate(age = fct_collapse(age,
                                   "young" = as.character(c(2, 3, 4, 5)),
                                   "middle" = as.character(c(6, 7, 8, 9)),
                                   "old" = as.character(c(10, 11, 12)))) %>% 
  dplyr::mutate(frstuse1 = fct_collapse(frstuse1, 
                                        "1" = as.character(c(1, 2, 3, 4)),
                                        "2" = as.character(c(5, 6)),
                                        "3" = as.character(c(7, 8)),
                                        "4" = as.character(c(9, 10)),
                                        "5" = as.character(c(11, 12)))) 


Maryland_merge$race    <- ifelse(Maryland_merge$race == 4, 0, 1) # 5-white-1, 4-black-0

# truncate daywait at 99% quantitle of daywait values
quan <- quantile(Maryland_merge$daywait, 0.99)
Maryland_merge$daywait[Maryland_merge$daywait > quan] <- quan

# new data set 
Maryland_teds <- as.data.frame(Maryland_merge) %>% rename(., A = race, Y = daywait)



#############################################################################################
########################## Tree 1: untruncate daywait, nonhonest ############################
#############################################################################################
gct_tedsa_untrunc <- Treebuild(data = Maryland_teds, form.outcome = "Y ~ .", form.pscore = "A ~ .", 
                               type.outcome = "continuous", train.frac = 0.8, 
                               modelFitBefore = TRUE, outcome.mthd = "GBM", pscore.mthd = "GBM", 
                               cp = qchisq(0.95, 1), minsplit = 100, minsize = 100)

gct_tedsa_untrunc$res.tab
opt_tedsa_untrunc <- gct_tedsa_untrunc$optree
fancyRpartPlot(opt_tedsa_untrunc, type = 4, palettes = "Set2", 
               sub = "optimal tree for TEDS-A (Maryland, 2015)")


#############################################################################################
############################ Tree 2: truncate daywait, nonhonest ############################
#############################################################################################
gct_tedsa_trunc <- Treebuild(data = Maryland_teds, form.outcome = "Y ~ .", form.pscore = "A ~ .", 
                                type.outcome = "continuous", train.frac = 0.8, 
                                modelFitBefore = TRUE, outcome.mthd = "GBM", pscore.mthd = "GBM", 
                                cp = qchisq(0.95, 1), minsplit = 100, minsize = 100)

gct_tedsa_trunc$res.tab
opt_tedsa_trunc <- gct_tedsa_trunc$optree
fancyRpartPlot(opt_tedsa_trunc, type = 4, palettes = "Greens",
               sub = "optimal tree for TEDS-A (Maryland, 2015)")


#############################################################################################
############################# Tree 3: truncate daywait, honest ##############################
#############################################################################################
# honest tree building
honest_gct_tedsa_trunc <- honest.Treebuild(data = Maryland_teds, form.outcome = "Y ~ .", form.pscore = "A ~ .", 
                                             type.outcome = "continuous", train.frac = 0.8,
                                             modelFitBefore = TRUE, est.modelFitBefore = TRUE,
                                             outcome.mthd = "GBM", pscore.mthd = "GBM", cp = qchisq(0.95, 1), 
                                             minsplit = 100, minsize = 100)
honest_gct_tedsa_trunc$res.tab
opt_honest_gct_tedsa_trunc <- honest_gct_tedsa_trunc$optree
subopt_honest_gct_tedsa_trunc <- honest_gct_tedsa_trunc$subtree.list[[34]] # &33
fancyRpartPlot(opt_honest_gct_tedsa_trunc, type = 4, palettes = "Greens", 
               sub = "honest optimal tree for TEDS-A (Maryland, 2015)")

fancyRpartPlot(subopt_honest_gct_tedsa_trunc, type = 4, palettes = "Greens",
               sub = "honest suboptimal tree for TEDS-A (Maryland, 2015)")


save(gct_tedsa_untrunc, gct_tedsa_trunc, honest_gct_tedsa_trunc,
     opt_tedsa_untrunc, opt_tedsa_trunc, opt_honest_gct_tedsa_trunc,
     file = "E:\\Project\\groupTest\\MethodCompare\\Data_analysis\\data_analysis_res.RData ")




##### plot 
for (g in 1:opt_tedsa_trunc$num.leaf) {
  index <- names(opt_tedsa_trunc$where)[opt_tedsa_trunc$where == which(opt_tedsa_trunc$frame$var == "<leaf>")[g]]
  index <- as.numeric(index)
  subgroup[index] <- g
}

pred <- predict(opt_tedsa_trunc, Maryland_teds)

y_dens <- data.frame(Y= Maryland_teds$Y, A = Maryland_teds$A, predicted = round(pred, 3))
y_dens$A <- factor(y_dens$A, levels = c(0, 1), labels = c("African American (A = 0)", "Whites (A = 1)"))
y_dens$predicted <- factor(y_dens$predicted, levels = unique(y_dens$predicted), 
                           labels = paste("ATE:", unique(y_dens$predicted), sep = ""))

ggplot(data = y_dens, aes(x = Y, fill = predicted, linetype = predicted)) +
  geom_density(alpha = 0.5) +
  facet_grid(predicted~A)

# ggplot(data = y_dens, aes(x = Y, fill = predicted, linetype = predicted), size = 1) +
#   geom_density(alpha = .3)

# ggplot(data = y_a1_dens, mapping = aes(x = predicted, y = Y)) + 
#   geom_boxplot(fill = 'steelblue', outlier.colour = 'red', outlier.shape = 15)
# 
# ggplot(data = y_a0_dens, mapping = aes(x = predicted, y = Y)) + 
#   geom_boxplot(fill = 'steelblue', outlier.colour = 'red', outlier.shape = 15)