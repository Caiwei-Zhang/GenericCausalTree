if (!require(forcats)) {install.packages("forcats"); library(forcats)}
if (!require(tidyverse)) {install.packages("tidyverse"); library(tidyverse)}
if (!require(ggplot2)) {intsall.packages("ggplot2"); library(ggplot2)}

setwd("..")
load("TEDSA_2015_PUF.Rda") 
# head(tedsa_2015, n = 10)
Maryland <- tedsa_2015 %>% dplyr::filter(STFIPS %in% c(24))

# colnames(Maryland)
Maryland_select <- Maryland %>% 
  dplyr::select(AGE, GENDER, RACE, MARSTAT, EDUC, EMPLOY, PREG, VET, DETNLF,
                LIVARAG, PRIMINC, ARRESTS, SERVSETA, METHUSE, DAYWAIT, 
                PSOURCE, DETCRIM, NOPRIOR, SUB1, ROUTE1, FREQ1,
                FRSTUSE1, NUMSUBS, PSYPROB, HLTHINS, PRIMPAY)
colnames(Maryland_select) <- str_to_lower(colnames(Maryland_select))
  
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

# merge levels of covariate age and frstuse1 (based on the distribution of variable and daywait)
ggplot(data = Maryland_filter) + geom_bar(mapping = aes(x = as.factor(age)))
ggplot(data = Maryland_filter) + geom_bar(mapping = aes(x = as.factor(frstuse1)))

Maryland_merge <- tibble(X.tmp, race = Maryland_filter$race, daywait = Maryland_filter$daywait) %>% 
  dplyr::mutate(age = fct_collapse(age,
                                   "young" = as.character(c(2, 3, 4, 5, 6)),
                                   "middle" = as.character(c(7, 8, 9)),
                                   "old" = as.character(c(10, 11, 12)))) %>% 
  dplyr::mutate(frstuse1 = fct_collapse(frstuse1, 
                                        "1" = as.character(c(1, 2, 3)),
                                        "2" = as.character(c(4, 5)),
                                        "3" = as.character(c(6, 7, 8, 9, 10, 11, 12))))                      

ggplot(data = Maryland_merge) + geom_bar(mapping = aes(x = age))
ggplot(data = Maryland_merge) + geom_bar(mapping = aes(x = frstuse1))
ggplot(data = Maryland_merge) + geom_density(mapping = aes(x = daywait)) + facet_grid(age ~.)
ggplot(data = Maryland_merge) + geom_density(mapping = aes(x = daywait)) + facet_grid(frstuse1 ~.)

Maryland_merge$race    <- ifelse(Maryland_merge$race == 4, 0, 1) # 5-white-1, 4-black-0

# Maryland_untrunc <- as.data.frame(Maryland_merge) %>% rename(., A = race, Y = daywait)
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
fancyRpartPlot(opt_honest_gct_tedsa_trunc, type = 4, palettes = "Greens", 
               sub = "honest optimal tree for TEDS-A (Maryland, 2015)", cex = 0.6)

save(gct_tedsa_untrunc, gct_tedsa_trunc, honest_gct_tedsa_trunc,
     opt_tedsa_untrunc, opt_tedsa_trunc, opt_honest_gct_tedsa_trunc,
     file = "..\\data_analysis_res.RData")

##### plot: honest
for (g in 1:opt_honest_gct_tedsa_trunc$num.leaf) {
  index <- names(opt_tedsa_trunc$where)[opt_tedsa_trunc$where == which(opt_tedsa_trunc$frame$var == "<leaf>")[g]]
  index <- as.numeric(index)
  subgroup[index] <- g
}

pred <- predict(opt_honest_gct_tedsa_trunc, Maryland_teds)

y_dens <- data.frame(daywait= Maryland_teds$Y, Race = Maryland_teds$A, STE = round(pred, 3))
y_dens$Race <- factor(y_dens$Race, levels = c(0, 1), labels = c("African American", "Whites"))
y_dens$STE <- factor(y_dens$STE, levels = unique(y_dens$STE), 
                           labels = paste("ATE:", unique(y_dens$STE), sep = ""))

honest_gct_density <- ggplot(data = y_dens, aes(x = daywait, fill = STE, linetype = STE)) +
                        geom_density(alpha = 0.5) +
                        facet_grid(Race~STE) +
                        theme(axis.title.x=element_text(size=10),
                              axis.title.y=element_text(size=10),
                              axis.text.x = element_text(size = 8),
                              legend.position = "bottom",
                              legend.key = element_rect(fill = "white"),
                              legend.key.size = unit(12, "pt"))


save(honest_gct_tedsa_trunc, opt_honest_gct_tedsa_trunc, y_dens,
     file = "\honest_gct_res.RData")
