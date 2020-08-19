# Check packages. 
#Install and load packages

if (!require(caret)) {install.packages("caret"); library(caret)} 
if (!require(causalTree)) {install_github("susanathey/causalTree"); library(causalTree)}
if (!require(devtools)) {install.packages("devtools"); library(devtools)}

#if (!require(causalForest)) {install_github("swager/causalForest"); library(causalForest)}
if (!require(dplyr)) {install.packages("dplyr"); library(dplyr)}
if (!require(grf)) {install.packages("grf"); library(grf)} 
if (!require(gbm)) {install.packages("gbm"); library(gbm)} 
if (!require(MASS)) {install.packages("MASS"); library(MASS)} 
if (!require(foreach)) {install.packages("foreach"); library(foreach)} 
if (!require(doParallel)) {install.packages("doParallel"); library(doParallel)} 
if (!require(rpart)) {install.packages("rpart"); library(rpart)} 
if (!require(randomForestSRC)) {install.packages("randomForestSRC"); library(randomForestSRC)}
if (!require(rattle)) {install.packages("rattle"); library(rattle)} 
if (!require(rpart.plot)) {install.packages("rpart.plot"); library(rpart.plot)} 
if (!require(tidyverse)) {install.packages("tidyverse"); library(tidyverse)}
if (!require(forcats)) {install.packages("forcats"); library(forcats)}
if (!require(RColorBrewer)) {install.packages("RcolorBrewer"); library(RcolorBrewer)}
