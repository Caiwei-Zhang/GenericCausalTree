
## caret
gbmCtrl <- trainControl(method = "cv", number = 10, selectionFunction = "best")
gbmGrid <- expand.grid(n.trees = c(100,200,300,500), interaction.depth= c(1,2), 
                       n.minobsinnode = 5, shrinkage = c(0.001, 0.05, 0.1))
t0 <- Sys.time()
set.seed(999)
modell <- train(Y ~ X1 + X2 + X3 + X4 + X5, data = data[1:500,], method = "gbm", 
                preProcess = c('center', 'scale'), trControl = gbmCtrl, tuneGrid = gbmGrid)
t1 <- Sys.time()
time.caret <- as.numeric(difftime(t1, t0, "secs"))

caret.predict <- predict(modell, newdata=data, type="raw")

rmse.caret<-rmse(data$Y, caret.predict)
rmse.caret

R2.caret <- cor(modell$finalModel$fit, data$Y[1:500])^2
R2.caret


t0 <- Sys.time()
set.seed(999)
data$A <- as.factor(data$A)
modele <- train(A ~ X1 + X2 + X3 + X4 + X5, data = data, method = "gbm", 
                preProcess = c('center', 'scale'), trControl = gbmCtrl, tuneGrid = gbmGrid)
t1 <- Sys.time()
time.caret <- as.numeric(difftime(t1, t0, "secs"))
caret.predict <- predict(modele, newdata=data, type="prob")[,2]

data$A <- as.numeric(as.character(data$A))
rmse.caret<-rmse(data$A, caret.predict)
rmse.caret



## gbm set.seed(999)
t0.gbm <- Sys.time()

model.gbm <- gbm(Y ~ X1 + X2 + X3 + X4 + X5, data = data[1:500,], distribution = "gaussian",
                     n.trees = 300, interaction.depth = 2,  # bag.fraction = 0.5, train.fraction = 0.8
                     n.minobsinnode = 5, shrinkage = 0.05, cv.folds = 10)
t1.gbm <- Sys.time()
time.gbm <- as.numeric(difftime(t1.gbm, t0.gbm, "secs"))

gbm.predict <- predict.gbm(model.gbm, newdata=data, 200)

rmse.gbm<-rmse(data$Y, gbm.predict)
rmse.gbm

R2.gbm <- cor(model.gbm$fit, data$Y[1:500])^2
R2.gbm

paste("diff.R2", R2.caret-R2.gbm)
paste("diff.RMSE", rmse.caret-rmse.gbm)

# gbm.fit test time
t0.fit <- Sys.time()
model.fit <- gbm.fit(x = data[1:500, !colnames(data) %in% c("Y", "A")], 
                     y = data[1:500,]$Y, distribution = "gaussian",
                     n.trees = 300, interaction.depth = 2,  # bag.fraction = 0.5,
                     n.minobsinnode = 5, shrinkage = 0.05,
                     nTrain = nrow(data))
t1.fit <- Sys.time()
time.fit <- as.numeric(difftime(t1.fit, t0.fit, "secs"))
# modell <- gbm.fit(x = data.xy[, colnames(data.xy) %in% form.r.terms], 
#                   y = data.xy$Y, distribution = "gaussian",
#                   n.trees = 100, interaction.depth = 2,  # bag.fraction = 0.5,
#                   n.minobsinnode = 5, shrinkage = 0.05, 
#                   nTrain = nrow(data.xy))



### 测评不同的gbm方法的速度与准确度
t0.out <- Sys.time()
y.hat <- est.outcome("Y ~ X1 + X2 + X3 + X4 + X5", data.node = data, method = "GBM", crossfit = TRUE, folds = 5)

t1.out <- Sys.time()
time.fit <- as.numeric(difftime(t1.out, t0.out, "secs"))
time.fit

mean(abs(data$Y - y.hat))
max(abs(data$Y - y.hat))
hist(data$Y - y.hat)



#### 1. gbmGrid <- expand.grid(n.trees = 300, interaction.depth= 2, n.minobsinnode = 5, shrinkage = 0.05)
# > time.fit
# [1] 4.322213
# > mean(abs(data$Y - y.hat))
# [1] 1.865802
# > max(abs(data$Y - y.hat))
# [1] 14.39048

#### 2. gbm(form.outcome, data = data.xy, distribution = "gaussian",
#           n.trees = 300, interaction.depth = 2,  # bag.fraction = 0.5, train.fraction = 0.8
#           n.minobsinnode = 5, shrinkage = 0.05, cv.folds = 10)
# > time.fit
# [1] 23.38977
# > mean(abs(data$Y - y.hat))
# [1] 1.622825
# > max(abs(data$Y - y.hat))
# [1] 9.186184
# > hist(data$Y - y.hat)


#### 3. gbm(form.outcome, data = data.xy, distribution = "gaussian",
#           n.trees = 300, interaction.depth = 2,  # bag.fraction = 0.5, train.fraction = 0.8
#           n.minobsinnode = 5, shrinkage = 0.05)

# > time.fit
# [1] 0.4530292
# > mean(abs(data$Y - y.hat))
# [1] 1.626184
# > max(abs(data$Y - y.hat))
# [1] 9.292727
# > hist(data$Y - y.hat)

#### 4. gbm.fit(x = aux.xy[, colnames(aux.xy) %in% form.r.terms],
#               y = aux.xy$Y, distribution = "gaussian",
#               n.trees = 300, interaction.depth = 2,  # bag.fraction = 0.5,
#               n.minobsinnode = 5, shrinkage = 0.05,
#               nTrain = nrow(aux.xy))
# > time.fit
# [1] 0.299
# > mean(abs(data$Y - y.hat))
# [1] 1.881433
# > max(abs(data$Y - y.hat))
# [1] 13.88836


#### 5.gbm.fit(x = aux.xy[, colnames(aux.xy) %in% form.r.terms],
#              y = aux.xy$Y, distribution = "gaussian",
#              n.trees = 300, interaction.depth = 1,  # bag.fraction = 0.5,
#              n.minobsinnode = 5, shrinkage = 0.1,
#              nTrain = nrow(aux.xy))
# > time.fit
# [1] 0.2840009
# > mean(abs(data$Y - y.hat))
# [1] 1.902335
# > max(abs(data$Y - y.hat))
# [1] 14.14553