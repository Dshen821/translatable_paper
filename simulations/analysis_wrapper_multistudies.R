# Analysis wrapper


library(pensim)
library(tidyverse)
library(knitr)
library(glmnet)
library(pheatmap)
source("util.R")
set.seed(10)

#simuname <- "simu1"

simudata <- get(simuname)

n.sim <- 10 # how many simulations to analyze - set to 10 for now to save time
n.rep <- 20 # how many CVs to run
response.type <- "gaussian"
test.trial.name <- "trial4"


####################
# ordinal cv
####################
cv.out <-  sapply(1:n.sim, function(s) {
  simu.out <- simudata[[s]] 
  x.names <- simu.out$x.names
  causal.names <- simu.out$causal.names
  if(!obs.subset){
    x.names <- setdiff(simu.out$x.names, "a.0")
    causal.names <- "a.1"
    simu.out$data <- simu.out$data %>% select(-a.0)
  }
  tmpdata <- simu.out$data%>% filter(trial != test.trial.name)
  boot.cv(x=tmpdata, x.names=x.names,
        response.type = response.type, causal.names = causal.names, topn = 5, n.rep = n.rep, replace = FALSE)
}, simplify=F)



# testing mse per cv
lasso.mse.te.cv<- sapply(cv.out, function(j)sapply(j, function(i)mean((i$lasso.err.2)^2) )) #MSE in each data set in each CV testing
uni.mse.te.cv <- sapply(cv.out, function(j)sapply(j, function(i)mean((i$uni.err.2[,1])^2)))
# rows: cv, cols: simulations
# average mse difference
lasso.mse.te.cv.mean <- apply(lasso.mse.te.cv,2,mean) # testing MSE in each data set (average across CVs)
uni.mse.te.cv.mean <- apply(uni.mse.te.cv, 2, mean)



# testing mse - training mse, per cv
lasso.diff.te.cv<- sapply(cv.out, function(j)sapply(j, function(i)mean(i$lasso.err.2^2) -mean(i$lasso.err.1^2)))  # in each data set in each CV testing
uni.diff.te.cv <- sapply(cv.out, function(j)sapply(j, function(i)mean(i$uni.err.2[,1]^2) - mean(i$uni.err.2[,1]^2)))
# rows: cv, cols: simulations
# average mse difference
lasso.diff.te.cv.mean <- apply(lasso.diff.te.cv,2,mean) # testing MSE in each data set (average across CVs)
uni.diff.te.cv.mean <- apply(uni.diff.te.cv, 2, mean)


####################
# leave one study out cv
####################
loo.cv.out <-  sapply(1:n.sim, function(s) {
  simu.out <- simudata[[s]] 
  x.names <- simu.out$x.names
  causal.names <- simu.out$causal.names
  if(!obs.subset){
    x.names <- setdiff(simu.out$x.names, "a.0")
    causal.names <- "a.1"
    simu.out$data <- simu.out$data %>% select(-a.0)
  }
  tmpdata <- simu.out$data%>% filter(trial != test.trial.name)
  boot.cv(x=tmpdata, x.names=x.names,
        response.type = response.type, causal.names = causal.names, 
        topn = 5, n.rep = n.rep, replace = FALSE,
        name.mat = sapply(unique(tmpdata$trial), 
                          function(i)which(tmpdata$trial!=i) , 
                          simplify = F))
}, simplify=F)

# testing mse per cv
lasso.mse.te.loo.cv<- sapply(loo.cv.out, function(j)sapply(j, function(i)mean((i$lasso.err.2)^2) )) #MSE in each data set in each CV testing
uni.mse.te.loo.cv <- sapply(loo.cv.out, function(j)sapply(j, function(i)mean((i$uni.err.2[,1])^2)))
# rows: cv, cols: simulations
# average mse difference
lasso.mse.te.loo.cv.mean <- apply(lasso.mse.te.loo.cv,2,mean) # testing MSE in each data set (average across CVs)
uni.mse.te.loo.cv.mean <- apply(uni.mse.te.loo.cv, 2, mean)


# testing mse - training mse, per cv
lasso.diff.te.loo.cv<- sapply(loo.cv.out, function(j)sapply(j, function(i)mean(i$lasso.err.2^2) -mean(i$lasso.err.1^2)))  # in each data set in each CV testing
uni.diff.te.loo.cv <- sapply(loo.cv.out, function(j)sapply(j, function(i)mean(i$uni.err.2[,1]^2) - mean(i$uni.err.2[,1]^2)))
# rows: cv, cols: simulations
# average mse difference
lasso.diff.te.loo.cv.mean <- apply(lasso.diff.te.loo.cv,2,mean) # testing MSE in each data set (average across CVs)
uni.diff.te.loo.cv.mean <- apply(uni.diff.te.loo.cv, 2, mean)



#####################################
# Build model using legacy and test in future trial
#####################################
# loop through simulated data sets
res.v <- sapply(1:n.sim, function(s) {
  simu.out <- simudata[[s]]
  x.names <- simu.out$x.names
  causal.names <- simu.out$causal.names
  if(!obs.subset){ # what if we don't observe the real causal signal (but only its surrogancies)
    x.names <- setdiff(simu.out$x.names, "a.0")
    causal.names <- "a.1"
    simu.out$data <- simu.out$data %>% select(-a.0)
  }
  run.multi(data.legacy=simu.out$data %>% filter(trial!=test.trial.name), 
            data.future=simu.out$data %>% filter(trial==test.trial.name), x.names=x.names,
            response.type = response.type, causal.names = causal.names)}, 
  simplify = F)



# future: true performance in future trial
# legancy: performance in legancy trial. expect to be overfitted

lasso.err.legacy<- sapply(res.v, function(i)i$lasso.err.1)
lasso.err.future<- sapply(res.v, function(i)i$lasso.err.2)
uni.1.err.legacy<- sapply(res.v, function(i)i$uni.err.1[,1])
uni.1.err.future<- sapply(res.v, function(i)i$uni.err.2[,1])

lasso.mse.legacy <- apply(lasso.err.legacy^2,2,mean)
lasso.mse.future <- apply(lasso.err.future^2,2,mean)
uni.1.mse.legacy <- apply(uni.1.err.legacy^2,2,mean)
uni.1.mse.future <- apply(uni.1.err.future^2,2,mean)

lasso.legacy.adj.cv <- lasso.mse.legacy + lasso.diff.te.cv.mean
lasso.legacy.adj.loo.cv <- lasso.mse.legacy + lasso.diff.te.loo.cv.mean

uni.1.legacy.adj.cv <- uni.1.mse.legacy + uni.diff.te.cv.mean
uni.1.legacy.adj.loo.cv <- uni.1.mse.legacy + uni.diff.te.loo.cv.mean


mse.mat <- cbind(lasso.mse.te.cv.mean, lasso.mse.te.loo.cv.mean, lasso.legacy.adj.cv, lasso.legacy.adj.loo.cv,
                 lasso.mse.future, lasso.mse.legacy, 
                 uni.mse.te.cv.mean, uni.mse.te.loo.cv.mean, uni.1.legacy.adj.cv, uni.1.legacy.adj.loo.cv
                 ,uni.1.mse.future, uni.1.mse.legacy)


saveRDS(mse.mat, file=paste0(simuname.out, ".mse.mat.rds"))
saveRDS(res.v, file=paste0(simuname.out, ".res.v.rds"))
saveRDS(cv.out, file=paste0(simuname.out, ".cv.out.rds"))
saveRDS(loo.cv.out, file=paste0(simuname.out, ".loo.cv.out.rds"))