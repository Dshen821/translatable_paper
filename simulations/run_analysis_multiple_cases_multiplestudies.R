library(pensim)
library(tidyverse)
library(knitr)
library(glmnet)
library(pheatmap)
library(doParallel)


source("/Users/u107586/translatable_paper/simulations/R/legacyfutureperformance.R")
source("/Users/u107586/translatable_paper/simulations/R/looCV.R")
source("/Users/u107586/translatable_paper/simulations/R/ordinalCV.R")

load("/Users/u107586/translatable_paper/simulations/simudata_multistudies.rda") #load simulated data: simu.1, ...., simu.4

# Set-up Cores, Backend
no_cores <- detectCores() - 1            
cl <- makeCluster(no_cores, type="FORK")
registerDoParallel(cl)

### Parallelized ###
foreach(i=1:4) %dopar% {
  set.seed(10)
  # Define simulation file
  simuname <- paste0("simu.",i)
  simuname.out <-  paste0("/Users/u107586/translatable_paper/simulations/out/simu_multistudies.",i)
  simudata <- get(simuname)
  obs.subset <- TRUE #causual in set- yes
  
  # 1) ORDINAL Cross Validation: For each of the 10 n.sim, create 20 CVs and calculate 4 errors
  cv.out <- ordinal.cv(simudata = simudata, obs.subset = obs.subset, response.type = "gaussian",
                       n.sim = 10, n.rep = 20, topn = 5, test.trial.name = "trial4")
  
  ## Calculate CV MSE: ROWS: CVs; COLs: Simulations
  
  # Test MSE: meansquared(err.2)
  lasso.mse.te.cv<- sapply(cv.out, function(j)sapply(j, function(i) mean((i$lasso.err.2) ^ 2))) 
  uni.mse.te.cv <- sapply(cv.out, function(j)sapply(j, function(i) mean((i$uni.err.2[,1]) ^ 2)))
  
  # Average MSE over CV folds: Average for each nsim
  lasso.mse.te.cv.mean <- apply(lasso.mse.te.cv,2,mean) # testing MSE in each data set (average across cross validations per sim
  uni.mse.te.cv.mean <- apply(uni.mse.te.cv, 2, mean)
  
  # Bias/optimism: Testing - Training MSE 
  lasso.diff.te.cv <- sapply(cv.out, function(j)sapply(j, function(i)mean(i$lasso.err.2^2) -mean(i$lasso.err.1^2)))  # in each data set in each CV testing
  uni.diff.te.cv <- sapply(cv.out, function(j)sapply(j, function(i)mean(i$uni.err.2[,1]^2) - mean(i$uni.err.1[,1]^2)))
  
  # Average Bias/Optimism over nsims:
  lasso.diff.te.cv.mean <- apply(lasso.diff.te.cv,2,mean) # testing MSE in each data set (average across sims)
  uni.diff.te.cv.mean <- apply(uni.diff.te.cv, 2, mean)
  
  ## 2) LOO Cross Validation: lasso
  loo.cv.out <- loo.cv(simudata = simudata, n.sim = 10, n.rep = 20, topn = 5,
                       obs.subset = obs.subset, response.type = "gaussian", test.trial.name = "trial4")
  # Test MSE: meansquared(err.2)
  lasso.mse.te.loo.cv <- sapply(loo.cv.out, function(j)sapply(j, function(i)mean((i$lasso.err.2)^2) )) #MSE in each data set in each CV testing
  uni.mse.te.loo.cv <- sapply(loo.cv.out, function(j)sapply(j, function(i)mean((i$uni.err.2[,1])^2)))
  # Average MSE over CV folds: Average for each nsim
  lasso.mse.te.loo.cv.mean <- apply(lasso.mse.te.loo.cv,2,mean) # testing MSE in each data set (average across CVs)
  uni.mse.te.loo.cv.mean <- apply(uni.mse.te.loo.cv, 2, mean)
  
  
  # Bias/optimism: Testing - Training MSE 
  lasso.diff.te.loo.cv <- sapply(loo.cv.out, function(j)sapply(j, function(i)mean(i$lasso.err.2^2) -mean(i$lasso.err.1^2)))  # in each data set in each CV testing
  uni.diff.te.loo.cv <- sapply(loo.cv.out, function(j)sapply(j, function(i)mean(i$uni.err.2[,1]^2) - mean(i$uni.err.1[,1]^2)))
  # Average Bias/Optimism over nsims:
  lasso.diff.te.loo.cv.mean <- apply(lasso.diff.te.loo.cv,2,mean) # testing MSE in each data set (average across CVs)
  uni.diff.te.loo.cv.mean <- apply(uni.diff.te.loo.cv, 2, mean)
  
  # 3) Build Models on Legacy + Model on Test
  res.v <- getLegacyTest(simudata)
  
  # Legacy/Future Error
  lasso.err.legacy <- sapply(res.v, function(i)i$lasso.err.1)
  lasso.err.future <- sapply(res.v, function(i)i$lasso.err.2)
  uni.1.err.legacy <- sapply(res.v, function(i)i$uni.err.1[,1])
  uni.1.err.future <- sapply(res.v, function(i)i$uni.err.2[,1])
  
  # MSE: mean(error ^2)
  lasso.mse.legacy <- apply(lasso.err.legacy^2,2,mean)
  lasso.mse.future <- apply(lasso.err.future^2,2,mean)
  uni.1.mse.legacy <- apply(uni.1.err.legacy^2,2,mean)
  uni.1.mse.future <- apply(uni.1.err.future^2,2,mean)
  
  # Adjusted MSE = legacy MSE + Bias
  lasso.legacy.adj.cv <- lasso.mse.legacy + lasso.diff.te.cv.mean # optimism from cv
  lasso.legacy.adj.loo.cv <- lasso.mse.legacy + lasso.diff.te.loo.cv.mean # optimism from loo
  
  uni.1.legacy.adj.cv <- uni.1.mse.legacy + uni.diff.te.cv.mean # optimism from cv
  uni.1.legacy.adj.loo.cv <- uni.1.mse.legacy + uni.diff.te.loo.cv.mean # optimism from loo
  
  # For Lasso,uni: ORDINAL CV testing mse,	LOO CV testing mse,  	ORDINAL CV adj mse,	LOO CV adj mse,	mse in future trial,	mse in legacy trial
  mse.mat <- cbind(lasso.mse.te.cv.mean,  lasso.mse.te.loo.cv.mean,  lasso.legacy.adj.cv, lasso.legacy.adj.loo.cv, lasso.mse.future,  lasso.mse.legacy, 
                   uni.mse.te.cv.mean, uni.mse.te.loo.cv.mean,  uni.1.legacy.adj.cv, uni.1.legacy.adj.loo.cv, uni.1.mse.future, uni.1.mse.legacy)
  
# 4) Save Files
  saveRDS(mse.mat, file=paste0(simuname.out, ".mse.mat.rds"))
  saveRDS(res.v, file=paste0(simuname.out, ".res.v.rds"))
  saveRDS(cv.out, file=paste0(simuname.out, ".cv.out.rds"))
  saveRDS(loo.cv.out, file=paste0(simuname.out, ".loo.cv.out.rds"))
  
  ### No Causal; Surragacy 
  set.seed(10)
  simuname.out <- paste0("/Users/u107586/translatable_paper/simulations/out/simu_multistudies.",i, "noa0")
  simudata <- get(simuname)
  obs.subset <- FALSE # surrogate in set

  # 1) ORDINAL Cross Validation: For each of the 10 n.sim, create 20 CVs and calculate 4 errors
  cv.out <- ordinal.cv(simudata = simudata, obs.subset = obs.subset, response.type = "gaussian",
                       n.sim = 10, n.rep = 20, topn = 5, test.trial.name = "trial4")
  
  ## Calculate CV MSE: ROWS: CVs; COLs: Simulations
  
  # Test MSE: meansquared(err.2)
  lasso.mse.te.cv<- sapply(cv.out, function(j)sapply(j, function(i) mean((i$lasso.err.2) ^ 2))) 
  uni.mse.te.cv <- sapply(cv.out, function(j)sapply(j, function(i) mean((i$uni.err.2[,1]) ^ 2)))
  
  # Average MSE over CV folds: Average for each nsim
  lasso.mse.te.cv.mean <- apply(lasso.mse.te.cv,2,mean) # testing MSE in each data set (average across cross validations per sim
  uni.mse.te.cv.mean <- apply(uni.mse.te.cv, 2, mean)
  
  # Bias/optimism: Testing - Training MSE 
  lasso.diff.te.cv <- sapply(cv.out, function(j)sapply(j, function(i)mean(i$lasso.err.2^2) -mean(i$lasso.err.1^2)))  # in each data set in each CV testing
  uni.diff.te.cv <- sapply(cv.out, function(j)sapply(j, function(i)mean(i$uni.err.2[,1]^2) - mean(i$uni.err.1[,1]^2)))
  
  # Average Bias/Optimism over nsims:
  lasso.diff.te.cv.mean <- apply(lasso.diff.te.cv,2,mean) # testing MSE in each data set (average across sims)
  uni.diff.te.cv.mean <- apply(uni.diff.te.cv, 2, mean)
  
  ## 2) LOO Cross Validation: lasso
  loo.cv.out <- loo.cv(simudata = simudata, n.sim = 10, n.rep = 20, topn = 5,
                       obs.subset = obs.subset, response.type = "gaussian", test.trial.name = "trial4")
  # Test MSE: meansquared(err.2)
  lasso.mse.te.loo.cv <- sapply(loo.cv.out, function(j)sapply(j, function(i)mean((i$lasso.err.2)^2) )) #MSE in each data set in each CV testing
  uni.mse.te.loo.cv <- sapply(loo.cv.out, function(j)sapply(j, function(i)mean((i$uni.err.2[,1])^2)))
  # Average MSE over CV folds: Average for each nsim
  lasso.mse.te.loo.cv.mean <- apply(lasso.mse.te.loo.cv,2,mean) # testing MSE in each data set (average across CVs)
  uni.mse.te.loo.cv.mean <- apply(uni.mse.te.loo.cv, 2, mean)
  
  
  # Bias/optimism: Testing - Training MSE 
  lasso.diff.te.loo.cv <- sapply(loo.cv.out, function(j)sapply(j, function(i)mean(i$lasso.err.2^2) -mean(i$lasso.err.1^2)))  # in each data set in each CV testing
  uni.diff.te.loo.cv <- sapply(loo.cv.out, function(j)sapply(j, function(i)mean(i$uni.err.2[,1]^2) - mean(i$uni.err.1[,1]^2)))
  # Average Bias/Optimism over nsims:
  lasso.diff.te.loo.cv.mean <- apply(lasso.diff.te.loo.cv,2,mean) # testing MSE in each data set (average across CVs)
  uni.diff.te.loo.cv.mean <- apply(uni.diff.te.loo.cv, 2, mean)
  
  # 3) Build Models on Legacy + Model on Test
  res.v <- getLegacyTest(simudata)
  
  # Legacy/Future Error
  lasso.err.legacy <- sapply(res.v, function(i)i$lasso.err.1)
  lasso.err.future <- sapply(res.v, function(i)i$lasso.err.2)
  uni.1.err.legacy <- sapply(res.v, function(i)i$uni.err.1[,1])
  uni.1.err.future <- sapply(res.v, function(i)i$uni.err.2[,1])
  
  # MSE: mean(error ^2)
  lasso.mse.legacy <- apply(lasso.err.legacy^2,2,mean)
  lasso.mse.future <- apply(lasso.err.future^2,2,mean)
  uni.1.mse.legacy <- apply(uni.1.err.legacy^2,2,mean)
  uni.1.mse.future <- apply(uni.1.err.future^2,2,mean)
  
  # Adjusted MSE = legacy MSE + Bias
  lasso.legacy.adj.cv <- lasso.mse.legacy + lasso.diff.te.cv.mean # optimism from cv
  lasso.legacy.adj.loo.cv <- lasso.mse.legacy + lasso.diff.te.loo.cv.mean # optimism from loo
  
  uni.1.legacy.adj.cv <- uni.1.mse.legacy + uni.diff.te.cv.mean # optimism from cv
  uni.1.legacy.adj.loo.cv <- uni.1.mse.legacy + uni.diff.te.loo.cv.mean # optimism from loo
  
  # For Lasso,uni: ORDINAL CV testing mse,	LOO CV testing mse,  	ORDINAL CV adj mse,	LOO CV adj mse,	mse in future trial,	mse in legacy trial
  mse.mat <- cbind(lasso.mse.te.cv.mean,  lasso.mse.te.loo.cv.mean,  lasso.legacy.adj.cv, lasso.legacy.adj.loo.cv, lasso.mse.future,  lasso.mse.legacy, 
                   uni.mse.te.cv.mean, uni.mse.te.loo.cv.mean,  uni.1.legacy.adj.cv, uni.1.legacy.adj.loo.cv, uni.1.mse.future, uni.1.mse.legacy)
  
  # 4) Save Files
  saveRDS(mse.mat, file=paste0(simuname.out, ".mse.mat.rds"))
  saveRDS(res.v, file=paste0(simuname.out, ".res.v.rds"))
  saveRDS(cv.out, file=paste0(simuname.out, ".cv.out.rds"))
  saveRDS(loo.cv.out, file=paste0(simuname.out, ".loo.cv.out.rds"))
}

stopCluster(cl)







# Sequential
# for (i in 1:4){
# simuname <- paste0("simu.",i) # simu.1, 2, 3, 4, ...
# simuname.out <-  paste0("/home/bceuser/shend9/TranslatableForked/simulations/out/simu_multistudies.",i) #simu_multistudies.i xxxxxx.rds
# obs.subset <- TRUE
# source("/home/bceuser/shend9/TranslatableForked/simulations/analysis_wrapper_multistudies.R", local = TRUE)
# obs.subset <- FALSE
# simuname.out <- paste0("/home/bceuser/shend9/TranslatableForked/simulations/out/simu_multistudies.",i, "noa0")
# source("/home/bceuser/shend9/TranslatableForked/simulations/analysis_wrapper_multistudies.R", local = TRUE)
# }


# Parallel
# Set cores, register backend
# no_cores <- detectCores() - 1
# cl <- makeCluster(no_cores, type="FORK")
# registerDoParallel(cl)
# 
# foreach(i=1:4) %dopar% {
#   simuname <- paste0("simu.",i)
#   simuname.out <-  paste0("/home/bceuser/shend9/TranslatableForked/simulations/out/simu_multistudies.",i)
#   obs.subset <- TRUE
#   source("/home/bceuser/shend9/TranslatableForked/simulations/analysis_wrapper_multistudies.R", local = TRUE)
#   obs.subset <- FALSE
#   simuname.out <- paste0("/home/bceuser/shend9/TranslatableForked/simulations/out/simu_multistudies.",i, "noa0")
#   source("/home/bceuser/shend9/TranslatableForked/simulations/analysis_wrapper_multistudies.R", local = TRUE)
# }
# 
# stopCluster(cl)

