load("simudata_multistudies.rda")


for(i in 1:4){
  simuname <- paste0("simu.",i)
  simuname.out <-  paste0("/home/bceuser/shend9/TranslatableForked/simulations/out/simu_multistudies.",i)
  obs.subset <- TRUE
  source("/home/bceuser/shend9/TranslatableForked/simulations/analysis_wrapper_multistudies.R")
  obs.subset <- FALSE
  simuname.out <- paste0("/home/bceuser/shend9/TranslatableForked/simulations/out/simu_multistudies.",i, "noa0")
  source("/home/bceuser/shend9/TranslatableForked/simulations/analysis_wrapper_multistudies.R")
}


# parallelize multiple simulations



run.multi.cases <- function(i){
  
  simuname <- paste0("simu.",i)
  simuname.out <-  paste0("/home/bceuser/shend9/TranslatableForked/simulations/out/simu_multistudies.",i)
  obs.subset <- TRUE
  source("/home/bceuser/shend9/TranslatableForked/simulations/analysis_wrapper_multistudies.R")
  obs.subset <- FALSE
  simuname.out <- paste0("/home/bceuser/shend9/TranslatableForked/simulations/out/simu_multistudies.",i, "noa0")
  source("analysis_wrapper_multistudies.R")
}

library(doParallel)  










no_cores <- detectCores() - 1
cl <- makeCluster(no_cores, type="FORK")  
registerDoParallel(cl)  

foreach(i=1:4) %dopar% {
  
  i <- i-1
  paste0("simu.",i)
}

stopCluster(cl)

