load("simudata_multistudies.rda")
for(i in 1:4){
  simuname <- paste0("simu.",i)
  simuname.out <-  paste0("simu_multistudies.",i)
  obs.subset <- TRUE
  source("analysis_wrapper_multistudies.R")
  obs.subset <- FALSE
  simuname.out <- paste0("simu_multistudies.",i, "noa0")
  source("analysis_wrapper_multistudies.R")
}