load("simudata.rda")
for(i in 1:4){
  
  simuname.out <- simuname <- paste0("simu.",i)
  obs.subset <- TRUE
  source("analysis_wrapper.R")
  obs.subset <- FALSE
  simuname.out <- paste0("simu.",i, "noa0")
  source("analysis_wrapper.R")
}