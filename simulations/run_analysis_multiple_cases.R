load("simudata.rda")
for(i in 1:4){
  simuname <- paste0("simu.",i)
  source("analysis_wrapper.R")
}