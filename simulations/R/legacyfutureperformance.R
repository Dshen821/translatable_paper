
#####################################
# Build model using legacy and test in future trial
#####################################
# loop through simulated data sets

getLegacyTest <- function(simudata = simudata,
                          n.sim = 10,
                          n.rep = 20,
                          obs.subset = TRUE,
                          response.type = "gaussian",
                          test.trial.name = "trial4",
                          topn = 5){
  res.v <- sapply(1:n.sim, function(s) {
    simu.out <- simudata[[s]]
    x.names <- simu.out$x.names
    causal.names <- simu.out$causal.names
    if(!obs.subset){ # what if we don't observe the real causal signal (but only its surrogancies)
      x.names <- setdiff(simu.out$x.names, "a.0")
      causal.names <- "a.1"
      simu.out$data <- simu.out$data %>% select(-a.0)
    }
    # Run models on legacy and trial data: returns top names, estimates, true estimates, predicted
    run.multi(data.legacy=simu.out$data %>% filter(trial!=test.trial.name), 
              data.future=simu.out$data %>% filter(trial==test.trial.name), x.names=x.names,
              response.type = response.type, causal.names = causal.names)}, 
    simplify = F)
  
  return(res.v)
  
}
# 
# res.v <- getLegacyTest(simudata)