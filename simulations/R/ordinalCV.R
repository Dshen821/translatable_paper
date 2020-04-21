####################
# ordinal cv: creates lists of CV results: nsim specifies number of cvs repeated, simudata datasets indicates repeats
####################

ordinal.cv <- function(simudata = simudata,
                       n.sim = 10,
                       n.rep = 20,
                       obs.subset = TRUE,
                       response.type = "gaussian",
                       test.trial.name = "trial4",
                       topn = 5){
  
  # applies over all numbers of datasets
  cv.out <-  sapply(1:n.sim, function(s) {
    simu.out <- simudata[[s]] 
    x.names <- simu.out$x.names #all covariates
    causal.names <- simu.out$causal.names #causals: a0,b0,etc.
    if(!obs.subset){ # determines if a0 is included in the subset
      x.names <- setdiff(simu.out$x.names, "a.0")
      causal.names <- "a.1"
      simu.out$data <- simu.out$data %>% select(-a.0)
    }
    # subset training set
    tmpdata <- simu.out$data%>% filter(trial != test.trial.name)
    boot.cv(x=tmpdata, x.names=x.names,
            response.type = response.type, causal.names = causal.names, topn = topn, n.rep = n.rep, replace = FALSE)
  }, simplify=F)
  
  return(cv.out)
}
# 
# cv.out <- ordinal.cv(simudata = simudata, #simu.1,2,3,4
#                      n.sim = 10,
#                      n.rep = 20,
#                      obs.subset = TRUE,
#                      response.type = "gaussian",
#                      test.trial.name = "trial4",
#                      topn = 5)