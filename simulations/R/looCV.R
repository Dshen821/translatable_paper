
####################
# leave one study out cv
####################

loo.cv <- function(simudata = simudata,
                   n.sim = 10,
                   n.rep = 20,
                   obs.subset = TRUE,
                   response.type = "gaussian",
                   test.trial.name = "trial4",
                   topn = 5){
  
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
            name.mat = sapply(unique(tmpdata$trial),              ########## extra in LOOCV
                              function(i)which(tmpdata$trial!=i) ,  ########
                              simplify = F))
  }, simplify=F)
  
  return(loo.cv.out)
}

# loo.cv.out <- loo.cv(simudata = simudata,
#                      n.sim = 10,
#                      n.rep = 20,
#                      obs.subset = TRUE,
#                      response.type = "gaussian",
#                      test.trial.name = "trial4",
#                      topn = 5)
