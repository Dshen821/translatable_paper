library("foreach")
library("doParallel")

cores <- detectCores()
cl <- makeCluster(7)
registerDoParallel(cl)

#dummy example

foreach(i=1:3) %dopar% sqrt(i)

registerDoParallel(cores=7)
foreach(i=1:3) %dopar% sqrt(i)

# To use multicore-like functionality, we would specify the number of cores to use instead (but
# note that on Windows, attempting to use more than one core with parallel results in an error):


# Bootstrapping Ex

cl <- makeCluster(7)
registerDoParallel(cl)

#parallel %dopar%
x <- iris[which(iris[,5] != "setosa"), c(1,5)]
trials <- 10000
ptime <- system.time({
  r <- foreach(icount(trials), .combine=cbind) %dopar% {
    ind <- sample(100, 100, replace=TRUE)
    result1 <- glm(x[ind,2]~x[ind,1], family=binomial(logit))
    coefficients(result1)
    }
  })
ptime #13.74


# sequential %do%
stime <- system.time({
   r <- foreach(icount(trials), .combine=cbind) %do% {
     ind <- sample(100, 100, replace=TRUE)
     result1 <- glm(x[ind,2]~x[ind,1], family=binomial(logit))
     coefficients(result1)
     }
   })
stime # 25.414


#Stopping your cluster
stopCluster(cl)


ind <- sample(100, 100, replace=TRUE)
result1 <- glm(x[ind,2]~x[ind,1], family=binomial(logit))
coefficients(result1)

cbind(coef(result1), coef(result1))