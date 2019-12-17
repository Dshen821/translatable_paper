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

n.sim <- 100
n.rep <- 50
response.type <- "gaussian"
res.v <- sapply(1:n.sim, function(s) {
  simu.out <- simu.1[[s]]
  run.multi(data.trial1=simu.out$data %>% filter(trial=="trial1"), data.trial2=simu.out$data %>% filter(trial=="trial2"), x.names=simu.out$x.names,
            response.type = response.type, causal.names = simu.out$causal.names)}, simplify = F)


true.diff <- sapply(res.v, function(i)i$true.est.1-i$true.est.2)
lasso.diff <- sapply(res.v, function(i)i$top.lasso.est.1-i$top.lasso.est.2)
lasso.diff <- unlist(lasso.diff[sapply(lasso.diff, function(i)!is.na(i[1]))])
sig.lasso.diff <- sapply(res.v, function(i)i$sig.lasso.est.1-i$sig.lasso.est.2)
sig.lasso.diff <- unlist(sig.lasso.diff[sapply(sig.lasso.diff, function(i)!is.na(i[1]))])
uni.diff <- sapply(res.v, function(i)i$top.uni.est.1-i$top.uni.est.2)

plot.df <- data.frame(diff=c(true.diff, uni.diff, lasso.diff, sig.lasso.diff), method=c(rep("true causal biomarker", length(true.diff)),
                                                                                        rep("uni variate picked", length(uni.diff)),
                                                                                        rep("lasso picked", length(lasso.diff)),
                                                                                        rep("lasso built signature", length(sig.lasso.diff))))
plot.df$method <- factor(plot.df$method, levels=unique(plot.df$method))




lasso.err.trial1<- sapply(res.v, function(i)i$lasso.err.1)
lasso.err.trial2<- sapply(res.v, function(i)i$lasso.err.2)
uni.1.err.trial1<- sapply(res.v, function(i)i$uni.err.1[,1])
uni.1.err.trial2<- sapply(res.v, function(i)i$uni.err.2[,1])

lasso.mse.trial1 <- apply(lasso.err.trial1^2,2,mean)
lasso.mse.trial2 <- apply(lasso.err.trial2^2,2,mean)
uni.1.mse.trial1 <- apply(uni.1.err.trial1^2,2,mean)
uni.1.mse.trial2 <- apply(uni.1.err.trial2^2,2,mean)

mse.plot.df <- data.frame(mse=c(lasso.mse.trial1, lasso.mse.trial2, uni.1.mse.trial1, uni.1.mse.trial2), 
                          method=rep(c("lasso.trial1", "lasso.trial2", "unitop.trial1","unitop.trial2"), each=length(lasso.mse.trial1)))
mse.plot.df$method <- factor(mse.plot.df$method, levels=unique(mse.plot.df$method))


boot.out <-  sapply(1:n.sim, function(s) {
  simu.out <- simu.1[[s]]
  boot.cv(x=simu.out$data, x.names=simu.out$x.names,
          response.type = response.type, causal.names = simu.out$causal.names, topn = 5, n.rep = n.rep, replace = TRUE)
}, simplify=F)

# mse diff per bootstrap
lasso.adj.boot<- sapply(boot.out, function(j)sapply(j, function(i)mean((i$lasso.err.1)^2) - mean((i$lasso.err.2)^2)))
uni.adj.boot <- sapply(boot.out, function(j)sapply(j, function(i)mean((i$uni.err.1[,1])^2) - mean((i$uni.err.2[,1])^2)))
# rows: bootstraps, cols: simulations
# average mse difference
lasso.adj.boot.mean <- apply(lasso.adj.boot,2,mean)
uni.adj.boot.mean <- apply(uni.adj.boot, 2, mean)

lasso.mse.trial2.boot <- t(t(lasso.mse.trial2) + lasso.adj.boot.mean)
uni.1.mse.trial2.boot <- t(t(uni.1.mse.trial2) + uni.adj.boot.mean)

mse.plot.df <- rbind(mse.plot.df, data.frame(mse=c(lasso.mse.trial2.boot, uni.1.mse.trial2.boot),
                                             method=rep(c( "lasso.trial2.bootcorrection", "unitop.trial2.bootcorrection"), each=length(lasso.mse.trial1))))

#ggplot(mse.plot.df, aes(x=mse, color=method)) + geom_density() + 


saveRDS(plot.df, file = paste0(simuname,".plot.df.rds"))
saveRDS(mse.plot.df, file = paste0(simuname,".mse.plot.df.rds"))
saveRDS(res.v, file=paste0(simuname, ".res.v.rds"))
saveRDS(res.v, file=paste0(simuname, ".boot.out.rds"))