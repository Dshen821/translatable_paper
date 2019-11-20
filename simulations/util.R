
simu <- function(nvars, cors, associations, nsamples, response, logisticintercept, 
                 n.trial1, n.trial2, response.use, eps, causal.v) {
  data.sim <- create.data(nvars=nvars,
                          cors=cors,
                          associations=associations,
                          firstonly=rep(TRUE, length(nvars)),
                          nsamples=nsamples,
                          response=response,
                          logisticintercept=logisticintercept)$data 
  
  if(response.use=="gaussian") data.sim$outcome <- colSums(t(data.sim[causal.v]) * associations) + rnorm(nsamples, 0, eps)
  if(response.use=="binomial") data.sim <- data.sim %>% dplyr::mutate(outcome=as.numeric(outcome)-1)
  
  cor.mat <- cor(data.sim[,c(causal.v, "outcome")])
  data.trial1 <- data.sim[1:n.trial1,]
  data.trial2 <- data.sim[(n.trial1+1): nsamples, ]
  rownames(data.trial1) <- paste0("trial1.", 1:n.trial1)
  rownames(data.trial2) <- paste0("trial2.", 1:n.trial2)
  x.names <- setdiff(names(data.sim), "outcome")
  
  out <- list(data.trial1=data.trial1, data.trial2=data.trial2, 
              x.names=x.names, cor.mat=cor.mat)
}



run.lasso <- function(x,y, alpha=1,top=NULL, family="binomial", ...){
  nona <- which(!is.na(rowMeans(x)))
  lasso.cv<- cv.glmnet(x=x[nona,], y=y[nona],alpha=alpha,family=family,...)
  lasso.pen <- lasso.cv$lambda.min #optimal lambda
  #lasso.pen #minimal shrinkage
  lasso.fit <-glmnet(x = x[nona,], y=y[nona],alpha=alpha,family=family,lambda = lasso.pen,...) #estimate the
  lasso.coef <- as.vector(matrix(coef(lasso.fit)))
  names(lasso.coef) <- coef(lasso.fit)@Dimnames[[1]]
  lasso.coef <- lasso.coef[-1]
  lasso.coef.od <- order(abs(lasso.coef),decreasing=T)
  lasso.coef.sort <- lasso.coef[lasso.coef.od]
  if(!is.null(top))nznames <- lasso.coef.sort[which(lasso.coef.sort!=0)][1:min(top, sum(lasso.coef.sort!=0))] else nznames <- lasso.coef.sort[which(lasso.coef.sort!=0)]
  out <- list(fit=lasso.fit, nonzero.names = nznames)
}


run.gaussian.uni <- function(x, y){
  vals <- sapply(colnames(x), function(i)coef(summary(lm(y~x[,i])))[2,c("Estimate","Pr(>|t|)")])
  vals.sort <- vals[,order(vals["Pr(>|t|)",])]
  vals.sort
}

run.glm.uni <- function(x, y, family="binomial"){
  vals <- sapply(colnames(x), function(i)coef(summary(glm(y~x[,i]), family=family))[2,c("Estimate","Pr(>|t|)")])
  vals.sort <- vals[,order(vals["Pr(>|t|)",])]
  vals.sort
}

run.multi <- function(data.trial1, data.trial2, x.names, response.use, causal.v, topn=5){
  
  lasso.out <- run.lasso(x = data.matrix(data.trial1[x.names]), y = data.trial1$outcome, family=response.use)
  lasso.top <- lasso.out$nonzero.names
  lasso.top.names <- names(lasso.top)
  if(length(lasso.top.names)==0) lasso.top.names <- lasso.top <-  NA
  uni.top <- run.glm.uni(x = data.matrix(data.trial1[x.names]), y = data.trial1$outcome, family=response.use)
  uni.top.names <-  colnames(uni.top)
  uni.trial2 <- run.glm.uni(x = data.matrix(data.trial2[uni.top.names]), y = data.trial2$outcome, family=response.use)
  out <- list(lasso.top.names = lasso.top.names, uni.top.names = uni.top.names[topn], 
              top.uni.est.1 = uni.top[1,1], top.uni.est.2 = uni.trial2[1, uni.top.names[1]],
              top.lasso.est.1 = uni.top[1,lasso.top.names[1]], top.lasso.est.2 = uni.trial2[1, lasso.top.names[1]],
              true.est.1 = uni.top[1, causal.v],true.est.2=uni.trial2[1, causal.v])
  
  }