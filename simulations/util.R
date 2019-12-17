
#' @param n.causal number of independent causal features to generate. they will be named as. c("a.1","b.1") (one feature within a group)
#' @param n.cor.causal number of correlated features to generate, for each group of features. e.g. c(100, 200).
#' feature names will be a.1, a.2....a.100, b.1,...b.200
#' @param cor.causal correlation structure. e.g.  c(.9, .1). could be a matrix of column 2 (if correlation in trial 1 and 2 are different)
#' @param coef.causal coefficient between causal features and outcome. e.g. c(1, .2). 
#' If it is a vector, assuming trial 1 and 2 share the same coefficients.
#' If trial 1 and 2 share different coefficients, it should be a matrix with 2 columns, where column 1 (2)
#' indicates trial 1(2)
#' @param n.trial1,n.trial2 number of samples in trial 1 and 2
#' @param n.noise response type, gaussian or binomial
#' @param outcome.sd when gaussian outcome is generated, noise follows N(0, eps)
#' @param response not used
#' @param logisticintercept not used
#' TODO: replace create.data with alpha
#' TODO: X shift between trial1 and trial2
simu <- function(n.trial1, n.trial2, 
                 n.causal = 2, n.cor.causal = c(10,20), cor.causal = c(.9, .8), coef.causal=c(1,0),
                 n.noise=50, outcome.sd=.1,
                 shift.mean = c(0, 0), shift.sd=c(.1, .1)) {
  
  nsamples <- n.trial1 + n.trial2
  coef.causal <- data.matrix(coef.causal)
  cor.causal <- data.matrix(cor.causal)
  if(n.causal>0)stopifnot(n.causal == length(n.cor.causal)) else stopifnot(is.null(n.cor.causal))
  stopifnot(n.causal+n.noise>0)
  
  
  # simulate causal features and noise feactures
  names.no.cor <- causal.names <- NULL
  if(n.causal>0) {
    causal.names <- paste0(letters[1:n.causal],".0")
    names.no.cor <- c(names.no.cor, causal.names) # causal are a.0, b.0...
  }
  if(n.noise>0) names.no.cor <- c(names.no.cor, paste0("noise.",1:n.noise-1)) # noises are noise.0, noise.1, ....
  x.no.cor <- sapply(1: (n.causal+n.noise), function(i)rnorm(nsamples))
  colnames(x.no.cor) <- names.no.cor
  x.no.cor <- as.data.frame(x.no.cor)
  x.no.cor <- x.no.cor %>% mutate(trial=c(rep("trial1", n.trial1), rep("trial2", n.trial2)))
  
  
  # population shift
  if(unique(shift.mean)!=0){
    for (fe in 1:n.causal){
      if(shift.mean[fe]>0){
        x.var <- paste0(letters[fe],".0")
        x.no.cor[which(x.no.cor$trial=='trail.2'), x.var] <- x.no.cor[which(x.no.cor$trial=='trail.2'), x.var] + rnorm(n.trial2 ,shift.mean[fe], shift.sd[fe])
      }
      }
  }
  
  # Simulate correlated features: allow different correlation matrix in trial 1 and 2
  x <- x.no.cor
  if(n.causal>0){
    for (fe in 1:n.causal){
      if(n.cor.causal[fe]>0){
        prefix <- letters[fe]
        cor.causal.use <- cor.causal[fe, ] # cor to use for each causal feature
        cor.v <- c(rep(cor.causal.use, n.trial1), rep(cor.causal.use[length(cor.causal.use)], n.trial2)) 
        # cor for each patient (all pts have same cor structure, or each trial has diff cor structure)
        x.cor.tmp <- sapply(1:n.cor.causal[fe], function(i) x.no.cor[, paste0(prefix,".0")] * sqrt(cor.v) + rnorm(nsamples) * sqrt(1-cor.v))
        colnames(x.cor.tmp) <- paste0(prefix, ".", 1:ncol(x.cor.tmp))
        x <- cbind(x, data.frame(x.cor.tmp))
      }
    }
    
  }
  
  
  
  # Simulate outcome
  if(n.causal==0) x$outcome <- rnorm(nsamples)
  
  if(n.causal>0){ 
    x$outcome <- NULL
    for(i in 1:2){
    coef.use <- coef.causal[,ifelse(ncol(coef.causal)==1, 1, i)]
    tmp <-which(x$trial==paste0("trial",i))
    x$outcome[tmp] <- colSums(t(x[tmp, causal.names]) * coef.use) + 
      rnorm(length(tmp), 0, outcome.sd)
  }
  }
  
  x <- x[order(names(x))]
  rownames(x) <- paste0('s.', 1:nrow(x))
  x.names <- setdiff(names(x), c("outcome","trial"))
  
  out <- list(data = x, x.names=x.names, causal.names=causal.names)
}



heat.fun <- function(simu.out){
  data.tmp <- simu.out$data 
  tmp <- data.tmp %>% select(outcome, trial)
  pheatmap(t(data.tmp[,simu.out$x.names]),
           cluster_cols=F, cluster_rows = F,
           annotation_row = NA,
           annotation_col=tmp)
  
}


run.lasso <- function(x,y, alpha=1,top=NULL, family="binomial", ...){
  nona <- which(!is.na(rowMeans(x)))
  lasso.cv<- cv.glmnet(x=x[nona,], y=y[nona],alpha=alpha,family=family,...)
  lasso.pen <- lasso.cv$lambda.min #optimal lambda
  #lasso.pen #minimal shrinkage
  lasso.fit <-glmnet(x = x[nona,], y=y[nona],alpha=alpha,family=family,lambda = lasso.pen,...) #estimate the
  lasso.coef <- as.vector(data.matrix(coef(lasso.fit)))
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

run.glm.uni <- function(x, y, family="binomial", save.model=FALSE){
  if(!save.model){
    vals <- sapply(colnames(x), function(i)coef(summary(glm(y~x[,i], family=family)))[2,c("Estimate","Pr(>|t|)")])
    vals.sort <- vals[,order(vals["Pr(>|t|)",])]
    out <- vals.sort
  }
  if(save.model){
    mod <- sapply(colnames(x), function(i){
      tmp <- data.frame(y0 = y, x0= x[,i]) 
      glm(y0 ~ x0, data=tmp, family=family)}
      , simplify = F)
    vals <- data.matrix(sapply(mod, function(i)coef(summary(i))[2,c("Estimate","Pr(>|t|)")]))
    vals.sort <- vals[,order(vals["Pr(>|t|)",])]
    pred <- data.matrix(sapply(mod, function(i)predict(i)))
    out <- list(model=mod, vals.sort=vals.sort, pred=pred)
  }
  out
}

run.multi <- function(data.trial1, data.trial2, x.names, response.type, causal.names, topn=5){
  
  lasso.out <- run.lasso(x = data.matrix(data.trial1[x.names]), y = data.trial1$outcome, family=response.type)
  lasso.top <- lasso.out$nonzero.names
  lasso.top.names <- names(lasso.top)
  if(length(lasso.top.names)==0) lasso.top.names <- lasso.top <-  NA
  uni.top <- run.glm.uni(x = data.matrix(data.trial1[x.names]), y = data.trial1$outcome, family=response.type, save.model = TRUE)
  uni.top.val <- uni.top$vals.sort
  uni.top.names <-  colnames(uni.top$vals.sort)
  uni.trial2 <- run.glm.uni(x = data.matrix(data.trial2[uni.top.names]), y = data.trial2$outcome, family=response.type)
  uni.pred.1 <- uni.top$pred[,uni.top.names]
  uni.pred.2 <- sapply(1:length(uni.top.names), function(i){
    tmp <- data.frame(y0 = data.trial2$outcome, x0= data.trial2[[uni.top.names[i]]]) 
    predict(uni.top$mod[[uni.top.names[i]]], newdata = tmp)
    })
  colnames(uni.pred.2) <- uni.top.names

  lasso.pred.trial1 <- predict(newx = data.matrix(data.trial1[x.names]), object=lasso.out$fit)
  sig.lasso.est.1 <- run.glm.uni(x = data.matrix(lasso.pred.trial1), y = data.trial1$outcome, family=response.type)["Estimate"]
  lasso.pred.trial2 <- predict(newx = data.matrix(data.trial2[x.names]), object=lasso.out$fit)
  sig.lasso.est.2 <- run.glm.uni(x = data.matrix(lasso.pred.trial2), y = data.trial2$outcome, family=response.type)["Estimate"]
  
  
  
  out <- list(lasso.top.names = lasso.top.names, uni.top.names = uni.top.names[1:topn], 
              top.uni.est.1 = uni.top.val[1,1], top.uni.est.2 = uni.trial2[1, uni.top.names[1]],
              top.lasso.est.1 = uni.top.val[1,lasso.top.names[1]], top.lasso.est.2 = uni.trial2[1, lasso.top.names[1]],
              sig.lasso.est.1 = sig.lasso.est.1, sig.lasso.est.2 = sig.lasso.est.2,
              true.est.1 = uni.top.val[1, causal.names],true.est.2=uni.trial2[1, causal.names],
              lasso.pred.1 = lasso.pred.trial1,
              lasso.pred.2 = lasso.pred.trial2,
              uni.pred.1 = uni.pred.1, uni.pred.2 = uni.pred.2,
              lasso.err.1 = lasso.pred.trial1 - data.trial1$outcome,
              lasso.err.2 = lasso.pred.trial2 - data.trial2$outcome,
              uni.err.1 = uni.pred.1-data.trial1$outcome, uni.err.2 = uni.pred.2 - data.trial2$outcome
              )
  
  
}

#' bootstrap or two fold CV

boot.cv <-  function(x, x.names, response.type, causal.names, topn = 5, n.rep = 100, replace = TRUE){
  name.mat <- data.matrix(sapply(1:n.rep, function(i)sample(1:nrow(x), nrow(x), replace = replace)))
  if(!replace) name.mat <- name.mat[1:round(nrow(x)/2),,drop=F] # if cv, only take the first half of patients
  cv.res <- sapply(1:n.rep, function(i){
    data.trial1 <- x[name.mat[,i],]
    data.trial2 <- x[setdiff(1:nrow(x), name.mat[,i]),]
    res <-  run.multi(data.trial1=data.trial1, data.trial2=data.trial2, x.names=x.names, 
                      response.type=response.type, causal.names=causal.names, topn=topn)
    out <- res[c("lasso.err.1", "lasso.err.2","uni.err.1", "uni.err.2")]
  }, simplify=F)
  
  
}