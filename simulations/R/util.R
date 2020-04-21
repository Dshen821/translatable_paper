#' Simulates Trial Data
#'
#' Simulates clinical trial data, given user specifications of sample size, number of causal features,
#' number of correlated causal features, noise features. Causual feature correlation structure,
#' distribution shifts can be specified. 
#' @param n.trial number of samples in each trial. should be a vector. e.g. c(200, 800) means 200 samples in trial 1
#' and 800 samples in future
#' @param n.causal number of independent causal features to generate in each trial. they will be named as c("a.0","b.0"). 
#' all trials will be generated with the same causal factors
#' @param n.cor.causal number of correlated features to generate, for each causal group. e.g. c(100, 200).
#' feature names will be a.1, a.2....a.100, b.1,...b.200. where a.1...a.100 are correlated to a.0 (causal group a), and 
#' b.1 ... b.200 are correlated to b.0 (causal group b)
#' @param cor.causal correlation coefficient within each "causal" group. 
#' If it is a vector, assuming all trials have the same correlation structure. e.g. c(.9, .2) indicates pairwise 
#' correlation within causal group 1 is .9 for all trials
#' and correlation within causal group 2 is .2 for all trials.
#' If trials are assumed to have different correlation structures, it should be a matrix with n.trial columns, n.causal rows. 
#' where column 1 row 1 indicates correlation of first causal group in legacy, column 2 row 1 indiactes correlation of the first causal group in trial 2. 
#' @param coef.causal coefficient between causal features and outcome (beta in Y = betaX). 
#' If it is a vector, assuming all trials follow the same coefficients. e.g. c(1, .2) indicates coef of causal group 1 is 1 for all trials
#' and coef of causal group 2 is .2 for all trials.
#' If trials are assumed to have different coefficients, it should be a matrix with n.trial columns, n.causal rows. 
#' where column 1 row 1 indicates coef of first causal group in legacy, column 2 row 1 indiactes coef of the first causal group in trial 2. 
#' @param n.noise number of noise features (Xs) to generate
#' @param outcome.sd when gaussian outcome is generated, noise follows N(0, eps)
#' @param shift.mean,shift.sd set to non-zero if want to simulate trials with population shifts. should be a vector with length of n.causal.
#' e.g. shift.mean = c(.2,.3) and shift.sd = c(.1, .15) indicates a.0 difference between any two trials follows N(.2, .1) and b.0 difference
#' between any two trials follows N(.3,.15)
#' 
#' @return List containing dataframe of trial data, covariate names, and causal variale name.
simu <- function(n.trial, 
                 n.causal = 2, n.cor.causal = c(10,20), cor.causal = c(.9, .8), coef.causal=c(1,0),
                 n.noise=50, outcome.sd=.1,
                 shift.mean = c(0, 0), shift.sd=c(.1, .1)) {
  
  nsamples <- sum(n.trial) # total patient= samples trial1 + samples trial 2 + ... + trial n
  nt <- length(n.trial) # number of trials
  coef.causal <- matrix(coef.causal, nrow=n.causal) # coefficients specified must match number specified
  cor.causal <- matrix(cor.causal, nrow=n.causal)
  if(n.causal>0)stopifnot(n.causal == length(n.cor.causal)) else stopifnot(is.null(n.cor.causal))
  stopifnot(n.causal+n.noise>0) # independent features + noise features > 0
  
  
  # BASIC DATA: simulate causal features and noise feactures
  names.no.cor <- causal.names <- NULL
  if(n.causal>0) {
    causal.names <- paste0(letters[1:n.causal],".0") # create name of causal features: a0, b0, ... ,
    names.no.cor <- c(names.no.cor, causal.names) # causal are a.0, b.0...
  }

  # CREATE NOISE
  if(n.noise>0) names.no.cor <- c(names.no.cor, paste0("noise.",1:n.noise-1)) # noises are noise.0, noise.1, ....
  x.no.cor <- sapply(1: (n.causal+n.noise), function(i)rnorm(nsamples)) # creates nsamples x n.causal + n.noise rnormals
  colnames(x.no.cor) <- names.no.cor
  x.no.cor <- as.data.frame(x.no.cor)
  
  # labels trial1, trial1,...,trial2, trial2,...trial3,trialnt, according to index count in n.trial[i]
  x.no.cor <- x.no.cor %>% 
    mutate(trial = unlist(sapply(1:nt, function(i)rep(paste0("trial",i), n.trial[i])))) # col for trial#
  
  
  # SCENARIO 1:  population shift
  # for each a0, a1,.... causal effect, within n.trials we randomly shift it by none, shift, 2shift, ...
  if(unique(shift.mean)!=0){ #when shift.mean is not all 0's
    for (fe in 1:n.causal){
      if(shift.mean[fe]>0){
        x.var <- paste0(letters[fe],".0")
        random.index <- sample(0:(nt-1), nt) # randomly asisgn n.trial trials to mean, mean+shift, mean+2shift, ... 
        
        for(i in 1:nt){
          x.no.cor[which(x.no.cor$trial==paste0('trial',i)), x.var] <- x.no.cor[which(x.no.cor$trial==paste0('trial',i)), x.var] + 
            rnorm(n.trial[i], shift.mean[fe]*random.index[i], shift.sd[fe])
            # rnorm(number in trial, shift.mean amount * random index (0,1,2,3), shift.sd)
        }
      }
      }
  }
  
  #  SCENARIO 2: Simulate correlated features: allow different correlation matrix in trial 1 and 2
  
  ## Creates a.1,....,a.10 correlated features to a.0 and binds them to x.no.cor
  x <- x.no.cor
  if(n.causal>0){
    for (fe in 1:n.causal){
      if(n.cor.causal[fe]>0){
        prefix <- letters[fe]
        cor.causal.use <- cor.causal[fe, ] # cor to use for each causal feature
        if(ncol(cor.causal)!=1)cor.v <- unlist(sapply(1:nt, function(i)rep(cor.causal.use[i], n.trial[i]))) # different cor c(0.98, 0.75) per trial1, trial2,..
        else cor.v <- rep(cor.causal.use, sum(n.trial)) # specify same cor for each obs
        
        # cor for each patient (all pts have same cor structure, or each trial has diff cor structure)
        # takes x.no.cor value of a.0, b.0,...  *sqrt(cor.v) + rnorm(nsamples ) * sqrt(1-cor.v) -- what eq is this ?
        x.cor.tmp <- sapply(1:n.cor.causal[fe], function(i) x.no.cor[, paste0(prefix,".0")] * sqrt(cor.v) + rnorm(nsamples) * sqrt(1-cor.v))
        colnames(x.cor.tmp) <- paste0(prefix, ".", 1:ncol(x.cor.tmp))
        
        x <- cbind(x, data.frame(x.cor.tmp))
        
        #x now has "a.1" "a.2"      "a.3"      "a.4"      "a.5"      "a.6"      "a.7"      "a.8"      "a.9"      "a.10" 
        
      }
    }
    
  }
  
  # Simulate outcome
  # nothing causal, so outcome is just rnorm vars
  if(n.causal==0) x$outcome <- rnorm(nsamples)
  # for each trial, take x[trialn, c(a.0, b.0,...)] * coef of (a.0, b.0,..), transpose and row sum
  # add rnorm(length of trialn, 0, outcome.sd) and add as noise 
  # gives outcome
  if(n.causal>0){ 
    x$outcome <- NULL
    for(i in 1:nt){
    coef.use <- coef.causal[,ifelse(ncol(coef.causal)==1, 1, i)]
    tmp <-which(x$trial==paste0("trial",i))
    x$outcome[tmp] <- colSums(t(x[tmp, causal.names]) * coef.use) + 
      rnorm(length(tmp), 0, outcome.sd)
  }
  }
  
  # order cols, label rows
  x <- x[order(names(x))]
  rownames(x) <- paste0('s.', 1:nrow(x))
  x.names <- setdiff(names(x), c("outcome","trial")) # appears in x, and not called outcome, trial
  out <- list(data = x, x.names=x.names, causal.names=causal.names)
}

#' CROSS-Validated LASSO Regression
#' Wrapper to run lasso penalized regression. Lambda chosen by CV to select minimum MSE lambda.
#' @inheritParams glmnet
#' @inheritParams cv.glmnet
#' @param top Numeric specifying "Top #" covariate chosen by lasso. Top is ordered by largest to smallest coefficient.
#' 
#' @return model fit and nonzero.covariates.
run.lasso <- function(x,y, alpha=1,top=NULL, family="binomial", ...){
  # filter out any rows of x that have missing
  nona <- which(!is.na(rowMeans(x))) # indexes non missing rows
  # CV determines best lambda for lasso
  lasso.cv<- cv.glmnet(x=x[nona,], y=y[nona],alpha=alpha,family=family,...)
  lasso.pen <- lasso.cv$lambda.min #optimal lambda; lasso.pen #minimal shrinkage
  
  #fit lasso lm based on best lambda
  lasso.fit <-glmnet(x = x[nona,], y=y[nona],alpha=alpha,family=family,lambda = lasso.pen,...) #estimate the
  lasso.coef <- as.vector(data.matrix(coef(lasso.fit)))
  names(lasso.coef) <- coef(lasso.fit)@Dimnames[[1]]
  lasso.coef <- lasso.coef[-1]# drop intercpt
  # Sort lasso coefficeints by magnitude
  lasso.coef.od <- order(abs(lasso.coef),decreasing=T)
  lasso.coef.sort <- lasso.coef[lasso.coef.od]
  
  # Return nonzeros (all or top # specified)
  if(!is.null(top))nznames <- lasso.coef.sort[which(lasso.coef.sort!=0)][1:min(top, sum(lasso.coef.sort!=0))] else nznames <- lasso.coef.sort[which(lasso.coef.sort!=0)]
  
  out <- list(fit=lasso.fit, nonzero.names = nznames) # lassofit and non zero coefs
  return(out)
}

# covs <- simu.1[[1]]$x.names
# run.lasso(x = data.matrix(simu.1[[1]]$data[, covs]), y = simu.1[[1]]$data[, "outcome"], family = "gaussian")



#' Univariate screening procedure.
#' 
#' @inheritParams glmnet
#' @inheritParams cv.glmnet
#' @param save.model Logical includes model and predictions along with pvalues
#' 
#' @return Returns covariates ordered by smallest p-value to largest. Can request model and predictions as well.
run.glm.uni <- function(x, y, family="binomial", save.model=FALSE){
  x <- data.frame(x)

  if(!save.model){ # singly regress y on each x, storing estimate and pval (row1 and row2) by cols (covariates)
    vals <- sapply(colnames(x), function(i)coef(summary(glm(y~x[,i], family=family)))[2,c("Estimate","Pr(>|t|)")])
    vals.sort <- vals[,order(vals["Pr(>|t|)",])] # smallest to largest based on pval (increasing)
    out <- vals.sort # just p value small to alrge
  }
  if(save.model){
    mod <- sapply(colnames(x), function(i){
      tmp <- data.frame(y0 = y, x0= x[,i]) # y, x[,i] data frame
      glm(y0 ~ x0, data=tmp, family=family)} # stores all model fits in list
      , simplify = F)
    # Store estimate and pval, ordering from smallest p val to greatest
    vals <- data.matrix(sapply(mod, function(i)coef(summary(i))[2,c("Estimate","Pr(>|t|)")]))
    vals.sort <- vals[,order(vals["Pr(>|t|)",])]
    pred <- data.matrix(sapply(mod, function(i)predict(i))) # use each model to predit on y
    out <- list(model=mod, vals.sort=vals.sort, pred=pred) # return model, pvals, prediciton
  }
  out
}

#run.glm.uni(x = simu.1[[1]]$data$a.0, y = simu.1[[1]]$data$outcome, family = "gaussian")

#' Run's LASSO and Uni-variate screening procedures. Returns Estimates, and calculates
#' MSE (Test, Train, Future)
#'
#' @param data.legacy
#' @param data.future
#' @param x.names names of the features to use
#' @param response.type e.g. gaussian
#' @param causal.names 
#' @param top5n
#' 
#' @return List of covariate names, estimates, predictions, MSE
run.multi <- function(data.legacy, data.future, x.names, response.type, causal.names, topn=5){
  
  lasso.out <- run.lasso(x = data.matrix(data.legacy[x.names]), y = data.legacy$outcome, family=response.type)
  lasso.top <- lasso.out$nonzero.names
  lasso.top.names <- names(lasso.top)
  if(length(lasso.top.names)==0) lasso.top.names <- lasso.top <-  NA
  uni.top <- run.glm.uni(x = data.matrix(data.legacy[x.names]), y = data.legacy$outcome, family=response.type, save.model = TRUE)
  uni.top.val <- uni.top$vals.sort
  uni.top.names <-  colnames(uni.top$vals.sort)
  uni.future <- run.glm.uni(x = data.matrix(data.future[uni.top.names]), y = data.future$outcome, family=response.type)
  uni.pred.1 <- uni.top$pred[,uni.top.names]
  uni.pred.2 <- sapply(1:length(uni.top.names), function(i){
    tmp <- data.frame(y0 = data.future$outcome, x0= data.future[[uni.top.names[i]]]) 
    predict(uni.top$mod[[uni.top.names[i]]], newdata = tmp)
    })
  colnames(uni.pred.2) <- uni.top.names

  lasso.pred.legacy <- predict(newx = data.matrix(data.legacy[x.names]), object=lasso.out$fit)
  sig.lasso.est.1 <- run.glm.uni(x = data.matrix(lasso.pred.legacy), y = data.legacy$outcome, family=response.type)["Estimate"]
  lasso.pred.future <- predict(newx = data.matrix(data.future[x.names]), object=lasso.out$fit)
  sig.lasso.est.2 <- run.glm.uni(x = data.matrix(lasso.pred.future), y = data.future$outcome, family=response.type)["Estimate"]
  
  
  out <- list(lasso.top.names = lasso.top.names, uni.top.names = uni.top.names[1:topn], 
              top.uni.est.1 = uni.top.val[1,1], top.uni.est.2 = uni.future[1, uni.top.names[1]],
              top.lasso.est.1 = uni.top.val[1,lasso.top.names[1]], top.lasso.est.2 = uni.future[1, lasso.top.names[1]],
              sig.lasso.est.1 = sig.lasso.est.1, sig.lasso.est.2 = sig.lasso.est.2,
              true.est.1 = uni.top.val[1, causal.names],true.est.2=uni.future[1, causal.names],
              lasso.pred.1 = lasso.pred.legacy,
              lasso.pred.2 = lasso.pred.future,
              uni.pred.1 = uni.pred.1, uni.pred.2 = uni.pred.2,
              lasso.err.1 = lasso.pred.legacy - data.legacy$outcome,
              lasso.err.2 = lasso.pred.future - data.future$outcome,
              uni.err.1 = uni.pred.1-data.legacy$outcome, uni.err.2 = uni.pred.2 - data.future$outcome
              )
}

#' Two fold CV , pre-defined CV, or Bootstrap
#' 
#' @param name.mat pre defined training Index for CV. name.mat should be a list. each element is for one CV run. 
#' Each column represents one repeat. 
#' Rows are index of data entries to in included in training set. (other entries will be automatically used as testing).
#' if name.mat is not NULL, n.rep and replace will be ignored
#' 
#' @n.rep Times to repeat Bootstrap or CV
#' @replace logical if TRUE then Bootstrap, else CV( 2-fold)
#' @inheritParams run.multi

boot.cv <-  function(x, x.names, response.type, causal.names, topn = 5, n.rep = 100, replace = TRUE, name.mat = NULL){
  if(is.null(name.mat)){ # if replace == TRUE, n.temp = nrow(x) else nrow(x) /2
    n.tmp <- ifelse(replace, nrow(x), round(nrow(x)/2)) # bootstrap vs CV
    name.mat <- sapply(1:n.rep, function(i)sample(1:nrow(x), nrow(x), replace = replace)[1:n.tmp], simplify=F)
    # if cv, only take the first half of patients as training
  }
  
  # Determine legacy and future set and runs run.multi
    cv.res <- sapply(1:length(name.mat), function(i){
      index <- name.mat[[i]]
      data.legacy <- x[index,]
      data.future <- x[setdiff(1:nrow(x), index),]
      res <-  run.multi(data.legacy=data.legacy, data.future=data.future, x.names=x.names, 
                      response.type=response.type, causal.names=causal.names, topn=topn)
    out <- res[c("lasso.err.1", "lasso.err.2","uni.err.1", "uni.err.2")]
  }, simplify=F)
    
    return(cv.res)
}



#' heatmap for visualization
heat.fun <- function(simu.out){
  data.tmp <- simu.out$data 
  tmp <- data.tmp %>% select(outcome, trial)
  pheatmap(t(data.tmp[,simu.out$x.names]),
           cluster_cols=F, cluster_rows = F,
           annotation_row = NA,
           annotation_col=tmp)
}

#' summarize empirical correlation estimates and shifts from simulated data
summary.cor <- function(simu.out){
  x.names <- simu.out$x.names[which(substr(simu.out$x.names,1,2)=="a.")]
  xcor.list <- sapply(unique(data$trial), function(i)
    cor(simu.out$data %>% filter(trial==i) %>% select(!!x.names, outcome)))
  x.shift <- simu.out$data %>% group_by(trial) %>% summarize(a0_mean=mean(a.0))
  model.list <- sapply(unique(data$trial), function(i) 
    coef(summary(lm(outcome~a.0, data=simu.out$data %>% filter(trial==i)))))
  out <- list(xcor.list=xcor.list, x.shift=x.shift,
              model.list=model.list)
}









#' run uni variate screening using lm
#' NOT USED FOR NOW  - run.glm.uni can cover this functionality
run.gaussian.uni <- function(x, y){
  vals <- sapply(colnames(x), function(i)coef(summary(lm(y~x[,i])))[2,c("Estimate","Pr(>|t|)")])
  vals.sort <- vals[,order(vals["Pr(>|t|)",])]
  vals.sort
}



# RF
# rf <- randomForest(outcome ~ ., data = simu.1[[1]]$data[, c(covs, "outcome")], importance = TRUE, ntree=1000)
# which.min(rf$mse)
# 
# imp <- as.data.frame(sort(importance(rf)[,1],decreasing = TRUE),optional = T)
# names(imp) <- "% Inc MSE"
# imp
# varImpPlot(rffit) 
# ###
