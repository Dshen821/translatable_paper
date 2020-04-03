library("glmnet")
data(QuickStartExample)

#' lasso wrapper
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
  
  out <- list(fit=lasso.fit, nonzero.names = nznames)
  return(out)
}


# 
# run.lasso(x,y, family = "gaussian")
# 
# # filter out any rows of x that have missing
# nonas <- which(!is.na(rowMeans(x)))
# # CV determines best lambda for lasso
# lasso.cvd <- cv.glmnet(x=x[nonas,], y=y[nonas],alpha=1,family="gaussian") # cross validated lasso via MSE
# lasso.opt <- lasso.cvd$lambda.min #optimal lambda that minimizes mse
# 
# #fit lasso lm based on best lambda
# lasso.fit <-glmnet(x = x[nonas,], y=y[nonas],alpha=1,family="gaussian",lambda = lasso.opt) 
# lasso.coef <- as.vector(data.matrix(coef(lasso.fit))) # return coeffs
# names(lasso.coef) <- coef(lasso.fit)@Dimnames[[1]]
# lasso.coef <- lasso.coef[-1] # drop int
# 
# 
# # Sort lasso coefficeints by magnitude
# lasso.coef.od <- order(abs(lasso.coef),decreasing=T) # index from largest to smallest
# lasso.coef.sort <- lasso.coef[lasso.coef.od]
# 
# # Return nonzeros
# if(!is.null(top)){
#   nznames <- lasso.coef.sort[which(lasso.coef.sort!=0)][1:min(top, sum(lasso.coef.sort!=0))] #returns nonzero coefficients limited by top # largest coefficeints
#   }else{
#     nznames <- lasso.coef.sort[which(lasso.coef.sort!=0)] # returns nonzero coefficeitns
#   }
# out <- list(fit=lasso.fit, nonzero.names = nznames) # returns the lasso.fit, and non coef names
# 
# 







#' run uni-variate screening
run.glm.uni <- function(x, y, family="binomial", save.model=FALSE){
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

# colnames(x) <- paste0(rep("cov", 20), 1:20)
# 
# 
# tmp <- data.frame(y0 = y, x0 = x[, 1])
# 
# 
# glm(y0 ~ x0, data=tmp, family="gaussian")
# # store linear reg fit of each x in a list
# mod <- sapply(colnames(x), function(i){
#   tmp <- data.frame(y0 = y, x0= x[,i]) # y, x[,i] data frame
#   glm(y0 ~ x0, data=tmp, family="gaussian")}, simplify = F)
# 
# vals <- data.matrix(sapply(mod, function(i)coef(summary(i))[2,c("Estimate","Pr(>|t|)")]))
# vals.sort <- vals[,order(vals["Pr(>|t|)",])] # order by increasing p value
# pred <- data.matrix(sapply(mod, function(i)predict(i))) # predict on y for each cov
# out <- list(model=mod, vals.sort=vals.sort, pred=pred)




run.multi <- function(data.legacy, data.future, x.names, response.type, causal.names, topn=5){
  
  lasso.out <- run.lasso(x = data.matrix(data.legacy[x.names]), y = data.legacy$outcome, family=response.type)
  lasso.top <- lasso.out$nonzero.names
  lasso.top.names <- names(lasso.top)
  if(length(lasso.top.names)==0) lasso.top.names <- lasso.top <-  NA
  uni.top <- run.glm.uni(x = data.matrix(data.legacy[x.names]), y = data.legacy$outcome, family=response.type, save.model = TRUE)
  uni.top.val <- uni.top$vals.sort
  uni.top.names <-  colnames(uni.top$vals.sort)
  uni.future <- run.glm.uni(x = data.matrix(data.future[uni.top.names]), y = data.future$outcome, family=response.type)
  uni.pred.1 <- uni.top$pred[,uni.top.names] #
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
              uni.err.1 = uni.pred.1-data.legacy$outcome,
              uni.err.2 = uni.pred.2 - data.future$outcome
  )
}

data.legacy <- cbind(data.frame(x), data.frame(outcome = y))
x.names <- colnames(x)

ind <- sample(1:nrow(data.legacy), 80, replace = FALSE)

data.legacy <- data.legacy[ind,]
data.future <- data.legacy[-ind,]

# lasso wrapper on y ~ all x's.
lasso.out <- run.lasso(x = data.matrix(data.legacy[x.names]), y = data.legacy$outcome, family="gaussian")
lasso.top <- lasso.out$nonzero.names #top covariates
lasso.top.names <- names(lasso.top) #top cov names
if(length(lasso.top.names)==0) lasso.top.names <- lasso.top <-  NA

# univariate screening wrapper: all y ~ each x singlely
uni.top <- run.glm.uni(x = data.matrix(data.legacy[x.names]), y = data.legacy$outcome, family=response.type, save.model = TRUE) #  return model, pvals, prediction
uni.top.val <- uni.top$vals.sort #estimates, and pvalues small to large 
uni.top.names <-  colnames(uni.top$vals.sort) # names of covariates

# univariate screening Future data: all y ~ each of the top x's above singlely
uni.future <- run.glm.uni(x = data.matrix(data.future[uni.top.names]), y = data.future$outcome, family=response.type)
uni.pred.1 <- uni.top$pred[,uni.top.names] # Using legacy model predictions on legacy daya
uni.pred.2 <- sapply(1:length(uni.top.names), function(i){
  tmp <- data.frame(y0 = data.future$outcome, x0= data.future[[uni.top.names[i]]])  # data.frame of outcome and one top x
  predict(uni.top$mod[[uni.top.names[i]]], newdata = tmp) # Using legacy model predict on future data
})
colnames(uni.pred.2) <- uni.top.names


# Lasso predictions on legacy data
lasso.pred.legacy <- predict(newx = data.matrix(data.legacy[x.names]), object=lasso.out$fit)
# legacy outcome ~ lasso prediction on legacy
sig.lasso.est.1 <- run.glm.uni(x = data.matrix(lasso.pred.legacy), y = data.legacy$outcome, family="gaussian")["Estimate"]

#Using lasso fit on legacy to predict on future data
lasso.pred.future <- predict(newx = data.matrix(data.future[x.names]), object=lasso.out$fit)
# lasso future predicted regressed: future y ~ pred y
sig.lasso.est.2 <- run.glm.uni(x = data.matrix(lasso.pred.future), y = data.future$outcome, family="gaussian")["Estimate"]



topn <- 3

out <- list(lasso.top.names = lasso.top.names, uni.top.names = uni.top.names[1:topn], 
            top.uni.est.1 = uni.top.val[1,1],                  # COV estimate
            top.uni.est.2 = uni.future[1, uni.top.names[1]],     #uni.top.names based on lowest p value
            top.lasso.est.1 = uni.top.val[1,lasso.top.names[1]], 
            top.lasso.est.2 = uni.future[1, lasso.top.names[1]],
            sig.lasso.est.1 = sig.lasso.est.1,
            sig.lasso.est.2 = sig.lasso.est.2,
            true.est.1 = uni.top.val[1, causal.names],
            true.est.2=uni.future[1, causal.names],
            lasso.pred.1 = lasso.pred.legacy,
            lasso.pred.2 = lasso.pred.future,
            uni.pred.1 = uni.pred.1, uni.pred.2 = uni.pred.2,
            lasso.err.1 = lasso.pred.legacy - data.legacy$outcome, #LASSO LEGACY - ACTUAL LEGACY = TRAINING E
            lasso.err.2 = lasso.pred.future - data.future$outcome, #LASSO FUTURE - ACTUAL FUTURE = TESTING E
            uni.err.1 = uni.pred.1 - data.legacy$outcome, # UNI LEGACY - ACTUAL LEGACY  = TRAINING E      
            uni.err.2 = uni.pred.2 - data.future$outcome  # UNI FUTURE - ACTUAL FUTURE = TESTING E
)




#' Two fold CV , pre-defined CV, or Bootstrap
#' @param name.mat pre difined training Index for CV. name.mat should be a list. each element is for one CV run. 
#' Each column represents one repeat. 
#' Rows are index of data entries to in included in training set. (other entries will be automatically used as testing).
#' if name.mat is not NULL, n.rep and replace will be ignored
#' @inheritParams run.multi

boot.cv <-  function(x, x.names, response.type, causal.names, topn = 5, n.rep = 100, replace = TRUE, name.mat = NULL){
  if(is.null(name.mat)){
    n.tmp <- ifelse(replace, nrow(x), round(nrow(x)/2)) # bootstrap vs CV
    name.mat <- sapply(1:n.rep, function(i)sample(1:nrow(x), nrow(x), replace = replace)[1:n.tmp], simplify=F)
    # if cv, only take the first half of patients as training
  }
  cv.res <- sapply(1:length(name.mat), function(i){
    index <- name.mat[[i]]
    data.legacy <- x[index,]
    data.future <- x[setdiff(1:nrow(x), index),]
    res <-  run.multi(data.legacy=data.legacy, data.future=data.future, x.names=x.names, 
                      response.type=response.type, causal.names=causal.names, topn=topn)
    out <- res[c("lasso.err.1", "lasso.err.2","uni.err.1", "uni.err.2")]
  }, simplify=F)
  
  
}





boot.cv(x=tmpdata, x.names=x.names,
        response.type = response.type, causal.names = causal.names, topn = 5, n.rep = n.rep, replace = FALSE)