#Estimate W using four-fold cross validation to find the optimal tuning parameter
optimalW_dantzig <- function(y, X, tune.par.seq){
  n <- length(y)
  sample.id <- sample(1:4, size=n, prob = rep(0.25, 4), replace=T)
  y.train <- lapply(1:4, function(m) y[sample.id!=m])
  y.test <- lapply(1:4, function(m) y[sample.id==m])
  X.train <- lapply(1:4, function(m) X[sample.id!=m,])
  X.test <- lapply(1:4, function(m) X[sample.id==m,])

  tune.par.slim <- slim(Y=y, X=X, method="dantzig", nlambda=length(tune.par.seq))$lambda

  par.grid <- expand.grid(1:4, tune.par.slim)
  names(par.grid) <- c("sample.id","tune.par.slim")
  W.hat.all <- apply(par.grid, 1, function(m)
    slim(Y=y.train[[as.numeric(m)[1]]], X=X.train[[as.numeric(m)[1]]], method="dantzig",lambda = as.numeric(m)[2])$beta)
  par.grid <- cbind(par.grid, id.expand = 1:nrow(par.grid))
  pred.error <- apply(par.grid, 1, function(m)
    sum(y.test[[as.numeric(m)[1]]] - X.test[[as.numeric(m)[1]]]%*%W.hat.all[,as.numeric(m[3])])^2)

  pred.error.mean <- tapply(pred.error, par.grid$tune.par.slim, function(x) mean(x,na.rm = T))
  tune.par.opt <- as.numeric(names(which.min(pred.error.mean)))
  W.hat <- slim(Y=y, X=X, method="dantzig",lambda = tune.par.opt)$beta
  return(list(W.hat = W.hat, tune.par.opt = tune.par.opt))
}
