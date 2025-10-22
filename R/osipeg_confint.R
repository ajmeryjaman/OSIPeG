#' Confidence intervals for one-step improved penalized G-estimator (OSIPeG)
#' @description This function implements a decorrelated score approach to provide confidence intervals
#' for the blip coefficients selected by penalized G-estimation in the context of a structural nested mean
#' model (SNMM) for repeated outcomes.
#' @param data A data frame containing the variables in longitudinal format. In the data, the
#' outcome should be continuous and the treatment/exposure should be binary. The data must have 
#' an additional column containing the propensity scores.
#' @param wc.str A character string specifying the working correlation structure. The
#' following are currently allowed: "independence", "exchangeable", "ar1", and "unstructured".
#' @param id.var The column name in data that corresponds to the variable id (unique identifier).
#' @param response.var The column name in data that corresponds to the response variable.
#' @param treat.var The column name in data that corresponds to the treatment/exposure variable.
#' @param tf.model A single formula object specifying the covariates of a (linear)
#' treatment-free model.
#' @param estimate The vector of blip parameter estimates obtained using penalized G-estimation.
#' @param alpha.hat The estimated correlation parameter(s) alpha(s) if the provided structure is either
#' "exchangeable", "ar1, or "unstructured". For unstructured, the elements of alpha.hat correspond
#' to the upper triangular portion of working correlation matrix having dimension equal to
#' the largest cluster size.
#' @param sigma2.hat The estimated variance parameter sigma^2.
#' @param lambda The optimal tuning parameter returned by penalized G-estimation.
#' @param test.size The significance level. For a 95\% confidence interval test.size should be 0.05.
#' @param tune.par.seq A sequence of values (in decreasing order) for the tuning parameter to execute
#' the dantzig selector. This is used to find the optimal weight vector.
#' @param continuous.covs A logical vector of TRUE/FALSE values identifying the potential
#' effect modifiers that are continuous.
#' @return A list containing the following:
#' \item{psi.hat.one.step.full}{A vector showing the one-step improved penalized G-estimates for the selected blip
#' coefficients, calculated using the full weight vector.}
#' \item{psi.hat.one.step.LASSO}{A vector showing the one-step improved penalized G-estimates for the selected blip
#' coefficients, using the sparse weight vector obtained via LASSO.}
#' \item{psi.hat.one.step.dantzig}{A vector showing the one-step improved penalized G-estimates for the selected blip
#' coefficients, using the sparse weight vector obtained via the Datzig selector.}
#' \item{CI.one.step.full}{A matrix showing confidence interval estimates for the selected blip
#' coefficients, calculated using the full weight vector.}
#' \item{CI.one.step.LASSO}{A matrix showing confidence interval estimates for the selected blip
#' coefficients, using the sparse weight vector obtained via LASSO.}
#' \item{CI.one.step.dantzig}{A matrix showing confidence interval estimates for the selected blip
#' coefficients, using the sparse weight vector obtained via the Datzig selector.}
#' @export
#'
#' @examples
#' library(mvtnorm)
#'
#' expit <- function(x) exp(x)/(1+exp(x))
#'
#' ## data.gen is a function that generates a longitudinal data set for a specific correlation
#' ## structure. Available structures in this function are: independence, exchangeable and AR1.
#'
#' ## Arguments(data.gen):
#' #     n = Number of subjects
#' #     ni = A vector containing number of time points for each subject
#' #     sigma2.e = Error variance
#' #     alpha = Correlation parameter
#' #     corstr = The correlation structure among the repeated outcomes
#' #     autocorr.coef = The autocorrelation coefficient for inducing correlation among the
#' #     continuous confounders and the noise covariates
#'
#' data.gen <- function(n, ni, sigma2.e, alpha, corstr, autocorr.coef){
#'  ncovs  <-  2+4+45 # 2+4=6 confounders and 45 noise covariates
#'   beta <- c(0, 1, -1.1, 1.2, 0.75, -0.9, 1.2) # treatment model parameters
#'   delta <- c(1, 1, 1.2, 1.2, -0.9, 0.8, -1,
#'              rep(1, 20), rep(0, ncovs-6-20),
#'              -0.8, 1, 1.2, -1.5) # treatment-free model parameters
#'   psi <- c(1, 1, -1, -0.9, 0.8, 1, 0, rep(0, 20), rep(0, ncovs-6-20)) # blip parameters
#'
#'   # generating two continuous baseline covariates
#'   l1 <- rnorm(n, 0, 1)
#'   l2 <- rnorm(n, 0, 1)
#'
#'   # V is the covariance matrix of the time-varying confounders (l3,..,l6) and
#'   # noise covariates (x1,...)
#'   V <- toeplitz(autocorr.coef^(0:(ncovs-2-1)))
#'
#'   lx <- a <- y <- vector(mode="list", length=n)
#'   lx.mat <- NULL
#'   for(i in 1:n){
#'     a[[i]] <- y[[i]] <- rep(NA, ni[i])
#'     lx[[i]] <- matrix(NA, ni[i], ncovs)
#'     lx[[i]][,1] <- rep(l1[i], ni[i])
#'     lx[[i]][,2] <- rep(l2[i], ni[i])
#'
#'     corr.mat <- switch(corstr,
#'                        "exchangeable" = toeplitz(c(1, rep(alpha, ni[i]-1))),
#'                        "ar1" = toeplitz(alpha^(0:(ni[i]-1))),
#'                        "independence" = diag(ni[i])
#'     )
#'     cov.mat <- diag(sqrt(sigma2.e), ni[i]) %*% corr.mat %*% diag(sqrt(sigma2.e), ni[i])
#'     e <- rmvnorm(1, sigma = cov.mat)
#'
#'     # j=1
#'     mu.lx <- c(NA, NA, # for l1 and l2
#'                rep(0, 4), rep(0, ncovs-6)) #rep(0.3*lx[[i]][1, 1]+0.3*lx[[i]][1, 2],4)
#'     lx[[i]][1,3:ncovs] <- rmvnorm(1, mean=mu.lx[3:ncovs], sigma = V)
#'     a[[i]][1] <- rbinom(1, 1, expit(sum(c(1,lx[[i]][1,1:6])*beta)))  ## no correlation
#'     tf.mean <- sum(c(1, lx[[i]][1,1:ncovs],
#'                      lx[[i]][1,1]*lx[[i]][1,5], lx[[i]][1,3]*lx[[i]][1,4], sin(lx[[i]][1,3]-lx[[i]][1,4]),
#'                      cos(2*lx[[i]][1,5])) * delta)
#'     blip <- sum(c(a[[i]][1], a[[i]][1]*c(lx[[i]][1,1:ncovs])) * psi)
#'     y[[i]][1] <- (tf.mean + blip) + e[1]
#'
#'
#'     # j=2:ni
#'     for(j in 2:ni[i]){
#'       mu.lx <- c(NA, NA, # for l1 and l2
#'                  0.3*lx[[i]][j-1, 3:6] + 0.3*a[[i]][j-1], 0.5*lx[[i]][j-1,7:ncovs])
#'       lx[[i]][j,3:ncovs] <- rmvnorm(1, mean=mu.lx[3:ncovs], sigma = V)
#'       a[[i]][j] <- rbinom(1, 1, expit(sum(c(1,lx[[i]][j,1:6])*beta)))  ## no correlation
#'       tf.mean <- sum(c(1, lx[[i]][j,1:ncovs],
#'                        lx[[i]][j,1]*lx[[i]][j,5], lx[[i]][j,3]*lx[[i]][j,4], sin(lx[[i]][j,3]-lx[[i]][j,4]),
#'                        cos(2*lx[[i]][j,5])) * delta)
#'       blip <- sum(c(a[[i]][j], a[[i]][j]*c(lx[[i]][j,1:ncovs])) * psi)
#'       y[[i]][j] <- (tf.mean + blip) + e[j]
#'     }
#'     lx.mat <- rbind(lx.mat, lx[[i]])
#'   }
#'
#'   colnames(lx.mat) <- c(paste("l", 1:6, sep=""), paste("x", 1:(ncovs-6), sep=""))
#'   data <- data.frame(id=rep(1:n, times=ni), a=unlist(a), lx.mat, y=round(unlist(y),3))
#'   return(data)
#' }
#'
#' data.s <- data.gen(n = 500, ni = rep(6, 500), sigma2.e = 1, alpha = 0.8,
#'                    corstr = "exchangeable", autocorr.coef = 0.25)
#'
#' ncovs  <-  2+4+45
#' # treatment-free model is misspecified
#' tf.model <- as.formula(paste("~",paste(c(paste("l",1:6,sep=""), paste("x",c(1:9,11:(ncovs-6)),sep="")),
#'                                        collapse = "+"), collapse=""))
#' treat.model <- ~l1+l2+l3+l4+l5+l6
#'
#' lam_max <- 1
#' lam_min <- 0.01*lam_max
#' lambda.seq <- sort(seq(lam_min,lam_max,(lam_max-lam_min)/99), decreasing=TRUE)
#'
#' # library(devtools) # if already installed, otherwise need to install it first
#' # install_github("ajmeryjaman/penalizedG") # run this command if this package is not already installed
#' library(penalizedG)
#'
#' # Run penalized G-estimation for a sequence of tuning parameters (lambda)
#' out <- penalizedG(data = data.s, wc.str = "exchangeable", id.var="id", response.var="y",
#'                    treat.var="a", tf.model=tf.model, treat.model = treat.model,
#'                    lambda.seq = lambda.seq, maxitr = 100, penalty = "SCAD")
#' out$selected.EMs
#'
#' data <- out$data # have an additional column containing the propensity scores
#'
#' # confidence interval for OSIPeG (using the decorrelated score approach)
#' library(flare)
#' CI.osipeg <- osipeg_confint(data = data, wc.str = "exchangeable", id.var = "id", response.var = "y",
#'                             treat.var = "a", tf.model = tf.model, estimate = out$estimate,
#'                             alpha.hat = out$alpha.hat, sigma2.hat = out$sigma2.hat,
#'                             lambda = out$lambda.optimal, test.size = 0.05,
#'                            tune.par.seq = sort(seq(0,0.2,(0.2-0)/100), decreasing = TRUE),
#'                            continuous.covs = rep(TRUE, length(all.vars(tf.model))))
#' CI.osipeg$CI.one.step.full
#' CI.osipeg$CI.one.step.LASSO
#' CI.osipeg$CI.one.step.dantzig
#'
#' # Naive CI based on sandwich variance
#' p <- length(model.matrix(tf.model, data=data)[1,])
#' psi.hat <- out$estimate[1:p]
#' effect.modification <- psi.hat[2:p]
#' M <- c(TRUE, abs(effect.modification) > 0.001)
#'
#' CI.naive <- cbind(psi.hat[M] - qnorm(0.975)*sqrt(diag(out$asymp.var.psi[M,M])),
#'                   psi.hat[M] + qnorm(0.975)*sqrt(diag(out$asymp.var.psi[M,M])))
#' CI.naive

osipeg_confint <- function(data, wc.str, id.var, treat.var, response.var, tf.model, estimate, alpha.hat,
                           sigma2.hat, lambda, test.size, tune.par.seq, continuous.covs){

  names(data)[names(data)==id.var] <- "id"
  names(data)[names(data)==treat.var] <- "a"
  names(data)[names(data)==response.var] <- "y"

  for(k in all.vars(tf.model)[continuous.covs]){
    data[,k] <- (data[,k] - mean(data[,k]))/sd(data[,k])
  }

  n <- length(split(data, data$id))
  dat <- split(data, data$id)
  l.mat <- data.frame(id=data$id,model.matrix(tf.model, data)) # cov+treat history
  l.mat.split <- split(l.mat, l.mat$id)
  ni <- lapply(1:n, function(i) length(dat[[i]]$id))
  l <- lapply(1:n, function(i) as.matrix(l.mat.split[[i]][,-1]))
  a <- lapply(1:n, function(i) as.matrix(dat[[i]]$a))
  y <- lapply(1:n, function(i) as.matrix(dat[[i]]$y))
  E.a <- lapply(1:n, function(i) as.matrix(dat[[i]]$E.a))
  p <- dim(l.mat[,-1])[2]

  ## Find which coeffs are selected
  coef.sel <- which(c(TRUE, abs(estimate[2:p]) > 0.001))

  En.mat <- diag(c(0,q_scad(lambda,abs(estimate[2:p]))/(10^(-6)+abs(estimate[2:p])),
                   rep(0,p)))
  e <- lapply(1:n, function(i)
    y[[i]]-(c(a[[i]])*l[[i]])%*%estimate[1:p]-l[[i]]%*%estimate[(p+1):(2*p)])
  V <- lapply(1:n, function(i) sigma2.hat*corMatF(alpha = alpha.hat, ni = ni[[i]], corstr = wc.str))
  V.inv <- lapply(1:n, function(i) solve(V[[i]]))

  sum3.n <- Reduce("+", lapply(1:n, function(i)
    t(c(a[[i]]-E.a[[i]])*l[[i]])%*%V.inv[[i]]%*%(c(a[[i]])*l[[i]])))
  ss4 <- Reduce("+", lapply(1:n, function(i)
    t(c(a[[i]]-E.a[[i]])*l[[i]])%*%V.inv[[i]]%*%e[[i]]%*%t(e[[i]])%*%
      V.inv[[i]]%*%(c(a[[i]]-E.a[[i]])*l[[i]])))

  I.psi.hat <- ss4
  B.hat <- sum3.n+n*En.mat[1:p,1:p]
  B.hat.inv <- solve(B.hat)
  var <- diag(B.hat.inv%*%I.psi.hat%*%B.hat.inv)

  K <- p
  psi.hat.sparse <- estimate[1:K]
  S.psi.i <- lapply(1:n, function(i)
    t(c(a[[i]]-E.a[[i]])*l[[i]])%*%V.inv[[i]]%*%e[[i]])
  I.psi.i <- lapply(1:n, function(i) S.psi.i[[i]]%*%t(S.psi.i[[i]]))
  S.psi <- (1/n)*Reduce("+", S.psi.i)
  I.psi <- (1/n)*Reduce("+", I.psi.i)
  ## S.psi.mat is different from S.psi; S.psi is p*1 but S.psi.mat is n*p
  S.psi.mat <- sapply(1:n, function(i)
    t(c(a[[i]]-E.a[[i]])*l[[i]])%*%V.inv[[i]]%*%e[[i]], simplify = TRUE)
  S.psi.mat <- t(S.psi.mat)

  tune.par.opt <- tune.par.opt.dantzig <- NULL
  psi.hat.decorr.OS <- psi.hat.decorr.OS.LASSO <- psi.hat.decorr.OS.dantzig <- NULL
  CI1 <- CI2 <- CI3 <- NULL
                      
  for(k in coef.sel){
    I.psi.partial <- c(I.psi[k,k] - I.psi[k,-k]%*%solve(I.psi[-k,-k])%*%I.psi[-k,k])

    ## using full weight vector
    W <- I.psi[k,-k]%*%solve(I.psi[-k,-k])
    weight.all <- rep(1, K)
    weight.all[-k] <- -W
    S.psi.decorr <- c(S.psi[k] - sum(W*S.psi[-k]))
    psi.hat.decorr <- psi.hat.sparse[k] - S.psi.decorr/I.psi.partial
    sigma.hat <- c(t(as.matrix(weight.all))%*%I.psi%*%as.matrix(weight.all))
    CI1 <- rbind(CI1, psi.hat.decorr + c(-qnorm(1-test.size/2), qnorm(1-test.size/2))*
                   sqrt(sigma.hat)/(sqrt(n)*I.psi.partial))
    psi.hat.decorr.OS <- c(psi.hat.decorr.OS, psi.hat.decorr)

    ## W estimated with LASSO
    out <- optimalW_lasso(y=S.psi.mat[,k], X=S.psi.mat[,-k], tune.par.seq=tune.par.seq)
    tune.par.opt <- c(tune.par.opt, out$tune.par.opt)
    W <- out$W.hat
    weight.all <- rep(1, K)
    weight.all[-k] <- -W
    S.psi.decorr <- c(S.psi[k] - sum(W*S.psi[-k]))
    psi.hat.decorr <- psi.hat.sparse[k] - S.psi.decorr/I.psi.partial
    sigma.hat <- c(t(as.matrix(weight.all))%*%I.psi%*%as.matrix(weight.all))
    CI2 <- rbind(CI2, psi.hat.decorr + c(-qnorm(1-test.size/2), qnorm(1-test.size/2))*
                   sqrt(sigma.hat)/(sqrt(n)*I.psi.partial))
    psi.hat.decorr.OS.LASSO <- c(psi.hat.decorr.OS.LASSO, psi.hat.decorr)

    ## W estimated with dantzig
    out <- optimalW_dantzig(y=S.psi.mat[,k], X=S.psi.mat[,-k], tune.par.seq=tune.par.seq)
    tune.par.opt.dantzig <- c(tune.par.opt.dantzig, out$tune.par.opt)
    W <- out$W.hat
    weight.all <- rep(1, K)
    weight.all[-k] <- -W
    S.psi.decorr <- c(S.psi[k] - sum(W*S.psi[-k]))
    psi.hat.decorr <- psi.hat.sparse[k] - S.psi.decorr/I.psi.partial
    sigma.hat <- c(t(as.matrix(weight.all))%*%I.psi%*%as.matrix(weight.all))
    CI3 <- rbind(CI3, psi.hat.decorr + c(-qnorm(1-test.size/2), qnorm(1-test.size/2))*
                   sqrt(sigma.hat)/(sqrt(n)*I.psi.partial))
    psi.hat.decorr.OS.dantzig <- c(psi.hat.decorr.OS.dantzig, psi.hat.decorr)
    #print(k)
  }
  CI1 <- as.data.frame(CI1)
  CI2 <- as.data.frame(CI2)
  CI3 <- as.data.frame(CI3)
  names(CI1) <- names(CI2) <- names(CI3) <- c("lower", "upper")

  M <- c(TRUE, abs(estimate[2:p]) > 0.001)
  selected.EMs <- colnames(model.matrix(tf.model, data))[M][-1]
  row.names(CI1) <- row.names(CI2) <- row.names(CI3) <- c("a", paste("a*", selected.EMs, sep=""))

  return(list(psi.hat.one.step.full = psi.hat.decorr.OS, psi.hat.one.step.LASSO =  psi.hat.decorr.OS.LASSO,
              psi.hat.one.step.dantzig = psi.hat.decorr.OS.dantzig,
              CI.one.step.full = CI1, CI.one.step.LASSO = CI2, CI.one.step.dantzig = CI3))
}

