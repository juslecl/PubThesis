library(readxl)
library(stats)
library(optimx)
library(BiocParallel)
library(MASS)
library(VGAM)
library(mixtools)

# Upload data and visualize it
data(Boston)
set.seed(1234)
model <- lm(medv ~ ., data = Boston)
summary(model)
plot(Boston$lstat, Boston$medv)
abline(lm(medv ~ lstat, data = Boston), col = "red")

#-------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------

# Residuals handling and density estimation.
# Normalization
resids <- residuals(model)
resids_scaled <- resids/sd(resids)

# Density estimation with mixtools::normalmixEM
set.seed(1234)
est_param <- normalmixEM(resids_scaled, k=3)

# Estimated mixed normal function
mixed_dens <- function(x){
  return(est_param$lambda[1]*dnorm(x, mean=est_param$mu[1], sd=est_param$sigma[1])+est_param$lambda[2]*dnorm(x, mean=est_param$mu[2], sd=est_param$sigma[2])+est_param$lambda[3]*dnorm(x, mean=est_param$mu[3], sd=est_param$sigma[3]))
}

par(mfrow=c(1,2))
plot(density(resids), main="Residuals' density",col="blue")
plot(density(resids_scaled), main="Scaled residuals' density",col="blue")
curve(mixed_dens,from=-4,to=6,add=T,col="red")

#-------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------

bootids <- function(data,iter=100){
  ids <- matrix(NA,ncol=iter,nrow=floor(nrow(data)*.6))
  for (i in seq_len(iter)){
    ids[,i] <- sample(seq_len(nrow(data)), floor(nrow(data)*.6))
  }
  return(ids)
}

ols <- function(data) {
  return(solve(t(data$X)%*%data$X)%*%t(data$X)%*%data$Y)
}

neg_loglik <- function(param, data) { #param is a vector of length(beta)+1
  t <- (data$Y-data$X%*%param[1:(length(param)-1)])/param[length(param)]
  s <- param[length(param)]
  log_f <- -log(s) + log(mixed_dens(t))
  return(-sum(log_f))
}

compute_mle_beta <- function(data) {
  ols_est <- ols(data)
  s_start <- sqrt(var(data$Y - data$X%*%ols_est))
  neg_loglik_wrapped <- function(param) {
    neg_loglik(param, data = data)
  }
  mle <- optimx(par = c(as.numeric(ols_est), s_start), fn = neg_loglik_wrapped, method = "BFGS")
  resvec <- mle[, 1:(ncol(data$X)+1)]
  return(mle=t(resvec))
}

mleols <- function(data){
  param_cfc <- BatchtoolsParam(workers = 32)
  id <- bootids(data)
  res <- bptry({bplapply(seq_len(ncol(id)), function(boots){
    rdata_train <- list(X=as.matrix(cbind(rep(1,nrow(data[id[,boots],])),data[id[,boots],-ncol(data)])),Y=as.matrix(data[id[,boots],ncol(data)]))
    mle <- compute_mle_beta(rdata_train)
    olse <- c(ols(rdata_train),NA)
    return(cbind.data.frame(mle=mle,ols=olse))
  },BPPARAM = param_cfc)})
  return(list(res,id))
}

pred <- function(mleols_obj){
  ols_tab <- matrix(unlist(lapply(mleols_obj[[1]], function(i) return(i[-nrow(i),2]))), ncol=length(mleols_obj[[1]]), byrow=F)
  mle_tab <- matrix(unlist(lapply(mleols_obj[[1]], function(i) return(i[-nrow(i),1]))), ncol=length(mleols_obj[[1]]), byrow=F)
  return(list(OLSE=ols_tab,MLE=mle_tab))
}

test <- function(mleols_obj,data){
  olse <- pred(mleols_obj)$OLSE
  mle <- pred(mleols_obj)$MLE
  param_test <- BatchtoolsParam(workers = 32)
  param_res <- BatchtoolsParam(workers = 32)
  id <- mleols_obj[[2]]
  res <- bptry({bplapply(seq_len(ncol(id)), function(boots){
    rdata_test <- list(X=as.matrix(cbind(rep(1,nrow(data[-id[,boots],])), data[-id[,boots],-ncol(data)])),Y=as.matrix(data[-id[,boots],ncol(data)]))
    fittedy_mle <- rdata_test$X%*%mle[,boots]
    fittedy_olse <- rdata_test$X%*%olse[,boots]
    res_mle <- fittedy_mle - rdata_test$Y
    res_olse <- fittedy_olse - rdata_test$Y
    return(cbind.data.frame(fitted_mle=fittedy_mle,fitted_olse=fittedy_olse,resid_mle=res_mle,resid_olse=res_olse))
  },BPPARAM=param_test)})
  mean_res <- bptry({bplapply(res, function(data){
    return(apply(data[,3:4],2, function(i) mean(i^2)))
  }, BPPARAM=param_test)})
  mean_res <- matrix(unlist(mean_res),ncol=2,byrow=T)
  colnames(mean_res) <- c("MLE","OLSE")
  return(list(res,mean_res))
}

set.seed(1234)
betatry1 <- mleols(Boston)
beta1 <- pred(betatry1)
test1 <- test(betatry1,Boston)
apply(test1[[2]],2,mean)
