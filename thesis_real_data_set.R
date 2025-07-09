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

neg_loglik <- function(param, data, g, type = 1) { #param is a vector of length(beta)+1
  t <- (data$Y-data$X%*%param[1:(length(param)-1)])/param[length(param)]
  s <- param[length(param)]
  print(s)
  
  if (type == 1) {
    # Type 1: f(t) = d/s * exp(-c * |t|^gamma)
    c_gamma <- (gamma(3 / g) / gamma(1 / g))^(g/ 2)
    d_gamma <- (g/ 2) * (gamma(3 / g) / gamma(1 / g)^3)^(1/2)
    log_f <- -log(s) + log(d_gamma) - c_gamma * abs(t)^g
    
  } else if (type == 2) {
    # Type 2: f(t) = d/s * |t|^{gamma - 1} * exp(-c * |t|), gamma ≥ 1
    c_gamma <- (gamma(g + 2) / gamma(g))^(1/2)
    d_gamma <- (1 / (2 * gamma(g))) * (gamma(g+ 2) / gamma(g))^(g/ 2)
    log_f <- -log(s) + log(d_gamma) + (g - 1) * log(abs(t)) - c_gamma * abs(t)
    
  } else if (type == 3) {
    # Type 3: f(t) = d/s * |t|^{gamma - 1} * exp(-c * |t|^gamma), gamma ≥ 1
    c_gamma <- gamma(2 / g + 1)^(g / 2)
    d_gamma <- (g / 2) * gamma(2 / g + 1)^(g / 2)
    log_f <- -log(s) + log(d_gamma) + (g - 1) * log(abs(t)) - c_gamma * abs(t)^g
  }
  return(-sum(log_f))
}

compute_mle_beta <- function(data, g, type = 1) {
  ols_est <- ols(data)
  s_start <- sqrt(var(data$Y - data$X%*%ols_est))
  neg_loglik_wrapped <- function(param) {
    neg_loglik(param, data = data, g = g, type = type)
  }
  mle <- optimx(par = c(as.numeric(ols_est), s_start), fn = neg_loglik_wrapped, method = "BFGS")
  resvec <- mle[, 1:(ncol(data$X)+1)]
  return(mle=t(resvec))
}

mleols <- function(data, g, type){
  param_cfc <- BatchtoolsParam(workers = 32)
  id <- bootids(data)
  res <- bptry({bplapply(seq_len(ncol(id)), function(boots){
    rdata_train <- list(X=as.matrix(data[id[,boots],-ncol(data)]),Y=as.matrix(data[id[,boots],ncol(data)]))
    mle <- compute_mle_beta(rdata_train, g, type)
    olse <- c(ols(rdata_train),NA)
    return(cbind.data.frame(mle=mle,ols=olse))
  },BPPARAM = param_cfc)})
  return(list(res,id))
}

pred <- function(mleols_obj){
  ols_tab <- matrix(unlist(lapply(mleols_obj[[1]], function(i) return(i[-nrow(i),2]))), ncol=length(mleols_obj[[1]]), byrow=F)
  mle_tab <- matrix(unlist(lapply(mleols_obj[[1]], function(i) return(i[-nrow(i),1]))), ncol=length(mleols_obj[[1]]), byrow=F)
  olseest <- apply(ols_tab,1,mean)
  mleest <- apply(mle_tab,1,mean)
  return(list(OLSE=ols_tab,MLE=mle_tab,OLSE_est=olseest, MLE_est=mleest))
}

test <- function(mleols_obj,data){
  olse <- pred(mleols_obj)$OLSE
  mle <- pred(mleols_obj)$MLE
  param_test <- BatchtoolsParam(workers = 32)
  param_res <- BatchtoolsParam(workers = 32)
  id <- mleols_obj[[2]]
  res <- bptry({bplapply(seq_len(ncol(id)), function(boots){
    rdata_test <- list(X=as.matrix(data[-id[,boots],-ncol(data)]),Y=as.matrix(data[-id[,boots],ncol(data)]))
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

# add an intercept.


set.seed(1234)
betatry1 <- mleols(hprd, g=2, type=1)
betatry2 <- mleols(hprd, g=2, type=2)
betatry3 <- mleols(hprd, g=2, type=3)

beta1 <- pred(betatry1)
beta2 <- pred(betatry2)
beta3 <- pred(betatry3)

test1 <- test(betatry1,hprd)
test2 <- test(betatry2,hprd)
test3 <- test(betatry3,hprd)