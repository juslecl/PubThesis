library(readxl)
library(stats)
library(optimx)
RE_data_set <- read_excel("Desktop/Real estate valuation data set.xlsx")

# Reformat and split.
bootids <- function(data,iter=100){
  ids <- matrix(NA,ncol=iter,nrow=ncol(data)*0.6)
  for (i in seq_len(iter)){
    ids[,i] <- sample(seq_len(ncol(data)), ncol(sce)*0.6)
  }
  return(ids)
}

set.seed(123)
id_train <- sample(seq_len(nrow(RE_data_set)),nrow(RE_data_set)*.5)
RE_train <- list(Y=RE_data_set$`Y house price of unit area`,X=as.matrix(RE_data_set[2:7]))
RE
ols <- function(data) {
  return(solve(t(data$X[id_train,])%*%data$X[id_train,])%*%t(data$X[id_train,])%*%data$Y[id_train])
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
  res <- optimx(par = c(as.numeric(ols_est), s_start), fn = neg_loglik_wrapped, method = "BFGS")
  resvec <- res[, 1:(ncol(data$X)+1)]
  return(mle=t(resvec))
}

#-------------------------------------------------------------
#-------------------------------------------------------------
mleols <- function(data, g, type){
  mle <- compute_mle_beta(data, g, type)
  olse <- c(ols(data),NA)
  return(cbind.data.frame(mle=mle,ols=olse))
}

beta_OLS <- ols(RE_data)
beta_est1 <- mleols(RE_data, g=3, type=1)
