set.seed(123)
library(vcvComp)
library(stats)
library(optimx)
library(BiocParallel)

#0) Create functions to draw errors followign one of the distributions at random.
rmod1 <- function(n, g, s) {
  library(gnorm)
  c_gamma <- (gamma(3 / g) / gamma(1 / g))^(g / 2)
  x <- rgnorm(n, mu = 0, alpha = 1, beta = g)
  scale <- s * c_gamma^(-1 / g)
  return(x * scale)
}

rmod2 <- function(n, g, s) {
  c_gamma <- (gamma(g + 2) / gamma(g))^(1/2)
  t <- rgamma(n, shape = g, scale = c_gamma^(-1))
  signs <- s*sample(c(-1, 1), n, replace = TRUE)
  return(signs * t)
}

rmod3 <- function(n, g, s) {
  c_gamma <- gamma(2 / g + 1)^(g/ 2)
  scale <- c_gamma^(-1 / g)
  t <- rweibull(n, shape = g, scale = scale)
  signs <- s*sample(c(-1, 1), n, replace = TRUE)
  return(signs * t)
} 

#-------------------------------------------------------------
#-------------------------------------------------------------

#1) Create a dataset with X following a normal but different patterns for the error.
datagen <- function(n, param, rcdf, g, meanX, sdX,  M=1000) {
  beta <- param[-length(param)]
  s <- param[length(param)]
  d <- length(meanX)
  datasets <- vector("list", M) 
  for (m in 1:M) {  
    x <- matrix(data = NA, nrow = n, ncol = d)  
    for (j in 1:d)  x[,j] <- rnorm(n,mean=meanX[j],sd=sdX[j])
    eps <- rcdf(n,g,s)  
    y <- x%*%beta + eps  
    datasets[[m]] <- list(Y = y, X = x, P = param, E = eps)  
  }
  return(datasets)  
}

EXXT <- function(datasets){
  d <- ncol(datasets[[1]]$X)
  n <- nrow(datasets[[1]]$X)
  EXXT <- matrix(0, nrow = d, ncol = d)
  for (data in datasets) {
    X <- data$X
    EXXT <- EXXT + t(X) %*% X
  }
  EXXT <- EXXT / (length(datasets) * n)
  return(EXXT)
}

#-------------------------------------------------------------
#-------------------------------------------------------------

#2) Function to build OLS. => retourne une dx1 matrice.
ols <- function(data) {
  return(solve(t(data$X)%*%data$X)%*%t(data$X)%*%data$Y)
}

#-------------------------------------------------------------
#-------------------------------------------------------------

# Generic negative log-likelihood function
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

# Function to compute MLE of beta. 
compute_mle_beta <- function(data, g, type = 1) {
  ols_est <- ols(data)
  print(ols_est)
  s_start <- sqrt(var(data$Y - data$X%*%ols_est))
  print(s_start)# or MAD-based
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
  param <- BatchtoolsParam(workers=32)
  M <- length(data)
  res <- bptry({bplapply(seq_len(M), function(nb){
    mle <- compute_mle_beta(data[[nb]], g, type)
    olse <- c(ols(data[[nb]]),NA)
    return(cbind.data.frame(mle=mle,ols=olse,truth=data[[nb]]$P))
  },BPPARAM=param)})
  ordremle <- matrix(unlist(lapply(res,function(i) return(i[-nrow(i),1]))),ncol=5,byrow=T)
  ordreols <- matrix(unlist(lapply(res,function(i) return(i[-nrow(i),2]))),ncol=5,byrow=T)
  return(list(res,ordremle,ordreols))
}

#-------------------------------------------------------------
#-------------------------------------------------------------
metrics <- function(mleols_res){
  param <- BatchtoolsParam(workers=32)
  mleols_res <- mleols_res[[1]]
  M <- length(mleols_res)
  norms <- bptry({bplapply(seq_len(M), function(nb){
    dati <- mleols_res[[nb]][-nrow(mleols_res[[1]]),]
    norm2ols <- norm(dati[,2]-dati[,3],type="2")
    norm2mle <- norm(dati[,1]-dati[,3],type="2")
    return(list(norm2ols,norm2mle))
  },BPPARAM=param)})
  norms <- matrix(unlist(norms),nrow=M,ncol=2,byrow=T)
  colnames(norms) <- c("OLSE","MLE")
  empirical_means <- vector("list", 2)
  empirical_vars  <- vector("list", 2)
  for (j in 1:2) {
    mat_j <- sapply(mleols_res, function(m) m[-nrow(mleols_res[[1]]), j])
    empirical_means[[j]] <- rowMeans(mat_j)
    empirical_vars[[j]]  <- cov(t(mat_j))
  }
  empirical_means <- matrix(unlist(empirical_means),ncol=2,byrow=F)
  colnames(empirical_means) <- c("MLE","OLSE")
  return(list(NORMS=norms,MEANS=empirical_means,VARS=empirical_vars,MSE=apply(norms,2,function(i) mean(i^2))))
}

#-------------------------------------------------------------
#-------------------------------------------------------------

# SIMULATION 1
set.seed(123)
d <- 5
n <- 500
beta <- c(sample(seq(from=-3,to=3,by=.01),d),sample(seq(from=1,to=3,by=.01),1))
meanX <- sample(seq(from=-5,to=5,by=.01),d)
sdX <- sample(seq(from=0,to=2,by=.1),d)

for (i in 1:3) {
  rmod <- get(paste0("rmod", i))
  assign(paste0("data1", i), datagen(n, beta, rmod, 7, meanX, sdX))
}

for (i in 1:3) {
  data_obj <- get(paste0("data1", i))
  assign(paste0("mleols1", i), mleols(data_obj, 7, i))
}

for (i in 1:3) {
  mleols_obj <- get(paste0("mleols1", i))
  assign(paste0("metrics1", i), metrics(mleols_obj))
}

par(mfrow = c(1, 3))
for (i in 1:3) {
  metrics_obj <- get(paste0("metrics1", i))
  boxplot(metrics_obj[[1]], ylab = "Euclidean norm", main = paste("Model", i))
}

par(mfrow=c(1,3))
for (i in 1:3){
plot(density(get(paste("mleols1", i, sep = ""))[[2]][, 2]),ylab="Density", xlab=expression(hat(beta)),col="red",main="")
lines(density(get(paste("mleols1", i, sep = ""))[[3]][, 2]),ylab="Density", xlab=expression(hat(beta)),col="blue")
legend("topleft", legend=c("MLE","OLSE"),col=c("red","blue"),pch = "_")
}

# SIMULATION 1
set.seed(123)
d <- 5
n <- 500
beta <- c(sample(seq(from=-3,to=3,by=.01),d),sample(seq(from=1,to=3,by=.01),1))
meanX <- sample(seq(from=-5,to=5,by=.01),d)
sdX <- sample(seq(from=0,to=2,by=.1),d)

for (i in 1:3) {
  rmod <- get(paste0("rmod", i))
  assign(paste0("data1", i), datagen(n, beta, rmod, 5, meanX, sdX))
}

for (i in 1:3) {
  data_obj <- get(paste0("data1", i))
  assign(paste0("mleols1", i), mleols(data_obj, 5, i))
}

for (i in 1:3) {
  mleols_obj <- get(paste0("mleols1", i))
  assign(paste0("metrics1", i), metrics(mleols_obj))
}

par(mfrow = c(1, 3))
for (i in 1:3) {
  metrics_obj <- get(paste0("metrics1", i))
  boxplot(metrics_obj[[1]], ylab = "Euclidean norm", main = paste("Model", i))
}

par(mfrow=c(1,3))
for (i in 1:3){
  plot(density(get(paste("mleols1", i, sep = ""))[[2]][, 2]),ylab="Density", xlab=expression(hat(beta)),col="red",main="")
  lines(density(get(paste("mleols1", i, sep = ""))[[3]][, 2]),ylab="Density", xlab=expression(hat(beta)),col="blue")
  legend("topleft", legend=c("MLE","OLSE"),col=c("red","blue"),pch = "_")
}


# SIMULATION 1
set.seed(123)
d <- 5
n <- 1000
beta <- c(sample(seq(from=-3,to=3,by=.01),d),sample(seq(from=1,to=3,by=.01),1))
meanX <- sample(seq(from=-5,to=5,by=.01),d)
sdX <- sample(seq(from=0,to=2,by=.1),d)

for (i in 1:3) {
  rmod <- get(paste0("rmod", i))
  assign(paste0("data1", i), datagen(n, beta, rmod, 5, meanX, sdX))
}

for (i in 1:3) {
  data_obj <- get(paste0("data1", i))
  assign(paste0("mleols1", i), mleols(data_obj, 5, i))
}

for (i in 1:3) {
  mleols_obj <- get(paste0("mleols1", i))
  assign(paste0("metrics1", i), metrics(mleols_obj))
}

par(mfrow = c(1, 3))
for (i in 1:3) {
  metrics_obj <- get(paste0("metrics1", i))
  boxplot(metrics_obj[[1]], ylab = "Euclidean norm", main = paste("Model", i))
}

par(mfrow=c(1,3))
for (i in 1:3){
  plot(density(get(paste("mleols1", i, sep = ""))[[2]][, 2]),ylab="Density", xlab=expression(hat(beta)),col="red",main="")
  lines(density(get(paste("mleols1", i, sep = ""))[[3]][, 2]),ylab="Density", xlab=expression(hat(beta)),col="blue")
  legend("topleft", legend=c("MLE","OLSE"),col=c("red","blue"),pch = "_")
}

