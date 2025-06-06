set.seed(123)
library(stats)
library(optimx)
library(BiocParallel)
#0) Create functions to draw errors followign one of the distributions at random.
rmod1 <- function(n, g, s) {
  library(normalp)
  # Constants
  c_gamma <- (gamma(3 / g) / gamma(1 / g))^(g/ 2)
  # Draw from exp(-|x|^gamma) then rescale
  x <- rnormp(n, mu = 0, sigmap = 1, p = g)
  scale <- s*c_gamma^(-1 / g)
  return(x * scale)
}

rmod2 <- function(n, g, s) {
  c_gamma <- (gamma(g + 2) / gamma(g))^(1/2)
  # Draw from Gamma(gamma, rate = c_gamma)
  t <- rgamma(n, shape = g, rate = c_gamma)
  # Symmetrize
  signs <- s*sample(c(-1, 1), n, replace = TRUE)
  return(signs * t)
}

rmod3 <- function(n, g, s) {
  c_gamma <- gamma(2 / g + 1)^(g/ 2)
  # Scale parameter for Weibull: b = c_gamma^{-1/gamma}
  scale <- c_gamma^(-1 / g)
  t <- rweibull(n, shape = g, scale = scale)
  signs <- s*sample(c(-1, 1), n, replace = TRUE)
  return(signs * t)
}

#-------------------------------------------------------------
#-------------------------------------------------------------

#1) Create a dataset with X following a normal but different patterns for the error.
datagen <- function(n, param, rcdf, g, meanX, sdX,  M=100) {
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
  
  if (type == 1) {
    # Type 1: f(t) = d * exp(-c * |t|^gamma)
    c_gamma <- (gamma(3 / g) / gamma(1 / g))^(g/ 2)
    d_gamma <- (g/ 2) * (gamma(3 / g) / gamma(1 / g)^3)^(1/2)
    log_f <- -log(s) + log(d_gamma) - c_gamma * abs(t)^g
    
  } else if (type == 2) {
    # Type 2: f(t) = d * |t|^{gamma - 1} * exp(-c * |t|), gamma ≥ 1
    c_gamma <- (gamma(g + 2) / gamma(g))^(1/2)
    d_gamma <- (1 / (2 * gamma(g))) * (gamma(g+ 2) / gamma(g))^(g/ 2)
    log_f <- -log(s) + log(d_gamma) + (g - 1) * log(abs(t)) - c_gamma * abs(t)
    
  } else if (type == 3) {
    # Type 3: f(t) = d * |t|^{gamma - 1} * exp(-c * |t|^gamma), gamma ≥ 1
    c_gamma <- gamma(2 / g + 1)^(g / 2)
    d_gamma <- (g / 2) * gamma(2 / g + 1)^(g / 2)
    log_f <- -log(s) + log(d_gamma) + (g - 1) * log(abs(t)) - c_gamma * abs(t)^g
  }
  
  return(-sum(log_f))
}

# Function to compute MLE of beta. 
compute_mle_beta <- function(data, g, type = 1) {
  neg_loglik_wrapped <- function(param) {
    neg_loglik(param, data = data, g = g, type = type)
  }
  
  res <- optimx(par = c(rep(0, ncol(data$X)),1), fn = neg_loglik_wrapped, method = "BFGS")
  
  resvec <- res[, 1:(ncol(data$X)+1)]  # fixed res -> res
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
  return(res)
}

#-------------------------------------------------------------
#-------------------------------------------------------------
metrics <- function(mleols_res){
  param <- BatchtoolsParam(workers=32)
  M <- length(mleols_res)
  norms <- bptry({bplapply(seq_len(M), function(nb){
    dati <- mleols_res[[nb]][-nrow(mleols_res[[1]]),]
    norm2ols <- norm(dati[,2]-dati[,3],type="2")
    norm2mle <- norm(dati[,1]-dati[,3],type="2")
    return(list(norm2ols,norm2mle))
  },BPPARAM=param)})
  empirical_means <- vector("list", 2)
  empirical_vars  <- vector("list", 2)
  for (j in 1:2) {
    mat_j <- sapply(mleols_res, function(m) m[-nrow(mleols_res[[1]]), j])
    empirical_means[[j]] <- rowMeans(mat_j)
    empirical_vars[[j]]  <- apply(mat_j, 1, var)
  }
  names(empirical_means) <- names(empirical_vars)  <- paste0("vector_", c("mle","olse"))
  return(list(norms,empirical_means,empirical_vars))
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

data11 <- datagen(n, beta, rmod1 , 3, meanX, sdX, M=100)
data12 <- datagen(n, beta, rmod2 , 3, meanX, sdX, M=100)
data13 <- datagen(n, beta, rmod3 , 3, meanX, sdX, M=100)

mleols11 <- mleols(data11, 3, 1)
mleols12 <- mleols(data12, 3, 2)
mleols13 <- mleols(data13, 3, 3)

metrics11 <- metrics(mleols11)
metrics12 <- metrics(mleols12)
metrics13 <- metrics(mleols13)
