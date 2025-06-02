set.seed(123)
library(stats)
library(optimx)
#0) Create functions to draw errors followign one of the distributions at random.
rmod1 <- function(n, g) {
  library(normalp)
  # Constants
  c_gamma <- (gamma(3 / g) / gamma(1 / g))^(g/ 2)
  # Draw from exp(-|x|^gamma) then rescale
  x <- rnormp(n, mu = 0, sigmap = 1, p = g)
  scale <- c_gamma^(-1 / g)
  return(x * scale)
}

rmod2 <- function(n, g) {
  c_gamma <- (gamma(g + 2) / gamma(g))^(1/2)
  # Draw from Gamma(gamma, rate = c_gamma)
  t <- rgamma(n, shape = g, rate = c_gamma)
  # Symmetrize
  signs <- sample(c(-1, 1), n, replace = TRUE)
  return(signs * t)
}

rmod3 <- function(n, g) {
  c_gamma <- gamma(2 / g + 1)^(g/ 2)
  # Scale parameter for Weibull: b = c_gamma^{-1/gamma}
  scale <- c_gamma^(-1 / g)
  t <- rweibull(n, shape = g, scale = scale)
  signs <- sample(c(-1, 1), n, replace = TRUE)
  return(signs * t)
}

#-------------------------------------------------------------
#-------------------------------------------------------------


#1) Create a dataset with X following a normal but different patterns for the error.
datagen <- function(n, beta, rcdf, g, meanX, sdX,  M=100) {
  d <- length(meanX)
  datasets <- vector("list", M) 
  for (m in 1:M) {  
    x <- matrix(data = NA, nrow = n, ncol = d)  
    for (j in 1:d)  x[,j] <- rnorm(n,mean=meanX[j],sd=sdX[j])
    eps <- rcdf(n,g)  
    y <- x%*%beta + eps  
    datasets[[m]] <- list(Y = y, X = x, B = beta, E = eps)  
  }
  return(datasets)  
}

datatry <- datagen(500, 1:5, rmod1 , 3, rep(0,5), rep(1,5), M=100)
  
  
#-------------------------------------------------------------
#-------------------------------------------------------------



#2) Function to build OLS. => retourne une dx1 matrice.
ols <- function(data) {
  return(solve(t(data$X)%*%data$X)%*%t(data$X)%*%data$Y)
}

olstry <- ols(datatry[[1]])

#-------------------------------------------------------------
#-------------------------------------------------------------

# Generic negative log-likelihood function
neg_loglik <- function(beta, data, g, type = 1) {
  t <- data$Y-data$X%*%beta
  
  if (type == 1) {
    # Type 1: f(t) = d * exp(-c * |t|^gamma)
    c_gamma <- (gamma(3 / g) / gamma(1 / g))^(g/ 2)
    d_gamma <- (g/ 2) * (gamma(3 / g) / gamma(1 / g)^3)^(1/2)
    log_f <- log(d_gamma) - c_gamma * abs(t)^g
    
  } else if (type == 2) {
    # Type 2: f(t) = d * |t|^{gamma - 1} * exp(-c * |t|), gamma ≥ 1
    c_gamma <- (gamma(g + 2) / gamma(g))^(1/2)
    d_gamma <- (1 / (2 * gamma(g))) * (gamma(g+ 2) / gamma(g))^(g/ 2)
    log_f <- log(d_gamma) + (g - 1) * log(abs(t)) - c_gamma * abs(t)
    
  } else if (type == 3) {
    # Type 3: f(t) = d * |t|^{gamma - 1} * exp(-c * |t|^gamma), gamma ≥ 1
    c_gamma <- gamma(2 / g + 1)^(g / 2)
    d_gamma <- (g / 2) * gamma(2 / g + 1)^(g / 2)
    log_f <- log(d_gamma) + (g - 1) * log(abs(t)) - c_gamma * abs(t)^g
  }
  
  return(-sum(log_f))
}

# Function to compute MLE of beta. A PARTIR DE LA CA MARCHE PLUS...
compute_mle_beta <- function(data, g, type = 1) {
  res <- optimx(par=c(rep(0,length(data$B))), fn=neg_loglik,
                data = data, g = g, type = type, method="BFGS")
  resvec <- fit[ ,1:length(data$B)]
  return(t(resvec))
}

mletry <- compute_mle_beta(datatry[[3]],g=3,type=1)

## Example you should anzi try on monday...................
mllaplace <- function(data,be) {
  fit <- optimx(par=c(rep(0,length(data$B))), fn=neg_lllaplace, X=data$X, Y=data$Y, b=be, method="BFGS")
  resvec <- fit[ ,1:length(data$B)]
  return(t(resvec))
}

## Making the wole thing work...
dati <- datagen()

#-------------------------------------------------------------
#-------------------------------------------------------------

#2) Function to build OLS. => retourne une dx1 matrice.
ols <- function(data) {
  return(solve(t(data$X)%*%data$X)%*%t(data$X)%*%data$Y)
}

#-------------------------------------------------------------
#-------------------------------------------------------------

# Generic negative log-likelihood function
neg_loglik <- function(data, g, type = 1) {
  t <- data$Y-data$X%*%data$beta
  
  if (type == 1) {
    # Type 1: f(t) = d * exp(-c * |t|^gamma)
    c_gamma <- (gamma(3 / g) / gamma(1 / g))^(g/ 2)
    d_gamma <- (g/ 2) * (gamma(3 / g) / gamma(1 / g)^3)^(1/2)
    log_f <- log(d_gamma) - c_gamma * abs(t)^g

  } else if (type == 2) {
    # Type 2: f(t) = d * |t|^{gamma - 1} * exp(-c * |t|), gamma ≥ 1
    c_gamma <- (gamma(g + 2) / gamma(g))^(1/2)
    d_gamma <- (1 / (2 * gamma(g))) * (gamma(g+ 2) / gamma(g))^(g/ 2)
    log_f <- log(d_gamma) + (g - 1) * log(abs(t)) - c_gamma * abs(t)
    
  } else if (type == 3) {
    # Type 3: f(t) = d * |t|^{gamma - 1} * exp(-c * |t|^gamma), gamma ≥ 1
    c_gamma <- gamma(2 / g + 1)^(g / 2)
    d_gamma <- (g / 2) * gamma(2 / g + 1)^(g / 2)
    log_f <- log(d_gamma) + (g - 1) * log(abs(t)) - c_gamma * abs(t)^g
  }

  return(-sum(log_f))
}

# Function to compute MLE of beta
compute_mle_beta <- function(data, g, type = 1) {
  res <- optimize(neg_loglik,
                  data = data, gamma = g, type = type)
  return(res$minimum)
}

## Making the wole thing work...
dati <- datagen()