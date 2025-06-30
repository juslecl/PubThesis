set.seed(123)
library(vcvComp)
library(stats)
library(optimx)
library(BiocParallel)
library(ggplot2)
library(purrr)
library(dplyr)
library(tidyr)

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
datagen <- function(n, obj, rcdf, g, meanX, sdX,  M=1000) {
  bita <- obj[-length(obj)]
  s <- obj[length(obj)]
  d <- length(meanX)
  datasets <- vector("list", M) 
  for (m in 1:M) {  
    x <- matrix(data = NA, nrow = n, ncol = d)  
    for (j in 1:d)  x[,j] <- rnorm(n,mean=meanX[j],sd=sdX[j])
    eps <- rcdf(n,g,s)  
    y <- x%*%bita + eps  
    datasets[[m]] <- list(Y = y, X = x, P = obj, E = eps)  
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
neg_loglik <- function(obj, data, g, type = 1) { #param is a vector of length(bita)+1
  t <- (data$Y-data$X%*%obj[1:(length(obj)-1)])/obj[length(obj)]
  s <- obj[length(obj)]
  
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

# Function to compute MLE of bita. 
compute_mle_bita <- function(data, g, type = 1) {
  ols_est <- ols(data)
  s_start <- sqrt(var(data$Y - data$X%*%ols_est))
  neg_loglik_wrapped <- function(obj) {
    neg_loglik(obj, data = data, g = g, type = type)
  }
  res <- optimx(par = c(as.numeric(ols_est), s_start), fn = neg_loglik_wrapped, method = "BFGS")
  resvec <- res[, 1:(ncol(data$X))]
  return(mle=t(resvec))
}

#-------------------------------------------------------------
#-------------------------------------------------------------
mleols <- function(data, g, type){
  param <- BatchtoolsParam(workers=32)
  M <- length(data)
  res <- bptry({bplapply(seq_len(M), function(nb){
    mle <- compute_mle_bita(data[[nb]], g, type)
    olse <- ols(data[[nb]])
    return(cbind.data.frame(mle=mle,ols=olse,truth=data[[nb]]$P[-length(data[[nb]]$P)]))
  },BPPARAM=param)})
  ordremle <- matrix(unlist(lapply(res,function(i) return(i[,1]))),ncol=ncol(data[[1]]$X),byrow=T)
  ordreols <- matrix(unlist(lapply(res,function(i) return(i[,2]))),ncol=ncol(data[[1]]$X),byrow=T)
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
  return(list(Norm=norms,MEANS=empirical_means,VARS=empirical_vars,MSE=apply(norms,2,function(i) mean(i^2))))
}

#-------------------------------------------------------------
#-------------------------------------------------------------
run_simulation <- function(d, n, g) {
  set.seed(123)
  bita <- c(sample(seq(from = -3, to = 3, by = .01), d),
            sample(seq(from = 1, to = 3, by = .01), 1))
  meanX <- sample(seq(from = -5, to = 5, by = .01), d)
  sdX <- sample(seq(from = 0, to = 2, by = .1), d)
  
  model_list <- list(rmod1, rmod2, rmod3)
  model_names <- paste0("Model ", 1:3)
  
  results <- list()
  
  for (i in seq_along(model_list)) {
    rmod <- model_list[[i]]
    data <- datagen(n, bita, rmod, g, meanX, sdX)
    mleols_res <- mleols(data, g, i)
    metric <- metrics(mleols_res)
    
    results[[i]] <- list(
      model = model_names[i],
      metric = metric,
      mle = mleols_res[[2]],
      ols = mleols_res[[3]]
    )
  }
  
  metric_df <- lapply(seq_along(results), function(i) {
    data.frame(
      Norm_OLS = results[[i]]$metric[[1]][,1],
      Norm_MLE = results[[i]]$metric[[1]][,2],
      Model = model_names[i]
    )
  }) %>% bind_rows()
  
  metric_df_long <- pivot_longer(metric_df,
                                 cols = c(Norm_OLS, Norm_MLE),
                                 names_to = "Estimator",
                                 values_to = "Norm")
  
  metric_df_long$Estimator <- gsub("Norm_", "", metric_df_long$Estimator)
  
  p1 <- ggplot(metric_df_long, aes(x = Model, y = Norm, fill = Estimator)) +
    geom_boxplot(outlier.shape = NA, position = position_dodge(0.8)) +
    labs(y = "Euclidean norm", x = "", title = "Comparison of the distance to the true objeter") +
    theme_minimal(base_size = 14) +
    theme(legend.position = "top") +
    scale_fill_manual(values = c("OLS" = "blue", "MLE" = "red"))
  
  boxplot_filename <- sprintf("boxplot_d%d_n%d_gamma%.2f.pdf", d, n, g)
  ggsave(boxplot_filename, p1, width = 8, height = 6)
  
  bita_index <- 2
  density_df <- lapply(seq_along(results), function(i) {
    mle_density <- density(results[[i]]$mle[, bita_index])
    ols_density <- density(results[[i]]$ols[, bita_index])
    data.frame(
      x = c(mle_density$x, ols_density$x),
      y = c(mle_density$y, ols_density$y),
      Estimator = rep(c("MLE", "OLS"), each = length(mle_density$x)),
      Model = results[[i]]$model
    )
  }) %>% bind_rows()
  
  p2 <- ggplot(density_df, aes(x = x, y = y, color = Estimator)) +
    geom_line(linewidth = 1.1) +
    facet_wrap(~Model, scales = "free", nrow = 1) +
    labs(x = expression(hat(beta)[2]), y = "Density", title = "Comparison of the estimators'densities (Component 2)") +
    theme_minimal(base_size = 14) +
    scale_color_manual(values = c("MLE" = "red", "OLS" = "blue")) +
    theme(legend.position = "top")
  
  density_filename <- sprintf("density_d%d_n%d_gamma%.2f.pdf", d, n, g)
  ggsave(density_filename, p2, width = 10, height = 4)
  
  return(list(
    bita = bita,
    meanX = meanX,
    sdX = sdX,
    results = results,
    metric_df = metric_df,
    density_df = density_df,
    boxplot_file = boxplot_filename,
    density_file = density_filename
  ))
}


sim1 <- run_simulation(5,500,3)
sim2 <- run_simulation(5,500,5)
sim3 <- run_simulation(5,500,7)
sim4 <- run_simulation(5,1000,3)
sim5 <- run_simulation(5,1000,5)
sim6 <- run_simulation(5,1000,7)
sim7 <- run_simulation(15,2000,3)
sim8 <- run_simulation(15,2000,5)
sim9 <- run_simulation(15,2000,7)

#-------------------------------------------------------------
#-------------------------------------------------------------
mse <- c(0.236 , 0.214 , 0.236 , 0.156 , 0.236 , 0.119, 0.126 , 0.111 , 0.126 , 0.0768 , 0.126 , 0.057 ,11.8 , 10.0 , 11.8 , 6.64 , 11.9 , 4.79,0.259 , 0.126 , 0.246 , 0.0614 , 0.220 , 0.0264,0.117 , 0.0435 ,0.121 , 0.024 , 0.121 , 0.0122,10.9 , 4.47 , 10.7 , 1.67 , 11.5 , 1.18,0.247 , 0.0481 , 0.264 , 0.0135 , 0.256 ,  0.00603,0.120 , 0.0182 , 0.120 , 0.00555 , 0.122 ,0.00292,11.0 , 1.78 , 10.9 , 0.54  , 11.1 , 0.27)
est <- rep(c("OLSE","MLE"),27)
fam <- rep(1:3,each=18)
nd <- rep(rep(c(100,200,400/3),each=6),3)
gamma <- rep(rep(c(3,5,7),each=2),9)

# First, make a data frame
df1 <- data.frame(
  nd = as.numeric(nd),
  mse = mse,
  est = as.factor(est),
  gamma = gamma,
  fam = paste("Family", fam)
)

# ggplot version
p_mse_nd <- ggplot(df1, aes(x = nd, y = mse, color = est)) +
  geom_point(size = 2) +
  geom_text_repel(aes(label = gamma), size = 3, color = "black", max.overlaps = Inf) +
  facet_wrap(~fam, nrow = 1) +
  scale_color_brewer(palette = "Dark2") +
  labs(x = "n/d", y = "MSE", color = "Estimator",
       title = "MSE vs n/d for each family") +
  ylim(c(0,15)) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "top",
    strip.text = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5)
  )

ggsave("mse_vs_nd_faceted.pdf", p_mse_nd, width = 10, height = 4)

# --------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------

nd <- rep(rep(c(100,200,"400/3"),each=6),3)

df2 <- data.frame(
  nd = as.character(nd),
  mse = mse,
  est = as.factor(est),
  gamma = gamma,
  fam = paste("Family", fam)
)

p_mse_gamma <- ggplot(df2, aes(x = gamma, y = mse, color = est)) +
  geom_point(size = 2) +
  geom_text_repel(aes(label = nd), size = 3, color = "black", max.overlaps = Inf) +
  facet_wrap(~fam, nrow = 1) +
  scale_color_brewer(palette = "Dark2") +
  labs(x = expression(gamma), y = "MSE", color = "Estimator",
       title = "MSE vs gamma for each family") +
  ylim(c(0,15)) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "top",
    strip.text = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5)
  )

ggsave("mse_vs_gamma_faceted.pdf", p_mse_gamma, width = 10, height = 4)

# --------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------

# Plot the variance efficiency 
d <- 5

# Function 1: eta_1
eta1 <- function(gamma, d) {
  inside <- (gamma(1/gamma)^2) / (gamma^2 * gamma(3/gamma) * gamma(2 - 1/gamma))
  inside^(1/d)
}

# Function 2: eta_2, only defined for gamma > 2
eta2 <- function(gamma, d) {
  inside <- (gamma(gamma) * (gamma - 2)) / gamma(gamma + 2)
  inside^(1/d)
}

# Function 3: eta_3, only defined for gamma > 2
eta3 <- function(gamma, d) {
  inside <- 1 / (gamma(2/gamma + 1) * gamma(1 - 2/gamma) * (gamma - 1)^2)
  inside^(1/d)
}

# Create gamma sequences for respective valid domains
gamma_seq_eta1 <- seq(0.6, 10, length.out = 1000)
gamma_seq_eta2_eta3 <- seq(2.01, 10, length.out = 1000)

# Evaluate eta values
df_eta1 <- data.frame(
  gamma = gamma_seq_eta1,
  eta = eta1(gamma_seq_eta1, d)
)

df_eta2 <- data.frame(
  gamma = gamma_seq_eta2_eta3,
  eta = eta2(gamma_seq_eta2_eta3, d)
)

df_eta3 <- data.frame(
  gamma = gamma_seq_eta2_eta3,
  eta = eta3(gamma_seq_eta2_eta3, d)
)

# Plot eta_1
p1 <- ggplot(df_eta1, aes(x = gamma, y = eta)) +
  geom_line(color = "darkred", size = 1.1) +
  labs(
    title = "Asymptotic relative efficiency for family 1",
    subtitle = bquote(d == .(d)),
    x = expression(gamma), y = expression(eta[1](gamma))
  ) +
  theme_minimal(base_size = 14)

# Plot eta_2
p2 <- ggplot(df_eta2, aes(x = gamma, y = eta)) +
  geom_line(color = "darkblue", size = 1.1) +
  labs(
    title = "Asymptotic relative efficiency for family 2",
    subtitle = bquote(d == .(d)),
    x = expression(gamma), y = expression(eta[2](gamma))
  ) +
  theme_minimal(base_size = 14)

# Plot eta_3
p3 <- ggplot(df_eta3, aes(x = gamma, y = eta)) +
  geom_line(color = "forestgreen", size = 1.1) +
  labs(
    title = "Asymptotic relative efficiency for family 3",
    subtitle = bquote(d == .(d)),
    x = expression(gamma), y = expression(eta[3](gamma))
  ) +
  theme_minimal(base_size = 14)

# Save each plot
ggsave(sprintf("eta1_d%d.pdf", d), p1, width = 7, height = 5)
ggsave(sprintf("eta2_d%d.pdf", d), p2, width = 7, height = 5)
ggsave(sprintf("eta3_d%d.pdf", d), p3, width = 7, height = 5)





