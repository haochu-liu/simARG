library(ggplot2)
library(mvtnorm)


set.seed(100)
tree <- clonal_genealogy(15L)
tol <- 0.1
# s_obs <- c(df_obs$r, df_obs$g, mean(df_obs$s))
s_obs <- c(0.009291352, 0.005619992, 0.002516754,
           0.002000000, 0.008000000, 0.007500000,
           0.201166667)

# model sample
p_s <- function(theta) {
  rho_site <- theta[1]
  delta <- theta[2]
  theta_site <- theta[3]
  s_vec <- rep(NA, 7)

  v_r <- rep(NA, 2000)
  v_g3 <- rep(NA, 2000)
  v_s <- rep(NA, 6000)
  for (i in 1:2000) {
    ARG <- ClonalOrigin_pair_seq(tree, rho_site, 1e6L, delta, 50L)
    ARG_mutated <- FSM_mutation(ARG, theta_site, binary=TRUE)
    mat <- ARG_mutated$node_site[1:15, ]

    v_r[i] <- LD_r(mat)
    v_g3[i] <- G3_test(mat)
    v_s[i] <- any(as.logical(mat[, 1]))
  }
  s_vec[1] <- mean(v_r)
  s_vec[4] <- mean(v_g3)

  v_r <- rep(NA, 2000)
  v_g3 <- rep(NA, 2000)
  for (i in 1:2000) {
    ARG <- ClonalOrigin_pair_seq(tree, rho_site, 1e6L, delta, 200L)
    ARG_mutated <- FSM_mutation(ARG, theta_site, binary=TRUE)
    mat <- ARG_mutated$node_site[1:15, ]

    v_r[i] <- LD_r(mat)
    v_g3[i] <- G3_test(mat)
    v_s[i+2000] <- any(as.logical(mat[, 1]))
  }
  s_vec[2] <- mean(v_r)
  s_vec[5] <- mean(v_g3)

  v_r <- rep(NA, 2000)
  v_g3 <- rep(NA, 2000)
  for (i in 1:2000) {
    ARG <- ClonalOrigin_pair_seq(tree, rho_site, 1e6L, delta, 2000L)
    ARG_mutated <- FSM_mutation(ARG, theta_site, binary=TRUE)
    mat <- ARG_mutated$node_site[1:15, ]

    v_r[i] <- LD_r(mat)
    v_g3[i] <- G3_test(mat)
    v_s[i+4000] <- any(as.logical(mat[, 1]))
  }
  s_vec[3] <- mean(v_r)
  s_vec[6] <- mean(v_g3)
  s_vec[7] <- mean(v_s)

  return(s_vec)
}

# prior density
prior <- function(theta_0) {
  d_rho <- dunif(theta_0[1], min=0, max=0.2, log=TRUE)
  d_delta <- dunif(theta_0[2], min=1, max=2000, log=TRUE)
  d_theta <- dunif(theta_0[3], min=0, max=0.2, log=TRUE)

  return(d_rho + d_delta + d_theta)
}
sigma_s <- diag(7)
sigma_0 <- diag(c(0.002, 200, 0.002))

# initial state
repeat {
  theta_0 <- rep(NA, 3)
  theta_0[1] <- runif(1, min=0, max=0.2)
  theta_0[2] <- runif(1, min=1, max=2000)
  theta_0[3] <- runif(1, min=0, max=0.2)
  s_0 <- p_s(theta_0)

  k_0 <- gaussian_kernel(s_obs, s_0, tol, sigma_0)
  if (exp(k_0) > .Machine$double.eps) {break}
}

abc_mat <- abc_mcmc_adaptive(s_obs, tol, gaussian_kernel, p_s, prior,
                             theta_0, s_0, 100, 1000,
                             sigma_s, sigma_0, 0)

# hist of posterior
hist(abc_mat$theta_matrix[1000:2001, 2], probability = TRUE, main = "Histogram of mu|s_obs",
     breaks = 20, col = "gray", border = "black", xlab="mu")

# trace plot
df <- data.frame(x = 1:2001,
                 y = abc_mat$theta_matrix[, 2])
ggplot(df, aes(x = x, y = y)) +
  geom_line()
mean(abc_mat$accept_vec)
