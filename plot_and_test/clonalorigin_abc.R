library(ggplot2)
library(patchwork)
library(latex2exp)
library(mvtnorm)
library(truncnorm)
library(foreach)
library(doParallel)


set.seed(100)
tree <- clonal_genealogy(15L)

tol <- 0.015
# s_obs <- c(df_obs$r, df_obs$g, mean(df_obs$s))
s_obs <- c(0.009291352, 0.005619992, 0.002516754,
           0.002000000, 0.008000000, 0.007500000,
           0.201166667)

# model sample
ClonalOrigin_pair_seq <- ClonalOrigin_pair_seq
FSM_mutation <- FSM_mutation
LD_r <- LD_r
G3_test <- G3_test

p_s_parallel <- function(theta,
                         tree, ClonalOrigin_pair_seq, FSM_mutation, LD_r, G3_test) {
  # Set up a parallel backend with 5 cores
  # 'makeCluster' creates a cluster of R processes on the local machine
  cl <- makeCluster(10)
  registerDoParallel(cl)

  rho_site <- theta[1]
  delta <- theta[2]
  theta_site <- theta[3]
  s_vec <- rep(NA, 7)



  # Function to run a single iteration for a given loop
  run_simulation <- function(rho_site, delta, theta_site, l_val) {
    ARG <- ClonalOrigin_pair_seq(tree, rho_site, 1e6L, delta, l_val)
    ARG_mutated <- FSM_mutation(ARG, theta_site, binary=TRUE)
    mat <- ARG_mutated$node_site[1:15, ]

    return(list(
      LD_r = LD_r(mat),
      G3_test = G3_test(mat),
      s = any(as.logical(mat[, 1]))
    ))
  }

  # Define loop lengths
  l_values <- c(50L, 200L, 2000L)

  # Use foreach to parallelize the three groups of simulations
  results <- foreach(
    l = l_values,
    .combine = 'list',
    .multicombine = TRUE,
    .packages = c("foreach") # Add any packages needed inside the loop
  ) %dopar% {

    # Run the 2000 simulations for each l_value
    l_results <- foreach(i = 1:2000, .combine = 'rbind') %do% {
      res <- run_simulation(rho_site, delta, theta_site, l)
      c(res$LD_r, res$G3_test, res$s)
    }

    # Convert to a data frame for easier manipulation
    l_results <- as.data.frame(l_results)
    colnames(l_results) <- c("v_r", "v_g3", "v_s")
    l_results
  }

  # Shut down the cluster
  stopCluster(cl)

  # Process the results
  # Note that this part is different since the results are structured differently
  s_vec[1] <- mean(results[[1]]$v_r)
  s_vec[4] <- mean(results[[1]]$v_g3)

  s_vec[2] <- mean(results[[2]]$v_r)
  s_vec[5] <- mean(results[[2]]$v_g3)

  s_vec[3] <- mean(results[[3]]$v_r)
  s_vec[6] <- mean(results[[3]]$v_g3)

  # s_vec[7] is a bit trickier, as v_s is a single vector of length 6000
  # We need to collect all v_s values from all results
  all_s_values <- c(results[[1]]$v_s, results[[2]]$v_s, results[[3]]$v_s)
  s_vec[7] <- mean(all_s_values)

  return(s_vec)
}

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
set.seed(809)
repeat {
  theta_0 <- rep(NA, 3)
  theta_0[1] <- runif(1, min=0, max=0.2)
  theta_0[2] <- runif(1, min=1, max=2000)
  theta_0[3] <- runif(1, min=0, max=0.2)
  s_0 <- p_s(theta_0)

  k_0 <- gaussian_kernel(s_obs, s_0, tol, sigma_s)
  if (exp(k_0) > .Machine$double.eps) {break}
}

abc_mat <- abc_mcmc_adaptive_parallel(s_obs, tol, gaussian_kernel, p_s_parallel, prior,
                                      theta_0, s_0, 100, 500,
                                      sigma_s, sigma_0, 0, "1")


# hist of posterior
hist(abc_mat$theta_matrix[501:5001, 3], probability = TRUE, main = "Histogram of mu|s_obs",
     breaks = 20, col = "gray", border = "black", xlab="mu")

prior_1 <- data.frame(
  x = c(0, 0.2),
  y = c(1/(0.2-0), 1/(0.2-0))
)

prior_2 <- data.frame(
  x = c(1, 2000),
  y = c(1/(2000-1), 1/(2000-1))
)

theta_df <- as.data.frame(abc_mat$theta_matrix[501:5001, ])
colnames(theta_df) <- c("rho_s", "delta", "theta_s")

hist1 <- ggplot(theta_df, aes(x = rho_s)) +
  # Add the posterior histogram, mapping the fill aesthetic to a fixed value
  geom_histogram(
    aes(y = ..density.., fill = "Posterior"),
    bins = 25,
    # binwidth = 0.0025,
    color = "black",
    alpha = 0.7
  ) +
  # Add the Prior density line, mapping the color aesthetic to a fixed value
  geom_line(
    data = prior_1,
    aes(x = x, y = y, color = "Prior"),
    size = 0.8
  ) +
  # Add the True value vertical line, mapping the color aesthetic to a fixed value
  geom_vline(
    xintercept = 0.02,
    size = 0.8,
    linetype = "solid",
    color = "green"
  ) +
  # Manually define the colors for the mapped aesthetics
  scale_color_manual(
    name = "",
    values = c("Posterior" = "darkblue", "Prior" = "red", "True value" = "green")
  ) +
  # Manually define the fill color for the histogram
  scale_fill_manual(
    name = "",
    values = "darkblue"
  ) +
  # Customize the theme and labels
  labs(
    x = TeX("$\\rho_s$"),
    y = "Density"
  ) +
  theme_minimal() +
  theme(legend.position = "none") +
  # Set axis limits
  xlim(0, 0.2)

hist2 <- ggplot(theta_df, aes(x = delta)) +
  # Add the posterior histogram, mapping the fill aesthetic to a fixed value
  geom_histogram(
    aes(y = ..density.., fill = "Posterior"),
    bins = 25,
    # binwidth = 10,
    color = "black",
    alpha = 0.7
  ) +
  # Add the Prior density line, mapping the color aesthetic to a fixed value
  geom_line(
    data = prior_2,
    aes(x = x, y = y, color = "Prior"),
    size = 0.8
  ) +
  # Add the True value vertical line, mapping the color aesthetic to a fixed value
  geom_vline(
    xintercept = 300,
    size = 0.8,
    linetype = "solid",
    color = "green"
  ) +
  # Manually define the colors for the mapped aesthetics
  scale_color_manual(
    name = "",
    values = c("Posterior" = "darkblue", "Prior" = "red", "True value" = "green")
  ) +
  # Manually define the fill color for the histogram
  scale_fill_manual(
    name = "",
    values = "darkblue"
  ) +
  # Customize the theme and labels
  labs(
    x = TeX("$\\delta$"),
    y = "Density"
  ) +
  theme_minimal() +
  theme(legend.position = "none") +
  # Set axis limits
  xlim(1, 2000)

hist3 <- ggplot(theta_df, aes(x = theta_s)) +
  # Add the posterior histogram, mapping the fill aesthetic to a fixed value
  geom_histogram(
    aes(y = ..density.., fill = "Posterior"),
    bins = 25,
    # binwidth = 0.005,
    color = "black",
    alpha = 0.7
  ) +
  # Add the Prior density line, mapping the color aesthetic to a fixed value
  geom_line(
    data = prior_1,
    aes(x = x, y = y, color = "Prior"),
    size = 0.8
  ) +
  # Add the True value vertical line, mapping the color aesthetic to a fixed value
  geom_vline(
    xintercept = 0.05,
    size = 0.8,
    linetype = "solid",
    color = "green"
  ) +
  # Manually define the colors for the mapped aesthetics
  scale_color_manual(
    name = "",
    values = c("Posterior" = "darkblue", "Prior" = "red", "True value" = "green")
  ) +
  # Manually define the fill color for the histogram
  scale_fill_manual(
    name = "",
    values = "darkblue"
  ) +
  # Customize the theme and labels
  labs(
    x = TeX("$\\theta_s$"),
    y = "Density"
  ) +
  theme_minimal() +
  theme(legend.position = "right") +
  # Set axis limits
  xlim(0, 0.2)

combined_hist <- hist1 + hist2 + hist3
combined_hist <- combined_hist +
  plot_annotation(
    title = "Estimated marginal posterior densities of the parameters")
print(combined_hist)


# trace plot
theta_df <- as.data.frame(abc_mat$theta_matrix[1:5001, ])
colnames(theta_df) <- c("rho_s", "delta", "theta_s")
theta_df$x <- 1:5001
trace1 <- ggplot(theta_df, aes(x = x, y = rho_s)) +
  geom_line() +
  labs(
    x = "Iterations",
    y = TeX("$\\rho_s$")
  ) +
  theme_minimal()
trace2 <- ggplot(theta_df, aes(x = x, y = delta)) +
  geom_line() +
  labs(
    x = "Iterations",
    y = TeX("$\\delta$")
  ) +
  theme_minimal()
trace3 <- ggplot(theta_df, aes(x = x, y = theta_s)) +
  geom_line() +
  labs(
    x = "Iterations",
    y = TeX("$\\theta_s$")
  ) +
  theme_minimal()

combined_trace <- trace1 / trace2 / trace3
combined_trace <- combined_trace +
  plot_annotation(
    title = "Trace plot of ABC-MCMC with AM algorithm")
print(combined_trace)


quantile(abc_mat$theta_matrix[501:5001, 3], probs = c(0.025, 0.975))
mean(abc_mat$theta_matrix[501:5001, 3])
mean(abc_mat$accept_vec[1:5001])
