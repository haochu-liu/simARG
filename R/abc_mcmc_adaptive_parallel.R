#' Adaptive ABC-MCMC
#'
#' Run ABC-MCMC with adaptive Metropolis algorithm.
#'
#' @param obs The vector of observed data point.
#' @param tol A positive numeric value for the tolerance.
#' @param kernel_func A kernel function.
#' @param p_s_parallel A function to provide the sampled summary statistics by given parameters.
#' @param prior A function to provide the log density of prior.
#' @param theta_0 Initial theta from the prior.
#' @param s_0 A matrix of the initial summary statistics.
#' @param n_adapt Number of initial period.
#' @param n_iter Number of iterations.
#' @param sigma_s The covariance matrix for kernel function.
#' @param sigma_0 The covariance matrix for initial period.
#' @param sigma_epi Term to control the noise.
#' @param string A string for file saving.
#' @return Function value.
#' @export
abc_mcmc_adaptive_parallel <- function(obs, tol, kernel_func, p_s_parallel, prior,
                                       theta_0, s_0, n_adapt, n_iter,
                                       sigma_s, sigma_0, sigma_epi, string) {
  d <- length(theta_0)
  # Variables for output
  theta_matrix <- matrix(NA, nrow=(n_iter+1), ncol=d)
  s_matrix <- matrix(NA, nrow=(n_iter+1), ncol=length(obs))
  accept_vec <- rep(FALSE, n_iter+1)
  file_str <- paste0("abc_mat_", string, ".rda")

  # Initialization
  k_0 <- kernel_func(obs, s_0, tol, sigma_s)

  theta_matrix[1, ] <- theta_0
  s_matrix[1, ] <- s_0

  for (i in 2:n_adapt) {
    theta_1 <- rtruncnorm(n = 1, a = c(0, 1, 0), b = c(0.2, 2000, 0.2),
                          mean = theta_0, sd = diag(sigma_0))

    s_1 <- as.vector(p_s_parallel(theta_1,
                                  tree, ClonalOrigin_pair_seq, FSM_mutation, LD_r, G3_test))
    k_1 <- kernel_func(obs, s_1, tol, sigma_s)

    log_alpha <- k_1 + prior(theta_1) - k_0 - prior(theta_0)
    if (log(runif(1)) < log_alpha) {
      theta_0 <- theta_1
      s_0 <- s_1
      accept_vec[i] <- TRUE
    }

    theta_matrix[i, ] <- theta_0
    s_matrix[i, ] <- s_0
    print(paste("Complete", i, "iterations"))
    if (i%%10 == 0) {
      abc_mat <- list(theta_matrix=theta_matrix,
                      s_matrix=s_matrix,
                      accept_vec=accept_vec)
      save(abc_mat, file = file_str)
    }
  }

  s_d <- 2.38^2 / d
  mean_old <- as.matrix(colMeans(theta_matrix[1:n_adapt, ]))
  cov_sigma <- s_d*cov(theta_matrix[1:n_adapt, ]) + s_d * sigma_epi * diag(d)
  for (i in (n_adapt+1):(n_iter+1)) {
    theta_1 <- rtruncnorm(n = 1, a = c(0, 1, 0), b = c(0.2, 2000, 0.2),
                          mean = theta_0, sd = diag(cov_sigma))

    s_1 <- as.vector(p_s_parallel(theta_1,
                                  tree, ClonalOrigin_pair_seq, FSM_mutation, LD_r, G3_test))
    k_1 <- kernel_func(obs, s_1, tol, sigma_s)

    log_alpha <- k_1 + prior(theta_1) - k_0 - prior(theta_0)
    if (log(runif(1)) < log_alpha) {
      theta_0 <- theta_1
      s_0 <- s_1
      accept_vec[i] <- TRUE
    }

    theta_matrix[i, ] <- theta_0
    s_matrix[i, ] <- s_0

    m_theta <- as.matrix(theta_0, )
    mean_new <- mean_old*(i-1)/i + m_theta/i
    cov_sigma <- cov_sigma * (i-1) / i +
      s_d / i * (i * mean_old %*% t(mean_old) -
                   (i+1) * mean_new %*% t(mean_new) +
                   m_theta %*% t(m_theta) +
                   sigma_epi * diag(d))
    print(paste("Complete", i, "iterations"))
    if (i%%10 == 0) {
      abc_mat <- list(theta_matrix=theta_matrix,
                      s_matrix=s_matrix,
                      accept_vec=accept_vec)
      save(abc_mat, file = file_str)
    }
  }

  return(list(theta_matrix=theta_matrix,
              s_matrix=s_matrix,
              accept_vec=accept_vec))
}
