library(ggplot2)
library(mvtnorm)


mu <- c(1, -1)
# observation
s_obs <- colMeans(rmvnorm(n=500, mean=mu, sigma=diag(c(1, 1))))

# proposal function
p_theta <- function(theta_0) {
  rmvnorm(n=1, mean=theta_0, sigma=diag(c(0.5, 0.5)))
}
d_theta <- function(theta_1, theta_0) {
  dmvnorm(theta_1, mean=theta_0, sigma=diag(c(0.5, 0.5)), log=TRUE)
}
# model sample
p_s <- function(theta) {
  rmvnorm(n=1, mean=theta, sigma=diag(c(1, 1)))
}
# prior density
prior <- function(theta_0) {
  dmvnorm(theta_0, mean=c(0, 0), sigma=diag(c(10, 10)), log=TRUE)
}
# initial theta
mu_0 <- rmvnorm(n=1, mean=c(0, 0), sigma=diag(c(10, 10)))
# ABC-MCMC
matrix_list <- abc_mcmc(c(s_obs), 0.01, gaussian_kernel, p_theta, d_theta, p_s, prior, mu_0, 10000)
# hist of posterior
hist(matrix_list$theta_matrix[5000:10001, ], probability = TRUE, main = "Histogram of mu|s_obs",
     breaks = 20, col = "gray", border = "black", xlab="mu")
abline(v = s_obs, col = "red", lwd = 2, lty = 2)
# trace plot
df <- data.frame(x = 1:10001,
                 y = matrix_list$theta_matrix[, 2])
ggplot(df, aes(x = x, y = y)) +
  geom_line()


# adaptive ABC-MCMC
repeat {
  mu_0 <- rmvnorm(n=1, mean=c(0, 0), sigma=diag(c(10, 10)))
  s_0 <- as.vector(p_s(mu_0))

  k_0 <- gaussian_kernel(c(s_obs), s_0, 1, diag(2))
  if (exp(k_0) > .Machine$double.eps) {break}
}
mu_0 <- as.vector(mu_0)
matrix_list <- abc_mcmc_adaptive(c(s_obs), 1, gaussian_kernel, p_s, prior, mu_0, s_0, 100, 1000, diag(2), 0)
# hist of posterior
hist(matrix_list$theta_matrix[100:1001, ], probability = TRUE, main = "Histogram of mu|s_obs",
     breaks = 20, col = "gray", border = "black", xlab="mu")
abline(v = s_obs, col = "red", lwd = 2, lty = 2)
# trace plot
df <- data.frame(x = 1:1001,
                 y = matrix_list$theta_matrix[, 2])
ggplot(df, aes(x = x, y = y)) +
  geom_line()
mean(matrix_list$accept_vec)

