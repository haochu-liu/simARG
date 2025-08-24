for (i in 1:1000) {
  set.seed(i)
  theta_0 <- rep(NA, 3)
  theta_0[1] <- runif(1, min=0, max=0.2)
  theta_0[2] <- runif(1, min=1, max=2000)
  theta_0[3] <- runif(1, min=0, max=0.2)
  if (all(abs(theta_0 - c(0.02, 300, 0.05)) <= c(0.01, 100, 0.01))) {
    print(i)
  }
}

library(simARG)


rho_site <- 0.02
delta <- 300
theta_site <- 0.05
rho_site <- theta_0[1]
delta <- theta_0[2]
theta_site <- theta_0[3]
df_obs <- data.frame(r50=rep(NA, 10),
                     r200=rep(NA, 10),
                     r2000=rep(NA, 10),
                     g50=rep(NA, 10),
                     g200=rep(NA, 10),
                     g2000=rep(NA, 10),
                     s=rep(NA, 10))

set.seed(100)
tree <- clonal_genealogy(15L)

set.seed(1)
for (j in 1:10) {
  print(j)
  v_r <- rep(NA, 2000)
  v_g3 <- rep(NA, 2000)
  v_s <- rep(NA, 2000)
  obs_s <- rep(NA, 3)
  for (i in 1:2000) {
    ARG <- ClonalOrigin_pair_seq(tree, rho_site, 1e6L, delta, 50L)
    ARG_mutated <- FSM_mutation(ARG, theta_site, binary=TRUE)
    mat <- ARG_mutated$node_site[1:15, ]

    v_r[i] <- LD_r(mat)
    v_g3[i] <- G3_test(mat)
    v_s[i] <- any(as.logical(mat[, 1]))
  }
  df_obs$r50[j] <- mean(v_r)
  df_obs$g50[j] <- mean(v_g3)
  obs_s[1] <- mean(v_s)

  print("Complete k=50")

  v_r <- rep(NA, 2000)
  v_g3 <- rep(NA, 2000)
  v_s <- rep(NA, 2000)
  for (i in 1:2000) {
    ARG <- ClonalOrigin_pair_seq(tree, rho_site, 1e6L, delta, 200L)
    ARG_mutated <- FSM_mutation(ARG, theta_site, binary=TRUE)
    mat <- ARG_mutated$node_site[1:15, ]

    v_r[i] <- LD_r(mat)
    v_g3[i] <- G3_test(mat)
    v_s[i] <- any(as.logical(mat[, 1]))
  }
  df_obs$r200[j] <- mean(v_r)
  df_obs$g200[j] <- mean(v_g3)
  obs_s[2] <- mean(v_s)

  print("Complete k=200")

  v_r <- rep(NA, 2000)
  v_g3 <- rep(NA, 2000)
  v_s <- rep(NA, 2000)
  for (i in 1:2000) {
    ARG <- ClonalOrigin_pair_seq(tree, rho_site, 1e6L, delta, 2000L)
    ARG_mutated <- FSM_mutation(ARG, theta_site, binary=TRUE)
    mat <- ARG_mutated$node_site[1:15, ]

    v_r[i] <- LD_r(mat)
    v_g3[i] <- G3_test(mat)
    v_s[i] <- any(as.logical(mat[, 1]))
  }
  df_obs$r2000[j] <- mean(v_r)
  df_obs$g2000[j] <- mean(v_g3)
  obs_s[3] <- mean(v_s)
  df_obs$s[j] <- mean(obs_s)

  print("Complete k=2000")
  save(df_obs, file = "df_obs.rda")
}

save(df_obs, file = "df_obs.rda")

sigma_s <- diag(7)
k_0 <- kernel_func(obs, obs, tol, sigma_s)
print(exp(k_0))
tol <- 0.02
for (i in 1:10) {
  s_1 <- as.numeric(as.vector(df_obs[i, ]))
  k_1 <- kernel_func(obs, s_1, tol, sigma_s)
  print(exp(k_1 - k_0))
}

tol <- 0.02
for (i in 1:20) {
  s_1 <- as.numeric(as.vector(df_obs_0[i, ]))
  k_1 <- kernel_func(obs, s_1, tol, sigma_s)
  print(exp(k_1 - k_0))
}



