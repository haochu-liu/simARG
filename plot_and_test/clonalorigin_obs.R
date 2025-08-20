rho_site <- 0.02
delta <- 300
theta_site <- 0.05
df_obs <- data.frame(r=rep(NA, 3),
                     g=rep(NA, 3),
                     s=rep(NA, 3))

set.seed(100)
tree <- clonal_genealogy(15L)

v_r <- rep(NA, 2000)
v_g3 <- rep(NA, 2000)
v_s <- rep(NA, 2000)
for (i in 1:2000) {
  ARG <- ClonalOrigin_pair_seq(tree, rho_site, 1e6L, delta, 50L)
  ARG_mutated <- FSM_mutation(ARG, theta_site, binary=TRUE)
  mat <- ARG_mutated$node_site[1:15, ]

  v_r[j] <- LD_r(mat)
  v_g3[j] <- G3_test(mat)
  v_s[j] <- any(as.logical(mat[, 1]))
}
df_obs$r[1] <- mean(v_r)
df_obs$g[1] <- mean(v_g3)
df_obs$s[1] <- mean(v_s)

v_r <- rep(NA, 2000)
v_g3 <- rep(NA, 2000)
v_s <- rep(NA, 2000)
for (i in 1:2000) {
  ARG <- ClonalOrigin_pair_seq(tree, rho_site, 1e6L, delta, 200L)
  ARG_mutated <- FSM_mutation(ARG, theta_site, binary=TRUE)
  mat <- ARG_mutated$node_site[1:15, ]

  v_r[j] <- LD_r(mat)
  v_g3[j] <- G3_test(mat)
  v_s[j] <- any(as.logical(mat[, 1]))
}
df_obs$r[2] <- mean(v_r)
df_obs$g[2] <- mean(v_g3)
df_obs$s[2] <- mean(v_s)

v_r <- rep(NA, 2000)
v_g3 <- rep(NA, 2000)
v_s <- rep(NA, 2000)
for (i in 1:2000) {
  ARG <- ClonalOrigin_pair_seq(tree, rho_site, 1e6L, delta, 2000L)
  ARG_mutated <- FSM_mutation(ARG, theta_site, binary=TRUE)
  mat <- ARG_mutated$node_site[1:15, ]

  v_r[j] <- LD_r(mat)
  v_g3[j] <- G3_test(mat)
  v_s[j] <- any(as.logical(mat[, 1]))
}
df_obs$r[3] <- mean(v_r)
df_obs$g[3] <- mean(v_g3)
df_obs$s[3] <- mean(v_s)
