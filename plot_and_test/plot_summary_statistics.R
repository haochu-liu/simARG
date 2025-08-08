library(ggplot2)
library(patchwork)


tree <- clonal_genealogy(15L)
rho_seq <- seq(from=0, to=0.4, length.out=201)
delta_seq <- seq(from=20, to=4020, length.out=201)
theta_seq <- seq(from=0, to=0.3, length.out=201)
df_s50 <- data.frame(r=rep(NA, 201*3),
                     g3=rep(NA, 201*3),
                     s=rep(NA, 201*3),
                     x=c(rho_seq, delta_seq, theta_seq),
                     type=c(rep("rho", 201),
                            rep("delta", 201),
                            rep("theta", 201)))
df_s200 <- data.frame(r=rep(NA, 201*3),
                      g3=rep(NA, 201*3),
                      s=rep(NA, 201*3),
                      x=c(rho_seq, delta_seq, theta_seq),
                      type=c(rep("rho", 201),
                             rep("delta", 201),
                             rep("theta", 201)))
df_s2000 <- data.frame(r=rep(NA, 201*3),
                       g3=rep(NA, 201*3),
                       s=rep(NA, 201*3),
                       x=c(rho_seq, delta_seq, theta_seq),
                       type=c(rep("rho", 201),
                              rep("delta", 201),
                              rep("theta", 201)))


# k = 50
for (i in 1:nrow(df_s50)) {
  rho_site <- 0.02
  delta <- 300
  theta_site <- 0.05
  if (df_s50$type[i] == "rho") {
    rho_site <- df_s50$x[i]
  } else if (df_s50$type[i] == "delta") {
    delta <- df_s50$x[i]
  } else if (df_s50$type[i] == "theta") {
    theta_site <- df_s50$x[i]
  }
  ARG <- ClonalOrigin_pair_seq(tree, rho_site, 1e6L, delta, 50L)
  ARG_mutated <- FSM_mutation(ARG, theta_site, binary=TRUE)
  mat <- ARG_mutated$node_site[1:15, ]
  tree1 <- local_tree(ARG, 1L)
  height1 <- tree1$sum_time

  df_s50$r[i] <- LD_r(mat)
  df_s50$g3[i] <- G3_test(mat)
  df_s50$s[i] <- 1 - exp(-(theta_site * height1) / 2)

  # if (i%%100 == 0) {print(paste("Complete", i, "iterations"))}
  print(paste("Complete", i, "iterations"))
}


