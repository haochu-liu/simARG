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
  theta_site <- 0.1
  if (df_s50$type[i] == "rho") {
    rho_site <- df_s50$x[i]
  } else if (df_s50$type[i] == "delta") {
    delta <- df_s50$x[i]
  } else if (df_s50$type[i] == "theta") {
    theta_site <- df_s50$x[i]
  }

  v_LD <- rep(NA, 100)
  v_G3 <- rep(NA, 100)
  v_s <- rep(NA, 100)
  for (j in 1:100) {
    ARG <- ClonalOrigin_pair_seq(tree, rho_site, 1e6L, delta, 50L)
    ARG_mutated <- FSM_mutation(ARG, theta_site, binary=TRUE)
    mat <- ARG_mutated$node_site[1:15, ]
    tree1 <- local_tree(ARG, 1L)
    length1 <- sum(tree1$edge[, 3])

    v_LD[j] <- LD_r(mat)
    v_G3[j] <- G3_test(mat)
    v_s[j] <- 1 - exp(-(theta_site * length1) / 2)
  }

  df_s50$r[i] <- mean(v_LD)
  df_s50$g3[i] <- mean(v_G3)
  df_s50$s[i] <- mean(v_s)

  if (i%%20 == 0) {print(paste("Complete", i, "iterations"))}
}


