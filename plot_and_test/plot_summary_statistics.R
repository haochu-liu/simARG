library(ggplot2)
library(patchwork)


tree <- clonal_genealogy(15L)
rho_seq <- seq(from=0, to=0.4, length.out=2001)
delta_seq <- seq(from=20, to=4020, length.out=2001)
theta_seq <- seq(from=0.2, to=1.2, length.out=2001)
df_s50 <- data.frame(r=rep(NA, 2001*3),
                     g3=rep(NA, 2001*3),
                     s=rep(NA, 2001*3),
                     x=c(rho_seq, delta_seq, theta_seq),
                     type=c(rep("rho", 2001),
                            rep("delta", 2001),
                            rep("theta", 2001)))
df_s200 <- data.frame(r=rep(NA, 2001*3),
                      g3=rep(NA, 2001*3),
                      s=rep(NA, 2001*3),
                      x=c(rho_seq, delta_seq, theta_seq),
                      type=c(rep("rho", 2001),
                             rep("delta", 2001),
                             rep("theta", 2001)))
df_s2000 <- data.frame(r=rep(NA, 2001*3),
                       g3=rep(NA, 2001*3),
                       s=rep(NA, 2001*3),
                       x=c(rho_seq, delta_seq, theta_seq),
                       type=c(rep("rho", 2001),
                              rep("delta", 2001),
                              rep("theta", 2001)))


# k = 50
for (i in 1:nrow(df_s50)) {
  rho_site <- 0.02
  delta <- 300
  theta_site <- 0.5
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
  length1 <- sum(tree1$edge[, 3])

  df_s50$r[i] <- LD_r(mat)
  df_s50$g3[i] <- G3_test(mat)
  df_s50$s[i] <- 1 - exp(-(theta_site * length1) / 2)

  if (i%%20 == 0) {print(paste("Complete", i, "iterations"))}
}

# k = 200
for (i in 1:nrow(df_s200)) {
  rho_site <- 0.02
  delta <- 300
  theta_site <- 0.5
  if (df_s200$type[i] == "rho") {
    rho_site <- df_s200$x[i]
  } else if (df_s200$type[i] == "delta") {
    delta <- df_s200$x[i]
  } else if (df_s200$type[i] == "theta") {
    theta_site <- df_s200$x[i]
  }

  ARG <- ClonalOrigin_pair_seq(tree, rho_site, 1e6L, delta, 200L)
  ARG_mutated <- FSM_mutation(ARG, theta_site, binary=TRUE)
  mat <- ARG_mutated$node_site[1:15, ]
  tree1 <- local_tree(ARG, 1L)
  length1 <- sum(tree1$edge[, 3])

  df_s200$r[i] <- LD_r(mat)
  df_s200$g3[i] <- G3_test(mat)
  df_s200$s[i] <- 1 - exp(-(theta_site * length1) / 2)

  if (i%%20 == 0) {print(paste("Complete", i, "iterations"))}
}

# k = 2000
for (i in 1:nrow(df_s2000)) {
  rho_site <- 0.02
  delta <- 300
  theta_site <- 0.5
  if (df_s2000$type[i] == "rho") {
    rho_site <- df_s2000$x[i]
  } else if (df_s2000$type[i] == "delta") {
    delta <- df_s2000$x[i]
  } else if (df_s2000$type[i] == "theta") {
    theta_site <- df_s2000$x[i]
  }

  ARG <- ClonalOrigin_pair_seq(tree, rho_site, 1e6L, delta, 2000L)
  ARG_mutated <- FSM_mutation(ARG, theta_site, binary=TRUE)
  mat <- ARG_mutated$node_site[1:15, ]
  tree1 <- local_tree(ARG, 1L)
  length1 <- sum(tree1$edge[, 3])

  df_s2000$r[i] <- LD_r(mat)
  df_s2000$g3[i] <- G3_test(mat)
  df_s2000$s[i] <- 1 - exp(-(theta_site * length1) / 2)

  if (i%%20 == 0) {print(paste("Complete", i, "iterations"))}
}


df_s50$k <- "50"
df_s200$k <- "200"
df_s2000$k <- "2000"

df_r <- rbind(df_s50[, c(1, 4, 5, 6)], df_s200[, c(1, 4, 5, 6)])
df_r <- rbind(df_r, df_s2000[, c(1, 4, 5, 6)])
plot_r_rho <- ggplot(subset(df_r, type=="rho"), aes(x=x, y=r, color=k)) +
  geom_point(size = 2, shape = 4, stroke = 1.2, alpha=0.7) +
  labs(
    x = "rho",
    y = "r square",
    color = "distance"
  ) +
  theme_minimal()
plot_r_delta <- ggplot(subset(df_r, type=="delta"), aes(x=x, y=r, color=k)) +
  geom_point(size = 2, shape = 4, stroke = 1.2, alpha=0.7) +
  labs(
    x = "delta",
    y = "r square",
    color = "distance"
  ) +
  theme_minimal()
plot_r_theta <- ggplot(subset(df_r, type=="theta"), aes(x=x, y=r, color=k)) +
  geom_point(size = 2, shape = 4, stroke = 1.2, alpha=0.7) +
  labs(
    x = "theta",
    y = "r square",
    color = "distance"
  ) +
  theme_minimal()
plot_r_rho + plot_r_delta + plot_r_theta


df_g <- rbind(df_s50[, c(2, 4, 5, 6)], df_s200[, c(2, 4, 5, 6)])
df_g <- rbind(df_g, df_s2000[, c(2, 4, 5, 6)])
plot_g_rho <- ggplot(subset(df_g, type=="rho"), aes(x=x, y=g3, color=k)) +
  geom_point(size = 2, shape = 4, stroke = 1.2, alpha=0.7) +
  labs(
    x = "rho",
    y = "G4 test",
    color = "distance"
  ) +
  theme_minimal()
plot_g_delta <- ggplot(subset(df_g, type=="delta"), aes(x=x, y=g3, color=k)) +
  geom_point(size = 2, shape = 4, stroke = 1.2, alpha=0.7) +
  labs(
    x = "delta",
    y = "G4 test",
    color = "distance"
  ) +
  theme_minimal()
plot_g_theta <- ggplot(subset(df_g, type=="theta"), aes(x=x, y=g3, color=k)) +
  geom_point(size = 2, shape = 4, stroke = 1.2, alpha=0.7) +
  labs(
    x = "theta",
    y = "G4 test",
    color = "distance"
  ) +
  theme_minimal()
plot_g_rho + plot_g_delta + plot_g_theta


df_s <- rbind(df_s50[, c(3, 4, 5, 6)], df_s200[, c(3, 4, 5, 6)])
df_s <- rbind(df_s, df_s2000[, c(3, 4, 5, 6)])
plot_s_rho <- ggplot(subset(df_s, type=="rho"), aes(x=x, y=s, color=k)) +
  geom_point(size = 2, shape = 4, stroke = 1.2, alpha=0.7) +
  labs(
    x = "rho",
    y = "S",
    color = "distance"
  ) +
  theme_minimal()
plot_s_delta <- ggplot(subset(df_s, type=="delta"), aes(x=x, y=s, color=k)) +
  geom_point(size = 2, shape = 4, stroke = 1.2, alpha=0.7) +
  labs(
    x = "delta",
    y = "S",
    color = "distance"
  ) +
  theme_minimal()
plot_s_theta <- ggplot(subset(df_s, type=="theta"), aes(x=x, y=s, color=k)) +
  geom_point(size = 2, shape = 4, stroke = 1.2, alpha=0.7) +
  labs(
    x = "theta",
    y = "S",
    color = "distance"
  ) +
  theme_minimal()
plot_s_rho + plot_s_delta + plot_s_theta








