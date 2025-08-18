library(ggplot2)
library(patchwork)
library(latex2exp)
library(ape)


# load df_s50, df_s200, df_s2000
df_s50$k <- "50"
df_s200$k <- "200"
df_s2000$k <- "2000"

df_r <- rbind(df_s50[, c(1, 4, 5, 6)], df_s200[, c(1, 4, 5, 6)])
df_r <- rbind(df_r, df_s2000[, c(1, 4, 5, 6)])
plot_r_rho <- ggplot(subset(df_r, type=="rho"), aes(x=x, y=r, color=k)) +
  geom_point(size = 2, shape = 4, stroke = 1.2, alpha=0.7) +
  labs(
    x = TeX("$\\rho$"),
    y = TeX("Mean $r^2$"),
    color = "distance"
  ) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.title.x = element_blank())
plot_r_delta <- ggplot(subset(df_r, type=="delta"), aes(x=x, y=r, color=k)) +
  geom_point(size = 2, shape = 4, stroke = 1.2, alpha=0.7) +
  labs(
    x = TeX("$\\delta$"),
    y = TeX("Mean $r^2$"),
    color = "distance"
  ) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
plot_r_theta <- ggplot(subset(df_r, type=="theta"), aes(x=x, y=r, color=k)) +
  geom_point(size = 2, shape = 4, stroke = 1.2, alpha=0.7) +
  labs(
    x = TeX("$\\theta$"),
    y = TeX("Mean $r^2$"),
    color = "distance"
  ) +
  theme_minimal() +
  theme(legend.position = "right",
        axis.title.x = element_blank(),
        axis.title.y = element_blank())


df_g <- rbind(df_s50[, c(2, 4, 5, 6)], df_s200[, c(2, 4, 5, 6)])
df_g <- rbind(df_g, df_s2000[, c(2, 4, 5, 6)])
plot_g_rho <- ggplot(subset(df_g, type=="rho"), aes(x=x, y=g3, color=k)) +
  geom_point(size = 2, shape = 4, stroke = 1.2, alpha=0.7) +
  labs(
    x = TeX("$\\rho$"),
    y = "Mean G4",
    color = "distance"
  ) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.title.x = element_blank())
plot_g_delta <- ggplot(subset(df_g, type=="delta"), aes(x=x, y=g3, color=k)) +
  geom_point(size = 2, shape = 4, stroke = 1.2, alpha=0.7) +
  labs(
    x = TeX("$\\delta$"),
    y = "Mean G4",
    color = "distance"
  ) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
plot_g_theta <- ggplot(subset(df_g, type=="theta"), aes(x=x, y=g3, color=k)) +
  geom_point(size = 2, shape = 4, stroke = 1.2, alpha=0.7) +
  labs(
    x = TeX("$\\theta$"),
    y = "Mean G4",
    color = "distance"
  ) +
  theme_minimal() +
  theme(legend.position = "right",
        axis.title.x = element_blank(),
        axis.title.y = element_blank())


df_s <- rbind(df_s50[, c(3, 4, 5, 6)], df_s200[, c(3, 4, 5, 6)])
df_s <- rbind(df_s, df_s2000[, c(3, 4, 5, 6)])
plot_s_rho <- ggplot(subset(df_s, type=="rho"), aes(x=x, y=s, color=k)) +
  geom_point(size = 2, shape = 4, stroke = 1.2, alpha=0.7) +
  labs(
    x = TeX("$\\rho$"),
    y = "Mean S",
    color = "distance"
  ) +
  theme_minimal() +
  theme(legend.position = "none")
plot_s_delta <- ggplot(subset(df_s, type=="delta"), aes(x=x, y=s, color=k)) +
  geom_point(size = 2, shape = 4, stroke = 1.2, alpha=0.7) +
  labs(
    x = TeX("$\\delta$"),
    y = "Mean S",
    color = "distance"
  ) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.title.y = element_blank())
plot_s_theta <- ggplot(subset(df_s, type=="theta"), aes(x=x, y=s, color=k)) +
  geom_point(size = 2, shape = 4, stroke = 1.2, alpha=0.7) +
  labs(
    x = TeX("$\\theta$"),
    y = "Mean S",
    color = "distance"
  ) +
  theme_minimal() +
  theme(legend.position = "right",
        axis.title.y = element_blank())


combined_hist <- (plot_r_rho | plot_r_delta | plot_r_theta) /
  (plot_g_rho | plot_g_delta | plot_g_theta) /
  (plot_s_rho | plot_s_delta | plot_s_theta)
combined_hist <- combined_hist +
  plot_annotation(title = "Relationship between model parameters, site distances and the summary statistics") +
  plot_layout(guides = "collect")
print(combined_hist)

set.seed(100)
arg <- FSM_ARG(15L, 0, 10L, bacteria = TRUE, delta = 1)
localtree <- local_tree(arg, 1L)
ape_tree <- localtree_to_phylo(localtree)
plot(ape_tree)
