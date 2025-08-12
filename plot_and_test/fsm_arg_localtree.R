library(ggplot2)
library(sdprisk)
library(patchwork)
set.seed(100)


height_t_df <- data.frame(s1=rep(NA, 1000),
                          s50=rep(NA, 1000),
                          s80=rep(NA, 1000))

for (i in 1:1000) {
  ARG <- FSM_ARG(100L, 0.05, 100L, bacteria = TRUE, delta = 10,
                 node_max = 1000, optimise_recomb = TRUE)
  tree1 <- local_tree(ARG, 1L)
  tree50 <- local_tree(ARG, 50L)
  tree80 <- local_tree(ARG, 80L)
  height_t_df$s1[i] <- tree1$sum_time
  height_t_df$s50[i] <- tree50$sum_time
  height_t_df$s80[i] <- tree80$sum_time
  if (i%%100 == 0) {print(paste("Complete", i, "iterations"))}
}

n <- 100
height_rate <- n:2 * (n-1):1 / 2
# x <- seq(0, 10, length.out = 500)
# height_density <- dhypoexp(x, rate=height_rate)
height_density <- function(x) {
  dhypoexp(x, rate=height_rate)
}

hist1 <- ggplot(height_t_df, aes(x = s1)) +
  geom_histogram(aes(y = after_stat(density)),
                 binwidth = 0.5,
                 fill = "lightgreen",
                 color = "black",
                 alpha = 0.7) +
  stat_function(fun = height_density, # Use your custom function here
                color = "darkblue",
                linewidth = 1.2,
                linetype = "solid") +
  labs(title = "Locat tree at site 1",
       x = "Height",
       y = "Density") +
  theme_minimal()
hist2 <- ggplot(height_t_df, aes(x = s50)) +
  geom_histogram(aes(y = after_stat(density)),
                 binwidth = 0.5,
                 fill = "lightgreen",
                 color = "black",
                 alpha = 0.7) +
  stat_function(fun = height_density, # Use your custom function here
                color = "darkblue",
                linewidth = 1.2,
                linetype = "solid") +
  labs(title = "Locat tree at site 50",
       x = "Height",
       y = "Density") +
  theme_minimal()
hist3 <- ggplot(height_t_df, aes(x = s80)) +
  geom_histogram(aes(y = after_stat(density)),
                 binwidth = 0.5,
                 fill = "lightgreen",
                 color = "black",
                 alpha = 0.7) +
  stat_function(fun = height_density, # Use your custom function here
                color = "darkblue",
                linewidth = 1.2,
                linetype = "solid") +
  labs(title = "Locat tree at site 80",
       x = "Height",
       y = "Density") +
  theme_minimal()
combined_hist <- hist1 + hist2 + hist3
combined_hist <- combined_hist +
  plot_annotation(
    title = "Local tree at three different positions")
print(combined_hist)
