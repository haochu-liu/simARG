library(ggplot2)
set.seed(100)


height_df <- data.frame(local=rep(NA, 1000),
                        arg=rep(NA, 1000))

for (i in 1:1000) {
  ARG <- FSM_ARG(100L, 0.05, 100L, bacteria = TRUE, delta = 10,
                 node_max = 1000, optimise_recomb = TRUE)
  tree1 <- local_tree(ARG, 1L)
  height_df$local[i] <- tree1$sum_time
  height_df$arg[i] <- ARG$sum_time
  if (i%%100 == 0) {print(paste("Complete", i, "iterations"))}
}


plot_df <- data.frame(y=c(height_df$local, height_df$arg),
                      type=c(rep("Local tree", 1000), rep("ARG", 1000)))
ggplot(plot_df, aes(x=type, y=y)) +
  geom_violin() +
  labs(title = "# of leaf lineages = 100, rho = 5",
       x = "",
       y = "Graph height") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 16),
    axis.text.x = element_text(size = 12),
    axis.title.y = element_text(size = 16),
    axis.text.y = element_text(size = 12),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 12)
  )
