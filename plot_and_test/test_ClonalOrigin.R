library(ggplot2)
library(sdprisk)


height_df <- data.frame(s1=rep(NA, 2000),
                        s50=rep(NA, 2000),
                        s80=rep(NA, 2000),
                        clonal=rep(NA, 2000),
                        type=c(rep("1", 1000), rep("2", 1000)))

for (i in 1:1000) {
  set.seed(i)
  Clonaltree1 <- ClonalOrigin2(100L, 10, 100L, 30, TRUE)
  set.seed(i+1000)
  Clonaltree2 <- ClonalOrigin2(100L, 10, 100L, 30, FALSE)
  if (Clonaltree2$n_recomb != nrow(Clonaltree2$recomb_edge)) {
    print(i)
  }
  height_df$s1[i] <- local_height_ClonalOrigin(Clonaltree1, 1)
  height_df$s50[i] <- local_height_ClonalOrigin(Clonaltree1, 50)
  height_df$s80[i] <- local_height_ClonalOrigin(Clonaltree1, 80)
  height_df$clonal[i] <- Clonaltree1$clonal_time
  height_df$s1[i+1000] <- local_height_ClonalOrigin(Clonaltree2, 1)
  height_df$s50[i+1000] <- local_height_ClonalOrigin(Clonaltree2, 50)
  height_df$s80[i+1000] <- local_height_ClonalOrigin(Clonaltree2, 80)
  height_df$clonal[i+1000] <- Clonaltree2$clonal_time
  if (i%%100 == 0) {print(paste("Complete", i, "iterations"))}
}

n <- 100
height_rate <- n:2 * (n-1):1 / 2
height_density <- function(x) {
  dhypoexp(x, rate=height_rate)
}

hist1 <- ggplot(height_df, aes(x = s1, fill = type)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 30,
                 color = "black",
                 alpha = 0.3,
                 position = "identity") +
  stat_function(fun = height_density, # Use your custom function here
                color = "darkblue",
                linewidth = 1.2,
                linetype = "solid") +
  labs(title = "Local tree at site 1",
       x = "Height",
       y = "Density") +
  scale_fill_manual(values = c("1"="blue", "2"="red")) +
  theme_minimal()
hist2 <- ggplot(height_df, aes(x = s50, fill = type)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 30,
                 color = "black",
                 alpha = 0.3,
                 position = "identity") +
  stat_function(fun = height_density, # Use your custom function here
                color = "darkblue",
                linewidth = 1.2,
                linetype = "solid") +
  labs(title = "Local tree at site 50",
       x = "Height",
       y = "Density") +
  scale_fill_manual(values = c("1"="blue", "2"="red")) +
  theme_minimal()
hist3 <- ggplot(height_df, aes(x = s80, fill = type)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 30,
                 color = "black",
                 alpha = 0.3,
                 position = "identity") +
  stat_function(fun = height_density, # Use your custom function here
                color = "darkblue",
                linewidth = 1.2,
                linetype = "solid") +
  labs(title = "Local tree at site 80",
       x = "Height",
       y = "Density") +
  scale_fill_manual(values = c("1"="blue", "2"="red")) +
  theme_minimal()
histClonal <- ggplot(height_df, aes(x = clonal, fill = type)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 30,
                 color = "black",
                 alpha = 0.3,
                 position = "identity") +
  stat_function(fun = height_density, # Use your custom function here
                color = "darkblue",
                linewidth = 1.2,
                linetype = "solid") +
  labs(title = "Clonal tree",
       x = "Height",
       y = "Density") +
  scale_fill_manual(values = c("1"="blue", "2"="red")) +
  theme_minimal()



height_df <- data.frame(s1=rep(NA, 2000),
                        s50=rep(NA, 2000),
                        s80=rep(NA, 2000),
                        clonal=rep(NA, 2000),
                        type=c(rep("1", 1000), rep("2", 1000)))

for (i in 1:1000) {
  set.seed(i)
  ARG1 <- ClonalOrigin_ARG(100L, 10, 100L, 30, optimise_recomb=FALSE)
  if (sum(ARG1$node_clonal==FALSE) != ARG1$n_recomb) {
    print(i)
  }
  set.seed(i+1000)
  ARG2 <- ClonalOrigin_ARG(100L, 10, 100L, 30, optimise_recomb=TRUE)

  tree1 <- local_tree(ARG1, 1L)
  tree50 <- local_tree(ARG1, 50L)
  tree80 <- local_tree(ARG1, 80L)
  clonal_tree_ARG <- clonal_tree_FSM(ARG1)
  height_df$s1[i] <- tree1$sum_time
  height_df$s50[i] <- tree50$sum_time
  height_df$s80[i] <- tree80$sum_time
  height_df$clonal[i] <- clonal_tree_ARG$sum_time
  tree1 <- local_tree(ARG2, 1L)
  tree50 <- local_tree(ARG2, 50L)
  tree80 <- local_tree(ARG2, 80L)
  clonal_tree_ARG <- clonal_tree_FSM(ARG2)
  height_df$s1[i+1000] <- tree1$sum_time
  height_df$s50[i+1000] <- tree50$sum_time
  height_df$s80[i+1000] <- tree80$sum_time
  height_df$clonal[i+1000] <- clonal_tree_ARG$sum_time
  if (i%%100 == 0) {print(paste("Complete", i, "iterations"))}
}

n <- 100
height_rate <- n:2 * (n-1):1 / 2
height_density <- function(x) {
  dhypoexp(x, rate=height_rate)
}

hist1 <- ggplot(height_df, aes(x = s1, fill = type)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 30,
                 color = "black",
                 alpha = 0.3,
                 position = "identity") +
  stat_function(fun = height_density, # Use your custom function here
                color = "darkblue",
                linewidth = 1.2,
                linetype = "solid") +
  labs(title = "Local tree at site 1",
       x = "Height",
       y = "Density") +
  scale_fill_manual(values = c("1"="blue", "2"="red")) +
  theme_minimal()
hist2 <- ggplot(height_df, aes(x = s50, fill = type)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 30,
                 color = "black",
                 alpha = 0.3,
                 position = "identity") +
  stat_function(fun = height_density, # Use your custom function here
                color = "darkblue",
                linewidth = 1.2,
                linetype = "solid") +
  labs(title = "Local tree at site 50",
       x = "Height",
       y = "Density") +
  scale_fill_manual(values = c("1"="blue", "2"="red")) +
  theme_minimal()
hist3 <- ggplot(height_df, aes(x = s80, fill = type)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 30,
                 color = "black",
                 alpha = 0.3,
                 position = "identity") +
  stat_function(fun = height_density, # Use your custom function here
                color = "darkblue",
                linewidth = 1.2,
                linetype = "solid") +
  labs(title = "Local tree at site 80",
       x = "Height",
       y = "Density") +
  scale_fill_manual(values = c("1"="blue", "2"="red")) +
  theme_minimal()
histClonal <- ggplot(height_df, aes(x = clonal, fill = type)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 30,
                 color = "black",
                 alpha = 0.3,
                 position = "identity") +
  stat_function(fun = height_density, # Use your custom function here
                color = "darkblue",
                linewidth = 1.2,
                linetype = "solid") +
  labs(title = "Clonal tree",
       x = "Height",
       y = "Density") +
  scale_fill_manual(values = c("1"="blue", "2"="red")) +
  theme_minimal()



height_df <- data.frame(n_recomb=rep(NA, 2000),
                        ratio=rep(NA, 2000),
                        clonal=rep(NA, 2000),
                        type=c(rep("1", 1000), rep("2", 1000)))

for (i in 1:1000) {
  set.seed(i)
  ARG1 <- ClonalOrigin_ARG(100L, 10, 100L, 30, optimise_recomb=FALSE)

  set.seed(i+1000)
  Clonaltree2 <- ClonalOrigin2(100L, 10, 100L, 30, FALSE)

  clonal_tree_ARG <- clonal_tree_FSM(ARG1)
  height_df$n_recomb[i] <- ARG1$n_recomb
  height_df$ratio[i] <- ARG1$n_recomb / sum(clonal_tree_ARG$edge[, 3])
  height_df$clonal[i] <- clonal_tree_ARG$sum_time

  height_df$n_recomb[i+1000] <- Clonaltree2$n_recomb
  height_df$ratio[i+1000] <- Clonaltree2$n_recomb / sum(Clonaltree2$clonal_edge[, 3])
  height_df$clonal[i+1000] <- Clonaltree2$clonal_time
  if (i%%100 == 0) {print(paste("Complete", i, "iterations"))}
}

n <- 100
height_rate <- n:2 * (n-1):1 / 2
height_density <- function(x) {
  dhypoexp(x, rate=height_rate)
}

hist1 <- ggplot(height_df, aes(x = n_recomb, fill = type)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 30,
                 color = "black",
                 alpha = 0.3,
                 position = "identity") +
  scale_fill_manual(values = c("1"="blue", "2"="red")) +
  theme_minimal()

hist2 <- ggplot(height_df, aes(x = ratio, fill = type)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 30,
                 color = "black",
                 alpha = 0.3,
                 position = "identity") +
  geom_vline(aes(xintercept = 10/2), color = "darkblue",
             linetype = "dashed", linewidth = 1, show.legend = TRUE) +
  scale_fill_manual(values = c("1"="blue", "2"="red")) +
  theme_minimal()

histClonal <- ggplot(height_df, aes(x = clonal, fill = type)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 30,
                 color = "black",
                 alpha = 0.3,
                 position = "identity") +
  stat_function(fun = height_density, # Use your custom function here
                color = "darkblue",
                linewidth = 1.2,
                linetype = "solid") +
  labs(title = "Clonal tree",
       x = "Height",
       y = "Density") +
  scale_fill_manual(values = c("1"="blue", "2"="red")) +
  theme_minimal()


