library(ggplot2)
library(sdprisk)
set.seed(100)


height_df <- data.frame(s1=rep(NA, 2000),
                        s50=rep(NA, 2000),
                        s80=rep(NA, 2000),
                        arg=rep(NA, 2000),
                        clonal=rep(NA, 2000),
                        type=c(rep("ARG-based", 1000), rep("Tree-based", 1000)))

for (i in 1:1000) {
  ARG <- ClonalOrigin_ARG_based(100L, 10, 100L, 30, optimise_recomb=TRUE)
  Clonaltree <- ClonalOrigin_tree_based(100L, 10, 100L, 30)
  tree1 <- local_tree(ARG, 1L)
  tree50 <- local_tree(ARG, 50L)
  tree80 <- local_tree(ARG, 80L)
  clonal_tree_ARG <- clonal_tree_FSM(ARG)
  height_df$s1[i] <- tree1$sum_time
  height_df$s50[i] <- tree50$sum_time
  height_df$s80[i] <- tree80$sum_time
  height_df$clonal[i] <- clonal_tree_ARG$sum_time
  height_df$arg[i] <- ARG$sum_time
  height_df$s1[i+1000] <- local_height_ClonalOrigin(Clonaltree, 1)
  height_df$s50[i+1000] <- local_height_ClonalOrigin(Clonaltree, 50)
  height_df$s80[i+1000] <- local_height_ClonalOrigin(Clonaltree, 80)
  height_df$clonal[i+1000] <- Clonaltree$clonal_time
  height_df$arg[i+1000] <- Clonaltree$sum_time
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
  scale_fill_manual(values = c("ARG-based"="blue", "Tree-based"="red")) +
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
  scale_fill_manual(values = c("ARG-based"="blue", "Tree-based"="red")) +
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
  scale_fill_manual(values = c("ARG-based"="blue", "Tree-based"="red")) +
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
  scale_fill_manual(values = c("ARG-based"="blue", "Tree-based"="red")) +
  theme_minimal()


local_height1 <- data.frame(height=c(height_df$s1[1:1000],
                                     height_df$s50[1:1000],
                                     height_df$s80[1:1000]),
                            site=c(rep("1", 1000),
                                   rep("50", 1000),
                                   rep("80", 1000)))

histCompare1 <- ggplot(local_height1, aes(x = height, fill = site)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 30,
                 color = "black",
                 alpha = 0.3,
                 position = "identity") +
  stat_function(fun = height_density, # Use your custom function here
                color = "darkblue",
                linewidth = 1.2,
                linetype = "solid") +
  labs(title = "Local tree height values from ARG-based method",
       x = "Height",
       y = "Density") +
  theme_minimal()

local_height2 <- data.frame(height=c(height_df$s1[1001:2000],
                                     height_df$s50[1001:2000],
                                     height_df$s80[1001:2000]),
                            site=c(rep("1", 1000),
                                   rep("50", 1000),
                                   rep("80", 1000)))

histCompare2 <- ggplot(local_height2, aes(x = height, fill = site)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 30,
                 color = "black",
                 alpha = 0.3,
                 position = "identity") +
  stat_function(fun = height_density, # Use your custom function here
                color = "darkblue",
                linewidth = 1.2,
                linetype = "solid") +
  labs(title = "Local tree height values from tree-based method",
       x = "Height",
       y = "Density") +
  theme_minimal()
