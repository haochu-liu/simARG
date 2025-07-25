recomb_info1 <- function(tree_length, L, rho_s, delta) {
  rho <- rho_s * L
  n_recomb <- rpois(1, rho*tree_length/2)
  probstart <- rep(1, L)
  probstart[1] <- delta
  probstart <- probstart / sum(probstart)
  probstartcum <- cumsum(probstart)
  recomb_pos <- matrix(NA, nrow=n_recomb, ncol=2)

  for (i in 1:n_recomb) {
    recomb_pos[i, 1] <- which(runif(1) < probstartcum)[1]
    recomb_pos[i, 2] <- min(recomb_pos[i, 1] + rgeom(1, 1/delta), L)
  }

  return(recomb_pos)
}


recomb_info2 <- function(tree_length, L, rho_s, delta) {
  rho <- rho_s * L
  recomb_pos <- matrix(NA, nrow=1000, ncol=2)
  recomb_num <- 0
  for (i in 1:L) {
    if (i == 1) {
      R_new <- rpois(1, rho_s*delta*tree_length/2)
      R_old <- 0
    } else {
      survive_index <- which(recomb_pos[1:recomb_num, 2] == (i-1))
      R_new <- rpois(1, rho_s*tree_length/2)
      if (length(survive_index) >= 0) {
        R_old <- rbinom(1, length(survive_index), (1 - 1/delta))
        if (length(survive_index) == 1) {
          remain_index <- survive_index
        } else {
          remain_index <- sample(survive_index, R_old)
        }
      } else {
        R_old <- 0
      }
    }

    if (R_new > 0) {
      recomb_pos[(recomb_num+1):(recomb_num+R_new), ] <- i
    }

    if (R_old > 0) {
      recomb_pos[remain_index, 2] <- i
    }

    recomb_num <- recomb_num + R_new
  }

  return(recomb_pos[1:recomb_num, ])
}


tree_length <- 2
L <- 100
rho_s <- 10 / L
delta <- 20
start_pos <- data.frame(
  method1 = rep(0, 100),
  method2 = rep(0, 100)
)
segment_length1 <- c()
segment_length2 <- c()
for (i in 1:1000) {
  m1 <- recomb_info1(tree_length, L, rho_s, delta)
  m2 <- recomb_info2(tree_length, L, rho_s, delta)

  if (nrow(m1) > 0) {
    t1 <- as.data.frame(table(m1[, 1]))
    t1_index <- as.numeric(as.character(t1[, 1]))
    start_pos$method1[t1_index] <- start_pos$method1[t1_index] + t1[, 2]
    segment_length1 <- c(segment_length1, m1[, 2] - m1[, 1])
  }

  if (nrow(m2) > 0) {
    t2 <- as.data.frame(table(m2[, 1]))
    t2_index <- as.numeric(as.character(t2[, 1]))
    start_pos$method2[t2_index] <- start_pos$method2[t2_index] + t2[, 2]
    segment_length2 <- c(segment_length2, m2[, 2] - m2[, 1])
  }
}


library(ggplot2)

start_pos$y <- 1:100
ggplot() +
  geom_point(data = start_pos,
             aes(x = y, y = method1),
             color = "darkblue", size = 3, shape = 16, alpha = 0.5) +

  geom_point(data = start_pos,
             aes(x = y, y = method2),
             color = "darkred", size = 3, shape = 17, alpha = 0.5) +

  geom_line(data = start_pos,
            aes(x = y, y = method1),
            color = "darkblue", linetype = "solid", linewidth = 1.2, alpha = 0.5) +

  geom_line(data = start_pos,
            aes(x = y, y = method2),
            color = "darkred", linetype = "solid", linewidth = 1.2, alpha = 0.5) +

  labs(
    title = "start pos",
    x = "pos",
    y = "num"
  ) +
  theme_minimal()


segment_df <- data.frame(
  seg = c(segment_length1, segment_length2),
  method = c(rep("m1", length(segment_length1)), rep("m2", length(segment_length2)))
)
ggplot(segment_df, aes(x = seg, fill = method)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 30,
                 color = "black",
                 alpha = 0.3,
                 position = "identity") +
  labs(title = "segment length",
       x = "length",
       y = "Density") +
  scale_fill_manual(values = c("m1"="darkblue", "m2"="darkred")) +
  theme_minimal()
