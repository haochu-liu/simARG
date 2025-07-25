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
      R_old <- rbinom(length(survive_index), 1, (1 - 1/delta))
      if (length(R_old) == 0) {
        R_old <- 0
        remain_index <- numeric(0)
      } else {
        remain_index <- sample(survive_index, R_old)
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
}
