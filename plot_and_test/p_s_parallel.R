# Load necessary packages
library(foreach)
library(doParallel)


set.seed(100)
tree <- clonal_genealogy(15L)

ClonalOrigin_pair_seq <- ClonalOrigin_pair_seq
FSM_mutation <- FSM_mutation
LD_r <- LD_r
G3_test <- G3_test

p_s_parallel <- function(theta,
                         tree, ClonalOrigin_pair_seq, FSM_mutation, LD_r, G3_test) {
  # Set up a parallel backend with 5 cores
  # 'makeCluster' creates a cluster of R processes on the local machine
  cl <- makeCluster(5)
  registerDoParallel(cl)

  rho_site <- theta[1]
  delta <- theta[2]
  theta_site <- theta[3]
  s_vec <- rep(NA, 7)



  # Function to run a single iteration for a given loop
  run_simulation <- function(rho_site, delta, theta_site, l_val) {
    ARG <- ClonalOrigin_pair_seq(tree, rho_site, 1e6L, delta, l_val)
    ARG_mutated <- FSM_mutation(ARG, theta_site, binary=TRUE)
    mat <- ARG_mutated$node_site[1:15, ]

    return(list(
      LD_r = LD_r(mat),
      G3_test = G3_test(mat),
      s = any(as.logical(mat[, 1]))
    ))
  }

  # Define loop lengths
  l_values <- c(50L, 200L, 2000L)

  # Use foreach to parallelize the three groups of simulations
  results <- foreach(
    l = l_values,
    .combine = 'list',
    .multicombine = TRUE,
    .packages = c("foreach") # Add any packages needed inside the loop
  ) %dopar% {

    # Run the 2000 simulations for each l_value
    l_results <- foreach(i = 1:2000, .combine = 'rbind') %do% {
      res <- run_simulation(rho_site, delta, theta_site, l)
      c(res$LD_r, res$G3_test, res$s)
    }

    # Convert to a data frame for easier manipulation
    l_results <- as.data.frame(l_results)
    colnames(l_results) <- c("v_r", "v_g3", "v_s")
    l_results
  }

  # Shut down the cluster
  stopCluster(cl)

  # Process the results
  # Note that this part is different since the results are structured differently
  s_vec[1] <- mean(results[[1]]$v_r)
  s_vec[4] <- mean(results[[1]]$v_g3)

  s_vec[2] <- mean(results[[2]]$v_r)
  s_vec[5] <- mean(results[[2]]$v_g3)

  s_vec[3] <- mean(results[[3]]$v_r)
  s_vec[6] <- mean(results[[3]]$v_g3)

  # s_vec[7] is a bit trickier, as v_s is a single vector of length 6000
  # We need to collect all v_s values from all results
  all_s_values <- c(results[[1]]$v_s, results[[2]]$v_s, results[[3]]$v_s)
  s_vec[7] <- mean(all_s_values)

  return(s_vec)
}

time_vec <- rep(NA, 10)
for (i in 1:10) {
  time_result <- system.time(
    s_out <- p_s_parallel(theta_0,
                          tree, ClonalOrigin_pair_seq, FSM_mutation, LD_r, G3_test)
  )
  time_vec[i] <- time_result["elapsed"]
}

time_vec2 <- rep(NA, 10)
for (i in 1:10) {
  time_result <- system.time(
    s_out <- p_s(theta_0)
  )
  time_vec2[i] <- time_result["elapsed"]
  print(i)
}
