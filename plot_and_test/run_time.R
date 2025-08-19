library(ggplot2)
library(patchwork)
library(sdprisk)


rho_values <- seq(5, 35, by = 5)
opt <- c("Basic",
         "Rejection",
         "SimBac")
time_df <- expand.grid(
  t = NA,
  rho_site = rho_values / 1e5,
  opt = opt
)

set.seed(100)
for (i in 1:nrow(time_df)) {
  n <- 20L
  rho_site <- time_df$rho_site[i]
  delta <- 30
  opt_str <- time_df$opt[i]
  time_vec <- rep(NA, 10)
  if (opt_str == "Rejection") {
    for (j in 1:10) {
      time_result <- system.time(
        sim_FSM_ARG(n, rho_site, 1e5L, bacteria = TRUE, delta = delta,
                    node_max = 1000, optimise_recomb = TRUE)
      )
      time_vec[j] <- time_result["elapsed"]
    }
  } else if (opt_str == "Basic") {
    if (rho_site * 100 > 10) {next}
    for (j in 1:10) {
      time_result <- system.time(
        sim_FSM_ARG(n, rho_site, 1e5L, bacteria = TRUE, delta = delta,
                    node_max = 1000, optimise_recomb = FALSE)
      )
      time_vec[j] <- time_result["elapsed"]
    }
  } else if (opt_str == "SimBac") {
    for (j in 1:10) {
      time_result <- system.time(
        simbac_ARG(n, rho_site, 1e5L, delta, node_max = 1000)
      )
      time_vec[j] <- time_result["elapsed"]
    }
  }
  time_df$t[i] <- median(time_vec)
  print(paste("Complete", i, "iterations"))
}


time_plot <- time_df
time_plot$n <- paste0("n = ", time_plot$n)
time_plot$delta <- paste0("delta = ", time_plot$delta)

ggplot(time_plot, aes(x=rho, y=t, color=func)) +
  geom_line(aes(group = interaction(n, delta, func)), linewidth = 1.2, alpha=0.7) +
  geom_point(size = 2, shape = 21, fill = "white", stroke = 1.2, alpha=0.7) +
  scale_color_manual(values=c("sim_FSM_ARG(optimise=T)"="darkblue",
                              "sim_FSM_ARG(optimise=T, clonal=T)"="darkgreen",
                              "sim_FSM_ARG(optimise=F)"="darkred",
                              "simbac_ARG(optimise_site=F)"="orange",
                              "simbac_ARG(optimise_site=T)"="yellow")) +
  facet_grid(n ~ delta) +
  labs(
    title = "Running time for different ARG simulation methods",
    x = "rho",
    y = "log time",
    color = "function"
  ) +
  scale_y_continuous(trans='log10') +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 16),
    axis.text.x = element_text(size = 12),
    axis.title.y = element_text(size = 16),
    axis.text.y = element_text(size = 12),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 12),
    legend.position = "bottom"
  )
