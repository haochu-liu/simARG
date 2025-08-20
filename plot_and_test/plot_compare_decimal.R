library(ggplot2)
library(patchwork)
library(sdprisk)
library(latex2exp)


rho_values <-  c(5, 10, 20, 30)
func_vec <- c("FSM_ARG(optimise=T)", "FSM_ARG.decimal(optimise=T)",
              "simbac_ARG", "simbac_ARG.decimal")
time_df <- expand.grid(
  t = NA,
  rho = rho_values / 1e5,
  func = func_vec
)

set.seed(100)
for (i in 1:nrow(time_df)) {
  n <- 20L
  rho <- time_df$rho[i]
  delta <- 1000
  func_str <- time_df$func[i]
  time_vec <- rep(NA, 5)
  if (func_str == "FSM_ARG(optimise=T)") {
    for (j in 1:5) {
      time_result <- system.time(
        FSM_ARG(n, rho, 1e5L, bacteria = TRUE, delta = delta,
                node_max = 1000, optimise_recomb = TRUE)
      )
      time_vec[j] <- time_result["elapsed"]
      print(j)
    }
  } else if (func_str == "FSM_ARG.decimal(optimise=T)") {
    for (j in 1:5) {
      time_result <- system.time(
        FSM_ARG.decimal(n, rho, 1e5L, bacteria = TRUE, delta = delta,
                        node_max = 1000, optimise_recomb = TRUE)
      )
      time_vec[j] <- time_result["elapsed"]
      print(j)
    }
  } else if (func_str == "simbac_ARG") {
    for (j in 1:5) {
      time_result <- system.time(
        simbac_ARG(n, rho, 1e5L, delta, node_max = 1000)
      )
      time_vec[j] <- time_result["elapsed"]
      print(j)
    }
  } else if (func_str == "simbac_ARG.decimal") {
    for (j in 1:5) {
      time_result <- system.time(
        simbac_ARG.decimal(n, rho, 1e5L, delta, node_max = 1000)
      )
      time_vec[j] <- time_result["elapsed"]
      print(j)
    }
  }
  print(time_vec)
  time_df$t[i] <- median(time_vec)
  print(paste("Complete", i, "iterations"))
}

time_plot <- time_df
time_plot$rho_per_site <- time_df$rho
time_plot$rho <- time_df$rho * 1e5

# save(time_plot, file = "decimal_time.rda")

ggplot(time_plot, aes(x=rho, y=t, color=func)) +
  geom_line(, linewidth = 1.2, alpha=0.7) +
  geom_point(size = 2, shape = 21, fill = "white", stroke = 1.2, alpha=0.7) +
  scale_color_manual(values=c("FSM_ARG(optimise=T)"="darkblue",
                              "FSM_ARG.decimal(optimise=T)"="darkgreen",
                              "simbac_ARG"="orange",
                              "simbac_ARG.decimal"="yellow")) +
  labs(
    title = "Running time for different ARG simulation methods (n = 1e5)",
    x = "rho",
    y = "time",
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

time_plot <- time_plot[c(1:4, 9:12), ]
time_plot$func <- as.character(time_plot$func)
time_plot$func[1:4] <- "Rejection"
time_plot$func[5:8] <- "SimBac"
ggplot(time_plot, aes(x=rho, y=t, color=func)) +
  geom_line(, linewidth = 1.2, alpha=0.7) +
  geom_point(size = 2, shape = 21, fill = "white", stroke = 1.2, alpha=0.7) +
  scale_color_manual(values=c("Rejection"="darkblue",
                              "SimBac"="darkgreen")) +
  labs(
    title = "Running time for two optimisation methods",
    x = TeX("$\\rho$"),
    y = "Time (second)",
    color = "Optimisation"
  ) +
  # scale_y_continuous(trans='log10') +
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
