# This R script is for making the figures for the manuscript
# color-blind friendly palette
cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")

custom_theme = function() {
  theme_bw() +
    theme(legend.position = "bottom",
          axis.title = element_text(size = 8))
}


create_rt_estimation_error_line <- function(data) {
  data %>%
    pivot_longer(cols = starts_with("error"), names_to = "method", names_prefix = "error.", values_to = "error") %>%
    ggplot() +
    geom_line(aes(time, error, group = paste(simulation_run, method, sep="_"), color = method), alpha = 0.2) +
    ylab("Predicted Rt - Actual Rt") +
    custom_theme() +
    scale_color_manual(values=cbPalette) +
    scale_fill_manual(values=cbPalette)
}

create_rt_estimation_error_boxplot <- function(data) {
  data %>%
    pivot_longer(., cols = starts_with("rt."), names_to = c("Method", "parameter"),
                 names_pattern = "rt.([[:alpha:]]+)_([[:print:]]+)", values_to = "value") %>%
    pivot_wider(names_from = "parameter", values_from = "value") %>%
    mutate(error = log10(abs(mean_error))) %>%
    ggplot(., aes(Method, error)) +
    geom_violin(aes(fill = Method)) +
    geom_boxplot(width = 0.1, outlier.alpha = 0.5) +
    ylab("Absolute Estimation Error (log10)") +
    custom_theme() +
    #scale_y_log10(labels = scales::label_number(accuracy = 0.00001)) +
    scale_color_manual(values=cbPalette) +
    scale_fill_manual(values=cbPalette)
}

create_rt_estimation_error_by_rt <- function(data) {
  data %>%
    ggplot() +
    geom_point(aes(r_eff, rt.simple_mean_error), size=0.5, alpha=0.3, color = cbPalette[2]) +
    labs(y = "Estimation Error", x = "True Rt") +
    custom_theme() +
    scale_color_manual(values=cbPalette) +
    scale_fill_manual(values=cbPalette)
}

create_estimate_plot <- function(data) {
  data = data %>%
    pivot_longer(., cols = starts_with("rt."), names_to = c("Method", "parameter"),
                 names_pattern = "rt.([[:alpha:]]+)_([[:print:]]+)", values_to = "value") %>%
    pivot_wider(names_from = "parameter", values_from = "value")

  ggplot(data) +
    geom_line(aes(time, r_eff), color="black", size=2) +
    geom_line(aes(time, mean, group = Method, color = Method), size = 1) +
    geom_ribbon(aes(time, ymin = quantile_025, ymax = quantile_975, group = Method, fill = Method), alpha = 0.2) +
    #geom_text(data = unique(data[,"simulation_run"]), x = 150, y = 1.8, aes(label = simulation_run), size = 2) +
    ylab("Rt") +
    custom_theme() +
    scale_color_manual(values=cbPalette) +
    scale_fill_manual(values=cbPalette)
}
