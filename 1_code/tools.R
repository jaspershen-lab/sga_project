library(ggplot2)
trimester_color =
  c(
    "2nd trimester" = ggsci::pal_jama()(n = 10)[1],
    "3rd trimester" = ggsci::pal_jama()(n = 10)[2]
  )

sga_color =
  c(
    "AGA" = ggsci::pal_aaas()(n=10)[1],
    "SGA" = ggsci::pal_aaas()(n=10)[2]
  )

base_theme =
  theme_bw() +
  theme(
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    # legend.position = c(0, 1),
    # legend.justification = c(0, 1),
    legend.background = element_blank(),
    strip.background = element_rect(fill = "#0099B47F"),
    strip.text = element_text(color = "white", size = 13),
    panel.grid.minor = element_blank()
  )

volcano_plot <- function(fc,
                         p_value,
                         p.cutoff = 0.05,
                         fc.cutoff = 2,
                         theme = c("light", "dark"),
                         alpha = 0.7,
                         text = FALSE,
                         variable_id = NULL,
                         point.size = 1) {
  theme <- match.arg(theme)
  temp_data <- data.frame(
    fc = log(fc, 2),
    p_value = -log(p_value, 10),
    variable_id,
    stringsAsFactors = FALSE
  )
  
  temp_data <-
    temp_data %>%
    dplyr::mutate(
      class = case_when(
        p_value > -log(p.cutoff, 10) & fc > log(fc.cutoff, 2) ~ "Increase",
        p_value > -log(p.cutoff, 10) &
          fc < log(1 / fc.cutoff, 2) ~ "Decrease",
        TRUE ~ "No"
      )
    )
  
  plot <-
    temp_data %>%
    ggplot(aes(fc, p_value)) +
    geom_vline(xintercept = 0,
               color = "black",
               linetype = 2) +
    geom_hline(
      yintercept = -log(p.cutoff, 10),
      color = "#FB8072",
      linetype = 2
    ) +
    geom_point(
      shape = 16,
      aes(color = class,
          size = p_value),
      show.legend = FALSE,
      alpha = alpha,
    ) +
    scale_color_manual(
      values = c(
        "Increase" = ggsci::pal_aaas()(10)[2],
        "Decrease" = ggsci::pal_aaas()(10)[1],
        "No" = "black"
      )
    ) +
    theme_bw() +
    theme(
      axis.title = element_text(size = 13),
      axis.text = element_text(size = 12),
      legend.title = element_text(size = 13),
      legend.text = element_text(size = 12),
      legend.position = c(0, 1),
      legend.justification = c(0, 1),
      panel.background = element_rect(fill = "transparent", color = NA),
      plot.background = element_rect(fill = "transparent", color = NA),
      legend.background = element_rect(fill = "transparent", color = NA),
      strip.background = element_rect(fill = "#0099B47F"),
      strip.text = element_text(color = "white", size = 13),
      panel.grid.minor = element_blank()
    ) +
    xlab(label = expression(log[2](Fold ~ change))) +
    ylab(label = expression(log[10](FDR ~ adjusted ~ P ~ value)))

  
  if (text) {
    plot <-
      plot +
      ggrepel::geom_text_repel(mapping = aes(fc, p_value, 
                                             label = variable_id,
                                             color = class),
                               data = temp_data[temp_data$class != "No", ]) +
      scale_color_manual(
        values = c(
          "Increase" = ggsci::pal_aaas()(10)[2],
          "Decrease" = ggsci::pal_aaas()(10)[1],
          "No" = "black"
        )
      )
  }
  
  plot
  
}