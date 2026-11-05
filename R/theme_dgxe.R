## Plot style: one blue for the model, grey for context, burgundy for the comparison.
dgxe <- c(blue = "#014E81", cyan = "#00A6D6", sky = "#59C5F2", green = "#469C2A",
          burgundy = "#8B1E3F", red = "#F94144", grey = "#9CA3AF", ink = "#1F2937")
theme_dgxe <- function(base_size = 10.5) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      axis.line = ggplot2::element_line(colour = "#6B7280", linewidth = 0.3),
      axis.ticks = ggplot2::element_line(colour = "#6B7280", linewidth = 0.3),
      axis.text = ggplot2::element_text(colour = "#374151"),
      axis.title = ggplot2::element_text(colour = "#374151"),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(colour = "#1F2937", hjust = 0, face = "plain"),
      legend.position = "none",
      panel.spacing = ggplot2::unit(1.1, "lines"),
      plot.margin = ggplot2::margin(4, 8, 4, 4))
}
