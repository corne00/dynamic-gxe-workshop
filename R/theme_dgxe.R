## Plot style used throughout the workshop: one palette, one theme.
dgxe <- c(blue = "#014E81", cyan = "#00A6D6", sky = "#59C5F2", green = "#469C2A",
          orange = "#F8961E", red = "#F94144", grey = "#6B7280")
theme_dgxe <- function(base_size = 11) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(colour = "#e5e7eb", linewidth = 0.3),
      axis.title = ggplot2::element_text(colour = "#374151"),
      plot.title = ggplot2::element_text(colour = dgxe[["blue"]], face = "bold", size = base_size + 1),
      strip.text = ggplot2::element_text(colour = dgxe[["blue"]], face = "bold"),
      legend.position = "bottom", legend.title = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(6, 8, 6, 6))
}
