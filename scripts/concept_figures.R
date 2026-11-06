## Conceptual figures for the method and design pages (written to images/).
## Run from the repository root.
library(ggplot2)
library(grid)
source("R/pipeline.R")
source("R/theme_dgxe.R")

concept <- function() {
  theme_dgxe(10) +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          axis.line = element_line(colour = "#6B7280", linewidth = 0.4,
                                   arrow = arrow(length = unit(1.6, "mm"))),
          plot.title = element_text(hjust = 0, size = 10, colour = "#1F2937"))
}
## Panels side by side (base grid), written as a png.
save <- function(plots, name, w, h) {
  if (inherits(plots, "ggplot")) plots <- list(plots)
  png(file.path("images", name), width = w, height = h, units = "in",
      res = 220, bg = "white")
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(1, length(plots))))
  for (i in seq_along(plots)) {
    print(plots[[i]], vp = viewport(layout.pos.row = 1, layout.pos.col = i))
  }
  dev.off()
}
blue <- dgxe[["blue"]]; cyan <- dgxe[["cyan"]]; bur <- dgxe[["burgundy"]]
grey <- "#D1D5DB"

## 1. Increments of an observed trajectory (genotype g114, season 1998) -------
g <- match("g114", tr$genotypes); v <- match(1998, tr$years)
yy <- tr$Y[g, v, 1:(tr$nend[g, v] + 1)]
d <- data.frame(das = seq_along(yy) + tr$das0 - 1, y = yy)
d$dy <- c(diff(yy), NA)
i <- 85                                           # highlighted day
win <- (i - 3):(i + 4)                            # days shown in the zoom
z <- d[win, ]
pad <- 0.06 * diff(range(z$y))
p1 <- ggplot(d, aes(das, y)) +
  geom_line(colour = blue, linewidth = 0.7) +
  annotate("rect", xmin = min(z$das) - 1, xmax = max(z$das) + 1,
           ymin = min(z$y) - pad - 0.3, ymax = max(z$y) + pad + 0.3,
           fill = NA, colour = bur, linewidth = 0.5) +
  labs(title = "(a) Observed biomass", x = "days after sowing",
       y = expression(y[i] ~ "(t/ha)")) +
  concept()
p2 <- ggplot(z, aes(das, y)) +
  geom_line(colour = grey, linewidth = 0.8) +
  geom_point(colour = blue, size = 2) +
  geom_point(data = d[i:(i + 1), ], colour = bur, size = 2.6) +
  annotate("segment", x = d$das[i], xend = d$das[i + 1], y = d$y[i],
           yend = d$y[i], colour = bur, linewidth = 0.5, linetype = "22") +
  annotate("segment", x = d$das[i + 1], xend = d$das[i + 1], y = d$y[i],
           yend = d$y[i + 1], colour = bur, linewidth = 1) +
  annotate("text", x = d$das[i + 1] + 0.25, y = mean(d$y[i:(i + 1)]),
           label = "Delta*y[i]", parse = TRUE, colour = bur, size = 4,
           hjust = 0) +
  annotate("text", x = d$das[i] - 0.2, y = d$y[i] + pad, label = "y[i]",
           parse = TRUE, colour = bur, size = 3.4, hjust = 1) +
  annotate("text", x = d$das[i + 1] - 0.2, y = d$y[i + 1] + pad,
           label = "y[i+1]", parse = TRUE, colour = bur, size = 3.4,
           hjust = 1) +
  labs(title = expression("(b) Zoom: " * Delta * y[i] == y[i + 1] - y[i]),
       x = "days after sowing", y = "t/ha") +
  concept()
p3 <- ggplot(d[-nrow(d), ], aes(das, dy)) +
  geom_point(colour = blue, size = 1) +
  geom_point(data = d[i, ], colour = bur, size = 2.6) +
  annotate("label", x = d$das[i] + 2,
           y = d$dy[i] + 0.1 * diff(range(d$dy, na.rm = TRUE)),
           label = "Delta*y[i]", parse = TRUE, colour = bur, size = 3.6,
           fill = "white", label.size = 0, hjust = 0) +
  labs(title = "(c) Daily increments", x = "days after sowing",
       y = "t/ha/day") +
  concept()
save(list(p1, p2, p3), "concept_increments.png", 9, 2.5)

## 2. Solution: the fitted equation integrated over the season ---------------
P <- readRDS("precomputed/cv_stlsq.rds")
pred <- P[g, v, 1:(tr$nend[g, v] + 1)]
e <- data.frame(das = d$das, rate = c(diff(pred), NA), pred = pred, obs = yy)
e <- e[-nrow(e), ]
tm <- 100                                         # day shown in the figure
shade <- e[e$das <= tm, ]
q1 <- ggplot(e, aes(das, rate)) +
  geom_area(data = shade, fill = "#B8CBDC") +
  geom_line(colour = blue, linewidth = 0.7) +
  geom_vline(xintercept = tm, colour = bur, linewidth = 0.4,
             linetype = "22") +
  annotate("text", x = tm - 6, y = max(e$rate) * 0.35, label = "area",
           colour = "#1F2937", size = 3.4, hjust = 1) +
  labs(title = "(a) Predicted growth rate",
       x = "days after sowing", y = "t/ha/day") +
  concept()
yt <- e$pred[e$das == tm]
q2 <- ggplot(e, aes(das)) +
  geom_line(aes(y = obs), colour = grey, linewidth = 1) +
  geom_line(aes(y = pred), colour = blue, linewidth = 0.9) +
  geom_vline(xintercept = tm, colour = bur, linewidth = 0.4,
             linetype = "22") +
  annotate("point", x = e$das[1], y = e$pred[1], colour = bur, size = 2.4) +
  annotate("text", x = e$das[1] + 3, y = e$pred[1] + 1.2, label = "y[0]",
           parse = TRUE, colour = bur, size = 3.6) +
  annotate("point", x = tm, y = yt, colour = bur, size = 2.4) +
  annotate("text", x = tm - 3, y = yt + 0.9, label = "y(t)", parse = TRUE,
           colour = bur, size = 3.6, hjust = 1) +
  annotate("text", x = max(e$das), y = e$obs[nrow(e)] - 0.9, label = "APSIM",
           colour = "#6B7280", size = 3.4, hjust = 1) +
  annotate("text", x = max(e$das), y = e$pred[nrow(e)] + 0.9, label = "model",
           colour = blue, size = 3.4, hjust = 1) +
  labs(title = expression("(b) Biomass " * y(t) == y[0] + "area"),
       x = "days after sowing", y = "biomass (t/ha)") +
  concept()
save(list(q1, q2), "concept_solution.png", 8, 2.5)

## 3. Leave-one-season-out validation ----------------------------------------
## Seasons 1, 2, 3 and 31; dots stand for the seasons in between.
shown <- c(1, 2, 3, 5)                            # grid positions
cv <- expand.grid(col = shown, row = shown)
cv$role <- factor(ifelse(cv$col == cv$row, "left out and predicted",
                         "used for fitting"),
                  c("used for fitting", "left out and predicted"))
dots <- rbind(
  data.frame(x = 4 + rep(c(-0.2, 0, 0.2), 4), y = rep(shown, each = 3)),
  data.frame(x = rep(shown, each = 3), y = 4 + rep(c(-0.2, 0, 0.2), 4)),
  data.frame(x = 4 + rep(c(-0.2, 0, 0.2), 1), y = 4 + c(-0.2, 0, 0.2)))
dots$y <- dots$y
pv <- ggplot() +
  geom_tile(data = cv, aes(col, row, fill = role), colour = "white",
            linewidth = 1, width = 1, height = 1) +
  geom_point(data = dots, aes(x, y), colour = "#6B7280", size = 0.7) +
  scale_fill_manual(values = c(`used for fitting` = "#B8CBDC",
                               `left out and predicted` = bur)) +
  scale_x_continuous(breaks = shown, labels = c("1", "2", "3", "31"),
                     expand = c(0.02, 0)) +
  scale_y_reverse(breaks = shown, labels = c("fit 1", "fit 2", "fit 3",
                                             "fit 31"),
                  expand = c(0.02, 0)) +
  labs(x = "season", y = NULL) +
  theme_dgxe() +
  theme(legend.position = "right", legend.title = element_blank(),
        axis.line = element_blank(), axis.ticks = element_blank())
save(pv, "concept_validation.png", 6.5, 2.6)

## 4. Regression on the library: dy = X xi (genotype g114, as a table) -------
dd <- make_design(Theta, tr, "g114")
xi <- fit_stlsq(dd$X, dd$dy)
first <- cumsum(c(1, table(factor(dd$season))))[c(1, 2, 31)]  # first row of seasons 1, 2, 31
blocks <- lapply(1:3, function(s) first[s] + 70 + 0:2)   # mid-season days
ypos <- unlist(lapply(1:3, function(s) (s - 1) * 4 + 1:3))
rows <- unlist(blocks)
k <- ncol(dd$X)
xcol <- 2.4 + rep(0:3, each = 3) * 3.35 + rep(0:2, 4) * 1.05
cells <- data.frame(
  x = rep(xcol, each = length(rows)), y = rep(ypos, k),
  label = sub("^-(0\\.0+)$", "\\1", sprintf("%.2f", as.vector(dd$X[rows, ]))))
dycells <- data.frame(x = 0, y = ypos, label = sprintf("%.2f", dd$dy[rows]))
xicells <- data.frame(x = 2.4 + 4 * 3.35 + 0.9, y = (seq_len(k) - 0.5) * 11 / k + 0.5,
                      label = ifelse(xi[-1] == 0, "0", sub("^-(0\\.0+)$", "\\1", sprintf("%.2f", xi[-1]))),
                      zero = xi[-1] == 0)
hi <- 2                                           # highlighted row of block 1
box <- function(x, y, w, h, fill = "#E8EEF4", colour = "white", lw = 0.6) {
  geom_rect(data = data.frame(x = x, y = y), aes(xmin = x - w / 2,
            xmax = x + w / 2, ymin = y - h / 2, ymax = y + h / 2),
            fill = fill, colour = colour, linewidth = lw)
}
heads <- data.frame(x = xcol, label = colnames(dd$X))
groups <- data.frame(x = 2.4 + (0:3) * 3.35 + 1.05, label = c("TT", "DL", "RAW", "RAD"))
dots <- function(x, y) {                          # three vertical dots
  p <- expand.grid(i = c(-0.25, 0, 0.25), k = seq_along(x))
  geom_point(data = data.frame(x = x[p$k], y = y + p$i), aes(x, y),
             size = 0.6, colour = "#6B7280")
}
pd <- ggplot() +
  box(rep(0, 9), ypos, 1.4, 0.92) +
  box(rep(xcol, each = 9), rep(ypos, k), 1.0, 0.92) +
  box(rep(xicells$x[1], k), xicells$y, 1.5, 11 / k * 0.92,
      fill = ifelse(xicells$zero, "#F3F4F6", "#EBD3DA")) +
  geom_rect(data = data.frame(y = ypos[hi]),
            aes(xmin = -0.8, xmax = 2.4 + 3 * 3.35 + 2.55, ymin = y - 0.5,
                ymax = y + 0.5), fill = NA, colour = bur, linewidth = 0.6) +
  geom_text(data = dycells, aes(x, y, label = label), size = 2.5,
            colour = "#1F2937") +
  geom_text(data = data.frame(x = rep(xcol, each = 9), y = rep(ypos, k),
                              label = cells$label),
            aes(x, y, label = label), size = 2.2, colour = "#1F2937") +
  geom_text(data = xicells, aes(x, y, label = label), size = 2.5,
            colour = ifelse(xicells$zero, "#9CA3AF", bur)) +
  dots(c(0, xcol), 4) + dots(c(0, xcol), 8) +
  annotate("text", x = -1.3, y = c(2, 6, 10), label = c("season 1", "season 2",
           "season 31"), angle = 90, size = 2.8, colour = "#6B7280") +
  annotate("text", x = 1.2, y = 6, label = "''%~~%''", parse = TRUE,
           size = 5) +
  annotate("text", x = 2.4 + 4 * 3.35 + 0.1 - 0.15, y = 6,
           label = "''%*%''", parse = TRUE, size = 4.5) +
  geom_text(data = heads, aes(x, -0.1, label = label), parse = TRUE,
            size = 2.8, colour = "#374151") +
  annotate("text", x = c(0, 2.4 + 1.5 * 3.35 + 1.05, xicells$x[1]), y = -1.6,
           parse = TRUE, label = c("Delta*bold(y)", "bold(X)", "bold(xi)"),
           size = 4.4) +
  annotate("text", x = c(0, 2.4 + 1.5 * 3.35 + 1.05, xicells$x[1]), y = -2.7,
           size = 3, colour = "#6B7280",
           label = c("increments", "library terms (standardised)",
                     "coefficients")) +
  annotate("text", x = 2.4 + 1.5 * 3.35 + 1.05, y = 12.3,
           label = "'row i:'~Delta*y[i] == sum(X[ij]*xi[j], j == 1, K)",
           parse = TRUE, size = 3.3, colour = bur) +
  scale_y_reverse() +
  coord_cartesian(clip = "off", ylim = c(12.8, -3.2), xlim = c(-1.6, 17.6)) +
  theme_void() +
  theme(plot.margin = margin(4, 4, 4, 4))
save(pd, "concept_design.png", 7.4, 3.5)

## 5. The data at a glance: weather, one season, one genotype ----------------
v <- match(1998, tr$years)
g0 <- spread_genotypes(tr)
gm <- match(g0[2], tr$genotypes)
traj_long <- function(A, idx) {
  do.call(rbind, lapply(idx, function(k) {
    y <- A(k)
    data.frame(id = k, das = seq_along(y) + tr$das0 - 1, y = y)
  }))
}
season_rows <- traj_long(function(k) {
  y <- tr$Y[k, v, ]; y[!is.na(y)]
}, seq_along(tr$genotypes))
genotype_rows <- traj_long(function(k) {
  y <- tr$Y[gm, k, ]; y[!is.na(y)]
}, seq_along(tr$years))
w98 <- wx[wx$year == 1998, ]
data_theme <- theme_dgxe(10) +
  theme(plot.title = element_text(hjust = 0, size = 10, colour = "#1F2937"))
d1 <- ggplot(w98, aes(das)) +
  geom_ribbon(aes(ymin = Tmin, ymax = Tmax), fill = "#B8CBDC") +
  geom_line(aes(y = Tmax), colour = blue, linewidth = 0.4) +
  geom_line(aes(y = Tmin), colour = blue, linewidth = 0.4) +
  labs(title = "(a) Weather, season 1998", x = "days after sowing",
       y = expression("temperature range (" * degree * "C)")) +
  data_theme
three_rows <- season_rows[season_rows$id %in% match(g0, tr$genotypes), ]
d2 <- ggplot(season_rows, aes(das, y, group = id)) +
  geom_line(colour = grey, linewidth = 0.2) +
  geom_line(data = three_rows, aes(colour = factor(id)), linewidth = 0.8) +
  scale_colour_manual(values = c(bur, cyan, blue)) +
  labs(title = "(b) 199 genotypes, season 1998", x = "days after sowing",
       y = "biomass (t/ha)") +
  data_theme
d3 <- ggplot(genotype_rows, aes(das, y, group = id)) +
  geom_line(colour = grey, linewidth = 0.3) +
  geom_line(data = genotype_rows[genotype_rows$id == v, ], colour = cyan,
            linewidth = 0.8) +
  labs(title = paste0("(c) Genotype ", g0[2], ", 31 seasons"),
       x = "days after sowing", y = "biomass (t/ha)") +
  data_theme
save(list(d1, d2, d3), "concept_data.png", 8, 2.6)

## 6. Covariates constructed from the weather variables -----------------------
raw <- data.frame(name = c("Tmin", "Tmax", "Rain", "RAD", "VP", "DL"),
                  y = 6:1)
cov <- data.frame(name = c("TT", "RAW", "RAD", "DL"), y = c(5.5, 3.6, 2, 1),
                  desc = c("thermal time", "relative available water",
                           "radiation", "day length"))
edges <- rbind(
  data.frame(from = c("Tmin", "Tmax"), to = "TT"),
  data.frame(from = c("Tmin", "Tmax", "Rain", "RAD", "VP"), to = "RAW"),
  data.frame(from = "RAD", to = "RAD"), data.frame(from = "DL", to = "DL"))
xr <- 0; xc <- 7.2; wb <- 1.3; wc <- 5.8
edges$x1 <- xr + wb / 2
edges$y1 <- raw$y[match(edges$from, raw$name)]
edges$x2 <- xc - wc / 2
edges$y2 <- cov$y[match(edges$to, cov$name)]
pc <- ggplot() +
  geom_segment(data = edges, aes(x = x1, y = y1, xend = x2, yend = y2),
               colour = "#9CA3AF", linewidth = 0.4,
               arrow = arrow(length = unit(1.6, "mm"), type = "closed")) +
  geom_rect(data = raw, aes(xmin = xr - wb / 2, xmax = xr + wb / 2,
                            ymin = y - 0.32, ymax = y + 0.32),
            fill = "#E8EEF4", colour = NA) +
  geom_text(data = raw, aes(xr, y, label = name), size = 3.6,
            colour = "#1F2937") +
  geom_rect(data = cov, aes(xmin = xc - wc / 2, xmax = xc + wc / 2,
                            ymin = y - 0.42, ymax = y + 0.42),
            fill = "white", colour = blue, linewidth = 0.5) +
  geom_text(data = cov, aes(xc - wc / 2 + 0.15, y + 0.12, label = name),
            hjust = 0, size = 3.8, fontface = "bold", colour = blue) +
  geom_text(data = cov, aes(xc - wc / 2 + 0.15, y - 0.17, label = desc),
            hjust = 0, size = 3, colour = "#6B7280") +
  annotate("text", x = xc + wc / 2 - 0.15, y = 5.5, hjust = 1, size = 3.2,
           colour = "#1F2937",
           label = "cumulative sum of max((Tmin + Tmax) / 2, 0)") +
  annotate("text", x = xc + wc / 2 - 0.15, y = 3.6 + 0.08, hjust = 1,
           size = 3.3, colour = "#1F2937", parse = TRUE,
           label = "RAW == W / W[max]") +
  annotate("text", x = xc + wc / 2 - 0.15, y = 3.6 - 0.22, hjust = 1,
           label = "soil water W: rain in, demand out", size = 2.6,
           colour = "#6B7280") +
  annotate("text", x = xc + wc / 2 - 0.15, y = c(2, 1), hjust = 1,
           label = "used as given", size = 3, colour = "#6B7280") +
  annotate("text", x = xr, y = 6.9, label = "weather variables (given)",
           size = 3.3, colour = "#6B7280") +
  annotate("text", x = xc, y = 6.9, label = "covariates in the library",
           size = 3.3, colour = "#6B7280") +
  coord_cartesian(clip = "off", xlim = c(-0.8, 10.4), ylim = c(0.4, 7)) +
  theme_void() + theme(plot.margin = margin(4, 4, 4, 4))
save(pc, "concept_covariates.png", 8, 3)
