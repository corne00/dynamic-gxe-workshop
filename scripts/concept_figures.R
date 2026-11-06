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
  geom_line(colour = grey, linewidth = 0.6) +
  geom_point(colour = blue, size = 0.8) +
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
  geom_col(fill = blue, width = 0.8) +
  geom_col(data = d[i, ], fill = bur, width = 1.8) +
  labs(title = "(c) Daily increments", x = "days after sowing",
       y = "t/ha/day") +
  concept()
save(list(p1, p2, p3), "concept_increments.png", 9, 2.5)

## 2. Solution: the fitted equation integrated over the season -----------------------
P <- readRDS("precomputed/cv_stlsq.rds")
pred <- P[g, v, 1:(tr$nend[g, v] + 1)]
e <- data.frame(das = d$das, rate = c(diff(pred), NA), pred = pred, obs = yy)
q1 <- ggplot(e[-nrow(e), ], aes(das, rate)) +
  geom_col(fill = blue, width = 0.8) +
  labs(title = expression("(a) Predicted growth " * theta(u[i])^T * hat(xi)),
       x = "days after sowing", y = "t/ha/day") +
  concept()
q2 <- ggplot(e, aes(das)) +
  geom_line(aes(y = obs, colour = "APSIM"), linewidth = 0.9) +
  geom_line(aes(y = pred, colour = "cumulative sum"), linewidth = 0.9) +
  scale_colour_manual(values = c(APSIM = grey, `cumulative sum` = blue)) +
  annotate("point", x = e$das[1], y = e$pred[1], colour = bur, size = 2.4) +
  annotate("text", x = e$das[1] + 3, y = e$pred[1] + 1.2, label = "y[0]",
           parse = TRUE, colour = bur, size = 3.4) +
  labs(title = expression("(b) " * y[i + 1] == y[i] + theta(u[i])^T * hat(xi)),
       x = "days after sowing", y = "biomass (t/ha)") +
  concept() +
  theme(legend.position = c(0.28, 0.85), legend.title = element_blank(),
        legend.key.height = unit(3.5, "mm"))
save(list(q1, q2), "concept_solution.png", 8, 2.4)

## 3. Leave-one-season-out validation ----------------------------------------
n <- length(tr$years)
folds <- c(1, 2, 3, 31)
cv <- do.call(rbind, lapply(seq_along(folds), function(k) {
  data.frame(fold = k, season = 1:n,
             role = ifelse(1:n == folds[k], "predicted", "fitted"))
}))
cv$fold <- factor(cv$fold, 4:1, labels = paste("fit", rev(folds)))
pv <- ggplot(cv, aes(season, fold, fill = role)) +
  geom_tile(colour = "white", linewidth = 0.5, height = 0.8) +
  scale_fill_manual(values = c(fitted = "#B8CBDC", predicted = bur)) +
  scale_x_continuous(breaks = c(1, 10, 20, 31),
                     labels = tr$years[c(1, 10, 20, 31)],
                     expand = c(0.01, 0)) +
  labs(x = "season", y = NULL) +
  theme_dgxe() +
  theme(legend.position = "top", legend.title = element_blank(),
        axis.line.y = element_blank(), axis.ticks.y = element_blank())
save(pv, "concept_validation.png", 7, 1.9)

## 4. Regression on the library: dy ~ X xi (genotype g114) ------------------
dd <- make_design(Theta, tr, "g114")
xi <- fit_stlsq(dd$X, dd$dy)
rows <- 1:30
Xs <- scale(dd$X[rows, ])
k <- ncol(Xs)
tiles <- data.frame(row = rep(rows, k), col = rep(seq_len(k), each = length(rows)),
                    v = as.vector(Xs))
dyd <- data.frame(row = rows, v = dd$dy[rows] / sd(dd$dy[rows]))
xid <- data.frame(row = seq_len(k) * length(rows) / k, v = xi[-1])
xid$h <- length(rows) / k
x0 <- 3 + k + 3                                   # position of xi
fillscale <- scale_fill_gradient2(low = bur, mid = "white", high = blue,
                                  limits = c(-3, 3), oob = scales::squish,
                                  na.value = "#F3F4F6", guide = "none")
labs_x <- data.frame(x = seq_len(k) + 2, label = colnames(dd$X))
pd <- ggplot() +
  geom_tile(data = dyd, aes(0, row, fill = v), colour = "white", width = 0.95,
            linewidth = 0.2) +
  geom_tile(data = tiles, aes(col + 2, row, fill = v), colour = "white",
            linewidth = 0.2) +
  geom_tile(data = xid, aes(x0, row, fill = ifelse(v == 0, NA, v / max(abs(v)) * 3)),
            height = xid$h[1] * 0.95, width = 0.95, colour = "white") +
  annotate("text", x = 1, y = 15.5, label = "''%~~%''", parse = TRUE,
           size = 6) +
  annotate("text", x = k + 3.5, y = 15.5, label = "''%*%''", parse = TRUE,
           size = 5) +
  annotate("text", x = c(0, 2 + (k + 1) / 2, x0), y = -1.8, parse = TRUE,
           label = c("Delta*bold(y)", "bold(X)", "bold(xi)"), size = 4.4) +
  annotate("text", x = c(0, 2 + (k + 1) / 2, x0), y = -4, size = 3,
           colour = "#6B7280",
           label = c("increments", "library terms", "coefficients")) +
  geom_text(data = labs_x, aes(x, 31.4, label = label), angle = 90,
            hjust = 1, size = 2.6, colour = "#6B7280") +
  scale_y_reverse() +
  fillscale +
  coord_cartesian(clip = "off", ylim = c(37, -5)) +
  theme_void() +
  theme(plot.margin = margin(4, 4, 4, 4))
save(pd, "concept_design.png", 8, 3.4)
