## Daily mean temperature of every year for the toy figure on the method
## page (3-day moving average of (Tmin + Tmax) / 2); writes
## data/interactive/weather.json. Run from the repository root.
library(jsonlite)
wx <- read.csv("data/narrabri_weather.csv")
smooth3 <- function(x) stats::filter(x, rep(1 / 3, 3), sides = 2)
by_year <- lapply(split(wx, wx$year), function(w) {
  t <- (w$Tmin + w$Tmax) / 2
  s <- as.numeric(smooth3(t))
  s[is.na(s)] <- t[is.na(s)]
  round(s, 1)
})
write_json(list(years = as.integer(names(by_year)), das0 = min(wx$das),
                temperature = unname(by_year)),
           "data/interactive/weather.json", auto_unbox = TRUE)
