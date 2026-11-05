## Small helpers for the workshop pages.
## show_fun("a", "b"): print the source (with comments) of functions defined in R/functions.R
show_fun <- function(...) {
  ex  <- parse("R/functions.R", keep.source = TRUE)
  src <- attr(ex, "srcref")
  nm  <- vapply(ex, function(e) if (is.call(e) && identical(e[[1]], as.name("<-"))) as.character(e[[2]])[1] else "", "")
  for (f in c(...)) {
    i <- match(f, nm); stopifnot(!is.na(i))
    cat("```r\n", paste(as.character(src[[i]]), collapse = "\n"), "\n```\n\n", sep = "")
  }
}
## cached("name.rds", expr): compute once, then reuse (stored in precomputed/)
cached <- function(file, expr) {
  f <- file.path("precomputed", file)
  if (file.exists(f)) return(readRDS(f))
  v <- expr; dir.create("precomputed", showWarnings = FALSE); saveRDS(v, f); v
}
