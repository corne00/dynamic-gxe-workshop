## Small helpers for the workshop pages.
## show_fun("a", "b"): print the source (with comments) of functions defined
## in R/functions.R and python/dgxe.py, as synchronised R / Python tabs.
## Comment lines directly above a function belong to it (banners do not).
with_comments <- function(lines, start, end, comment) {
  while (start > 1 && grepl(comment, lines[start - 1]) &&
         !grepl("^(#|##) ----", lines[start - 1])) start <- start - 1
  paste(lines[start:end], collapse = "\n")
}

r_source <- function(names) {
  lines <- readLines("R/functions.R")
  ex  <- parse("R/functions.R", keep.source = TRUE)
  src <- attr(ex, "srcref")
  nm  <- vapply(ex, function(e) {
    if (is.call(e) && identical(e[[1]], as.name("<-"))) {
      as.character(e[[2]])[1]
    } else ""
  }, "")
  i <- match(names, nm)
  stopifnot(!anyNA(i))
  vapply(i, function(k) {
    with_comments(lines, src[[k]][1], src[[k]][3], "^##")
  }, "")
}

py_source <- function(names) {
  lines <- readLines("python/dgxe.py")
  one <- function(f) {
    i <- grep(paste0("^def ", f, "\\("), lines)
    stopifnot(length(i) == 1)
    nxt <- which(seq_along(lines) > i & grepl("^[^ \t]", lines))
    end <- if (length(nxt)) nxt[1] - 1 else length(lines)
    while (!nzchar(lines[end])) end <- end - 1
    with_comments(lines, i, end, "^#")
  }
  vapply(names, one, "")
}

show_fun <- function(..., level = 2) {
  f <- c(...)
  h <- strrep("#", level)
  cat('::: {.panel-tabset group="language"}\n\n', h, ' R\n\n```r\n',
      paste(r_source(f), collapse = "\n\n"), "\n```\n\n",
      h, " Python\n\n```python\n",
      paste(py_source(f), collapse = "\n\n\n"), "\n```\n\n:::\n\n", sep = "")
}

## cached("name.rds", expr): compute once, then reuse (stored in precomputed/)
cached <- function(file, expr) {
  f <- file.path("precomputed", file)
  if (file.exists(f)) return(readRDS(f))
  v <- expr; dir.create("precomputed", showWarnings = FALSE); saveRDS(v, f); v
}
