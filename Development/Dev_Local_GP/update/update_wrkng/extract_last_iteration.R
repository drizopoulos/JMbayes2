extract_last_iteration <- function(x) {
  namsX <- names(x)
  last_iter_b <- lapply(x[namsX %in% 'b'], function(x) lapply(x, function(x) x[, , dim(x)[length(dim(x))], drop = FALSE]))
  x <- x[!namsX %in% 'b']
  last_iter <- lapply(x, function(x) lapply(x, function(x) x[nrow(x), , drop = TRUE]))
  last_iter <- c(last_iter, last_iter_b)
  last_iter
}
