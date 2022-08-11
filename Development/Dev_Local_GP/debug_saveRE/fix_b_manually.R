stats <- statistics

fix_b <- function (stats) {
    x <- stats$b
    dim(x) <- sapply(dnames_b, length)
    dimnames(x) <- dnames_b
    stats$b <- x
    stats
}

x <- statistics$Mean$b
dim(x) <- sapply(dnames_b, length)

length(x)
length(idL[[2]])

1560 / 5
