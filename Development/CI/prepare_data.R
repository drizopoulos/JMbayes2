data("pbc2", package = "JM")
data("pbc2.id", package = "JM")
pbc2 <- pbc2[c("id", "year", "serBilir", "prothrombin", "years", "status",
               "age", "sex", "drug")]

pbc2.id <- pbc2.id[c("id", "years", "status", "age", "sex", "drug")]
# we take as intermediate event prothrombin >= 12
f <- function (x) {
    n <- length(x)
    out <- numeric(n)
    ind <- which(x >= 11.8)
    if (length(ind)) out[seq(ind[1L], n)] <- 1
    out
}
pbc2$IE <- with(pbc2, ave(prothrombin, id, FUN = f))
splt <- split(pbc2, pbc2$id)
f <- function (d) {
    x <- d$year[d$IE == 1]
    d$S <- if (length(x)) min(x) else 0
    d
}
pbc2 <- do.call("rbind", lapply(splt, f))
pbc2$status2 <- as.numeric(pbc2$status != "alive")
row.names(pbc2) <- 1:nrow(pbc2)


####
pbc2_CR <- pbc2[c("id", "years", "status", "status2", "IE", "year", "age",
                  "sex", "drug")]
pbc2_CR <- pbc2_CR[!duplicated(pbc2_CR[c("id", "years", "status", "IE")]), ]
pbc2_CR$start <- pbc2_CR$year
splt <- split(pbc2_CR[c("id", "start", "years")], pbc2_CR$id)
pbc2_CR$stop <- unlist(lapply(splt, function (d) c(d$start[-1], d$years[1])))
f <- function (x) {
    if (length(x) > 1) {x[seq(length(x) - 1)] <- 0; x} else x
}
pbc2_CR$event <- with(pbc2_CR, ave(status2, id, FUN = f))
splt <- split(pbc2_CR, pbc2_CR$id)
f <- function (d) {
    n <- nrow(d)
    d <- d[rep(seq_len(n), 2), ]
    d$CR <- factor(rep(1:2, each = n), labels = c("dead", "transplanted"))
    d$event <- (as.character(d$CR) == as.character(d$status)) * d$event
    d
}
pbc2_CR <- do.call("rbind", lapply(splt, f))
row.names(pbc2_CR) <- 1:nrow(pbc2_CR)

pbc2.idCR <- crisk_setup(pbc2.id, statusVar = "status", censLevel = "alive",
                    nameStrata = "CR")


rm(f, splt)

pbc2.id$status2 <- as.numeric(pbc2.id$status != "alive")

####################################################################
####################################################################


effectPlotData <- function (object, newdata, orig_data, ...) {
    if (inherits(object, "MixMod")) {
        return(GLMMadaptive::effectPlotData(object, newdata, ...))
    }
    form <- formula(object)
    namesVars <- all.vars(form)
    betas <- if (!inherits(object, "lme")) coef(object) else fixef(object)
    V <- if (inherits(object, "geeglm")) object$geese$vbeta else vcov(object)
    orig_data <- orig_data[complete.cases(orig_data[namesVars]), ]
    Terms <- delete.response(terms(form))
    mfX <- model.frame(Terms, data = orig_data)
    Terms_new <- attr(mfX, "terms")
    mfX_new <- model.frame(Terms_new, newdata, xlev = .getXlevels(Terms, mfX))
    X <- model.matrix(Terms_new, mfX_new)
    pred <- c(X %*% betas)
    ses <- sqrt(diag(X %*% V %*% t(X)))
    newdata$pred <- pred
    newdata$low <- pred - 1.96 * ses
    newdata$upp <- pred + 1.96 * ses
    newdata
}

