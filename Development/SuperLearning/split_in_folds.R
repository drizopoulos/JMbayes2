library("JMbayes2")
create_folds <- function (data, V = 5, id_var = "id", seed = 123L) {
    if (!exists(".Random.seed", envir = .GlobalEnv))
        runif(1L)
    RNGstate <- get(".Random.seed", envir = .GlobalEnv)
    on.exit(assign(".Random.seed", RNGstate, envir = .GlobalEnv))
    data <- as.data.frame(data)
    ids <- data[[id_var]]
    unq_ids <- unique(ids)
    n <- length(unq_ids)
    splits <- split(seq_len(n), sample(rep(seq_len(V), length.out = n)))
    training <- testing <- vector("list", V)
    for (i in seq_along(training)) {
        ind <- ids %in% unq_ids[splits[[i]]]
        training[[i]] <- data[!ind, ]
        testing[[i]] <- data[ind, ]
    }
    list("training" = training, "testing" = testing)
}

newdata <- create_folds(pbc2)
fit_models <- function (data) {
    lmeFit <- lme(log(serBilir) ~ year * sex, data = data,
                  random = ~ year | id)
    data_id <- data[!duplicated(data$id), ]
    CoxFit <- coxph(Surv(years, status2) ~ sex, data = data_id)
    jmFit1 <- jm(CoxFit, lmeFit, time_var = "year")
    jmFit2 <- jm(CoxFit, lmeFit, time_var = "year",
                 functional_forms = ~ slope(log(serBilir)))
    list(M1 = jmFit1, M2 = jmFit2)
}
object <- lapply(newdata$training, fit_models)


Tstart = 5
Thoriz = 7
##
# if newdata is a list with components 'training' and 'testing',
# Super Learning will be used
if (!is.data.frame(newdata) &&
    all(names(newdata) %in% c("training", "testing"))) {
    CV_data <- newdata
    newdata <- do.call("rbind", CV_data$testing)
    newdata[["fold_"]] <- rep(seq_along(CV_data$testing),
                              sapply(CV_data$testing, nrow))
}
# if Super Learning, object needs to be a list with length the
# number of folds. In each element of the list, we have a list of fitted
# models
obj <- if (inherits(object, "jm")) object else object[[1L]][[1L]]
id_var <- obj$model_info$var_names$idVar
time_var <- obj$model_info$var_names$time_var
Time_var <- obj$model_info$var_names$Time_var
event_var <- obj$model_info$var_names$event_var
if (is.null(newdata[[id_var]])) {
    stop("cannot find the '", id_var, "' variable in newdata.", sep = "")
}
if (is.null(newdata[[time_var]])) {
    stop("cannot find the '", time_var, "' variable in newdata.", sep = "")
}
if (any(sapply(Time_var, function (nmn) is.null(newdata[[nmn]])))) {
    stop("cannot find the '", paste(Time_var, collapse = ", "),
         "' variable(s) in newdata.", sep = "")
}
if (is.null(newdata[[event_var]])) {
    stop("cannot find the '", event_var, "' variable in newdata.", sep = "")
}
newdata <- newdata[newdata[[Time_var]] > Tstart, ]
newdata <- newdata[newdata[[time_var]] <= Tstart, ]
if (!nrow(newdata)) {
    stop("there are no data on subjects who had an observed event time after Tstart ",
         "and longitudinal measurements before Tstart.")
}
newdata[[id_var]] <- newdata[[id_var]][, drop = TRUE]
test1 <- newdata[[Time_var]] < Thoriz & newdata[[event_var]] == 1
if (!any(test1)) {
    stop("it seems that there are no events in the interval [Tstart, Thoriz).")
}
newdata2 <- newdata
newdata2[[Time_var]] <- Tstart
newdata2[[event_var]] <- 0

id <- newdata[[id_var]]
Time <- newdata[[Time_var]]
event <- newdata[[event_var]]
f <- factor(id, levels = unique(id))
Time <- tapply(Time, f, tail, 1L)
event <- tapply(event, f, tail, 1L)
names(Time) <- names(event) <- as.character(unique(id))

# subjects who had the event before Thoriz
ind1 <- Time < Thoriz & event == 1
# subjects who had the event after Thoriz
ind2 <- Time > Thoriz
# subjects who were censored in the interval (Tstart, Thoriz)
ind3 <- Time < Thoriz & event == 0

if (!inherits(object, "jm")) {
    # Super Learning
    V <- length(object) # number of folds
    L <- length(object[[1]]) # number of models
    ns <- sapply(CV_data$testing, nrow) # size of testing datasets
    ids <- tapply(newdata2[[id_var]], newdata2[["fold_"]], unique)
    ns <- sapply(ids, length) # number of test subjects per fold
    predictions <- weights <- vector("list", V)
    for (v in seq_len(V)) {
        temp_p <- temp_w <- vector("list", L)
        for (l in seq_len(L)) {
            preds <- predict(object[[v]][[l]], process = "event", times = Thoriz,
                             newdata = newdata2[newdata2$fold_ == v, ])
            temp_p[[l]] <- preds$pred[preds$times > Tstart]
            if (any(ind3)) {
                nams <- names(ind3[ind3])
                preds2 <- predict(object[[v]][[l]],
                                  newdata = newdata[id %in% nams, ],
                                  process = "event", times = Thoriz)
                weights <- preds2$pred
                f <- factor(preds2$id, levels = unique(preds2$id))
                names(weights) <- f
                temp_w[[l]] <- tapply(weights, f, tail, 1)
            }
        }
        predictions[[v]] <- do.call("cbind", temp_p)
    }
    predictions <- do.call("rbind", predictions)
} else {
    preds <- predict(object, newdata = newdata2, process = "event",
                     times = Thoriz)
    pi_u_t <- preds$pred
    names(pi_u_t) <- preds$id
    # cumulative risk at Thoriz
    pi_u_t <- pi_u_t[preds$times > Tstart]
    if (any(ind3)) {
        nams <- names(ind3[ind3])
        preds2 <- predict(object, newdata = newdata[id %in% nams, ],
                          process = "event", times = Thoriz, ...)
        weights <- preds2$pred
        f <- factor(preds2$id, levels = unique(preds2$id))
        names(weights) <- f
        weights <- tapply(weights, f, tail, 1)
    }
}
coefs <- rep(0.0, L)
varpi <- exp(coefs) / sum(exp(coefs))
pi_u_t <- rowSums(predictions * rep(varpi, each = nrow(predictions)))

loss <- function (x) x * x
events <- sum(loss(1 - pi_u_t[ind1]), na.rm = TRUE)
no_events <- sum(loss(pi_u_t[ind2]), na.rm = TRUE)
censored <- if (any(ind3)) {
    sum(weights * loss(1 - pi_u_t[ind3]) +
            (1 - weights) * loss(pi_u_t[ind3]), na.rm = TRUE)
} else 0.0
nr <- length(Time)
Brier <- (events + no_events + censored) / nr



