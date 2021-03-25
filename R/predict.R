predict.jm <- function (object, newdata = NULL, newdata2 = NULL,
              process = c("longitudinal", "event"),
              type_pred = c("response", "link"),
              type = c("subject_specific", "mean_subject"),
              level = 0.95, return_newdata = FALSE,
              n_samples = 200L, n_mcmc = 25L, cores = NULL, seed = 123L, ...) {
    process <- match.arg(process)
    type_pred <- match.arg(type_pred)
    type <- match.arg(type)
    if (is.null(cores)) {
        n <- length(unique(newdata[[object$model_info$var_names$idVar]]))
        cores <- if (n > 20) 4L else 1L
    }
    components_newdata <- get_components_newdata(object, newdata, n_samples,
                                                 n_mcmc, cores, seed)
    if (process == "longitudinal") {
        predict_Long(object, components_newdata, newdata, newdata2, type,
                     type_pred, level, return_newdata)
    } else {
        predict_Event(object, components_newdata, newdata, level,
                      return_newdata)
    }
}
