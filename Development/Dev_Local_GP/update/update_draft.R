object <- joint_model_fit_1
object <- joint_model_fit_2

update <- function(object, ...) {
  extra_args <- match.call(expand.dots = FALSE)$...
  extra_args
}

update(object)

?match.call

JMbayes2:::extract_log_sigmas
