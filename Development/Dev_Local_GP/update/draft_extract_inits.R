function (object) 
{
  if (!inherits(object, "JMbayes")) 
    stop("Use only with 'JMbayes' objects.\n")
  init <- object$postMeans
  if (!is.null(init$sigma)) {
    init$tau <- 1/object$postMeans$sigma^2
    init <- init[names(init) != "sigma"]
  }
  else {
    init$tau <- NA
  }
  init$invD <- solve(object$postMeans$D)
  init
}


extract_init_vals <- function(object) {
  if (!inherits(object, 'jm'))
    stop("Use only with 'jm' object.\n")
  #initial_values <- object$statistics$Mean
  
}


