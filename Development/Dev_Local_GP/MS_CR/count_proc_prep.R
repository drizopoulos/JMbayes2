ms_setup <- function (data, timevars, statusvars, transitionmat, id, covs = NULL) {
  # setup times matrix with NAs
  # First row is NA as this is starting state 
  timesmat <- matrix(NA, nrow(data), length(timevars))
  timecols_data <- match(timevars[!is.na(timevars)], names(data))
  timesmat[, -which(is.na(timevars))] <- as.matrix(data[, timecols_data]) 
  # setup status matrix with NAs
  # First row is NA as this is starting state 
  statusmat <- matrix(NA, nrow(data), length(statusvars))
  statuscols_data <- match(statusvars[!is.na(statusvars)], names(data))
  statusmat[, -which(is.na(statusvars))] <- as.matrix(data[, statuscols_data]) 
  # ensure convert to matrices
  timesmat <- as.matrix(timesmat)
  statusmat <- as.matrix(statusmat)
  # check dimesnions are the same
  if (any(dim(timesmat) != dim(statusmat))) 
    stop("Dimensions of \"time\" and \"status\" data should be equal")
  # components
  # number of unique subjects
  n_subj <- nrow(timesmat)
  # number of states
  n_states <- dim(transitionmat)[1]
  # set start state to 1 and start time to 0 for all subjects
  # ATTENTION: this needs to be adjusted to more flexible to allow subjects starting at different states
  # this could be achieved by a requesting a separate argument (vector with starting state)
  starting_state <- rep(1, n_subj)
  starting_time <- rep(0, n_subj)
  idnam <- id
  id <- data[[id]]
  order_id <- order(id)
  out <- ms_prepdat(timesmat = timesmat, statusmat = statusmat, id = id, 
                 starting_time = starting_time, starting_state = starting_state, 
                 transitionmat = transitionmat, 
                 original_states = (1:nrow(transitionmat)), longmat = NULL)
  out <- as.data.frame(out)
  names(out) <- c(idnam, "from_state", "to_state", "transition", 
                  "Tstart", "Tstop", "status")
  out$time <- out$Tstop - out$Tstart
  out <- out[, c(1:6, 8, 7)]
  ord <- order(out[, 1], out[, 5], out[, 2], out[, 3])
  out <- out[ord, ]
  row.names(out) <- 1:nrow(out)
  # Covariates
  if (!is.null(covs)) {
    n_covs <- length(covs)
    cov_cols <- match(covs, names(data))
    cov_names <- covs
    covs <- data[, cov_cols]
    if (!is.factor(out[, 1])) 
      out[, 1] <- factor(out[, 1])
    n_per_subject <- tapply(out[, 1], out[, 1], length)
    if (n_covs > 1) 
      covs <- covs[order_id, , drop = FALSE]
    if (n_covs == 1) {
      longcovs <- rep(covs, n_per_subject)
      longcovs <- longcovs[ord]
      longcovs <- as.data.frame(longcovs)
      names(longcovs) <- cov_names
    } else {
      longcovs <- lapply(1:n_covs, function(i) rep(covs[, i], n_per_subject))
      longcovs <- as.data.frame(longcovs)
      names(longcovs) <- cov_names
    }
    out <- cbind(out, longcovs)
  }
  # add attributes maybe
  # add specific class maybe
  # need to add functionality for covariates (e.g. like keep in mstate)
  return(out)
} 
