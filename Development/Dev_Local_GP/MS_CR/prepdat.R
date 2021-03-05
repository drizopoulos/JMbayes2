ms_prepdat <- function (timesmat, statusmat, id, starting_time, starting_state, transitionmat, 
                     original_states, longmat) {
  if (is.null(nrow(timesmat))) 
    return(longmat)
  if (nrow(timesmat) == 0) 
    return(longmat)
  from_states <- apply(!is.na(transitionmat), 2, sum)
  to_states <- apply(!is.na(transitionmat), 1, sum)
  absorbing_states <- which(to_states == 0)
  starts <- which(from_states == 0)
  new_states <- starting_state
  new_times <- starting_time
  rmv <- NULL
  for (i in 1:starts) {
    subjects <- which(starting_state == starts)
    n_start <- length(subjects)
    to_states_2 <- which(!is.na(transitionmat[starts, ]))
    trans_states <- transitionmat[starts, to_states_2]
    n_trans_states <- length(to_states_2)
    if (all(n_start > 0, n_trans_states > 0)) {
      Tstart <- starting_time[subjects]
      Tstop <- timesmat[subjects, to_states_2, drop = FALSE]
      Tstop[Tstop <= Tstart] <- Inf
      state_status <- statusmat[subjects, to_states_2, drop = FALSE]
      mintime <- apply(Tstop, 1, min)
      hlp <- Tstop * 1 / state_status
      hlp[Tstop == 0 & state_status == 0] <- Inf
      next_time <- apply(hlp, 1, min)
      censored <- which(is.infinite(apply(hlp, 1, min)))
      wh <- which(mintime < next_time)
      whminc <- setdiff(wh, censored)
      if (length(whminc) > 0) {
        whsubjs <- id[subjects[whminc]]
        whsubjs <- paste(whsubjs, collapse = " ")
        warning("Subjects ", whsubjs, " Have smaller transition time with status = 0, larger transition time with status = 1, 
                from starting state ", original_states[starting])
      }
      next_time[censored] <- mintime[censored]
      if (ncol(hlp) > 1) {
        hlpsrt <- t(apply(hlp, 1, sort))
        warn1 <- which(hlpsrt[, 1] - hlpsrt[, 2] == 0)
        if (length(warn1) > 0) {
          isw <- id[subjects[warn1]]
          isw <- paste(isw, collapse = " ")
          hsw <- hlpsrt[warn1, 1]
          hsw <- paste(hsw, collapse = " ")
          warning("simultaneous transitions possible for subjects ", isw, " at times ", hsw, 
                  " -> Smallest receiving state will be used")
        }
      }
      if (length(censored) > 0) {
        next_state <- apply(hlp[-censored, , drop = FALSE], 
                            1, which.min)
        absorbed <- (1:n_start)[-censored][which(to_states_2[next_state] %in% absorbing_states)]
      } else {
        next_state <- apply(hlp, 1, which.min)
        absorbed <- (1:n_start)[which(to_states_2[next_state] %in% absorbing_states)]
      }
      states_matrix <- matrix(0, n_start, n_trans_states)
      if (length(censored) > 0) {
        states_matrix_min <- states_matrix[-censored, , drop = FALSE]
      } else {
        states_matrix_min <- states_matrix
      }
      if (nrow(states_matrix_min) > 0) 
        states_matrix_min <- t(sapply(1:nrow(states_matrix_min), function(i) {
          x <- states_matrix_min[i, ]
          x[next_state[i]] <- 1
          return(x)
        }))
      if (length(censored) > 0) {
        states_matrix[-censored, ] <- states_matrix_min
      } else {
        states_matrix <- states_matrix_min
      }
      mm <- matrix(c(rep(id[subjects], rep(n_trans_states, n_start)), 
                     rep(original_states[starts], n_trans_states * n_start), 
                     rep(original_states[to_states_2], n_start), 
                     rep(trans_states, n_start), rep(Tstart, rep(n_trans_states, n_start)), 
                     rep(next_time, rep(n_trans_states, n_start)), as.vector(t(states_matrix))), 
                   n_trans_states * n_start, 7)
      longmat <- rbind(longmat, mm)
      rmv <- c(rmv, subjects[c(censored, absorbed)])
      if (length(censored) > 0) {
        new_states[subjects[-censored]] <- to_states_2[next_state]
      } else {
        new_states[subjects] <- to_states_2[next_state]
      }
      if (length(censored) > 0)  {
        new_times[subjects[-censored]] <- next_time[-censored]
      } else {
        new_times[subjects] <- next_time
      }
    }
  }
    if (length(rmv) > 0) {
      timesmat <- timesmat[-rmv, ]
      statusmat <- statusmat[-rmv, ]
      new_times <- new_times[-rmv]
      new_states <- new_states[-rmv]
      id <- id[-rmv]
    }
    n_states <- nrow(transitionmat)
    idx <- rep(1, n_states)
    idx[starts] <- 0
    idx <- cumsum(idx)
    new_states <- idx[new_states]
    Recall(timesmat = timesmat[, -starts], statusmat = statusmat[, -starts], 
           id = id, starting_time = new_times, starting_state = new_states, 
           transitionmat = transitionmat[-starts, -starts], original_states = original_states[-starts], 
           longmat = longmat)
}
