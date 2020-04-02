# Source Functions
source(file.path(getwd(), "R/jm.R"))
source(file.path(getwd(), "R/help_functions.R"))
source(file.path(getwd(), "Development/jm/R_to_Cpp.R"))
source(file.path(getwd(), "Development/jm/PBC_data.R"))

# load data
load(file = file.path(getwd(), "/Dev_Local/sample_case_env_02042020.RData"))

# note that for missing subjects for a long outcome the vector of random-effects
# equals zero
# lapply(idL, FUN = function(x) length(unique(x)))
# lapply(b, nrow)
# which(!unique(idL[[1]]) %in% unique(idL[[2]]))
# which(b[[2]][, 2] == 0)

#---------------------------------------------------------------
#                        FUNCTIONS
#---------------------------------------------------------------
# sample random values from multivariate normal distribution
mvrnorm_gp <- function(n, S) {
  A <- chol(S)
  z <- matrix(rnorm(n * ncol(S)), nrow = n, ncol = ncol(S))
  Y <- z %*% t(A)
  Y
}

# sample random values from t distribution

# MCMC
M <- 3000
bs <- list(NULL)
b.rows <- max(do.call(c, lapply(b, nrow)))
b.cols <- do.call(c, lapply(b, ncol))

for (m in seq_len(M)) {
  accept_b <- matrix(0.0, nrow = b.rows, ncol = M)
  
}