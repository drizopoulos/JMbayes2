library(Rcpp)
library(RcppArmadillo)
library(microbenchmark)
library(ggplot2)

sourceCpp(file = 'D:\\gpapageorgiou\\Github_Repos\\JMbayes2\\Dev_Local\\Cpp_funs\\propose_rnorm.cpp')
sourceCpp(file = 'D:\\gpapageorgiou\\Github_Repos\\JMbayes2\\Dev_Local\\Cpp_funs\\propose_rnorm_3.cpp')
sourceCpp(file = 'D:\\gpapageorgiou\\Github_Repos\\JMbayes2\\Dev_Local\\Cpp_funs\\propose_rnorm_2.cpp')
sourceCpp(file = 'D:\\gpapageorgiou\\Github_Repos\\JMbayes2\\Dev_Local\\Cpp_funs\\chol_cube_S.cpp')
load(file = 'D:\\gpapageorgiou\\Github_Repos\\JMbayes2\\Dev_Local\\Cpp_funs\\propose_rnorm.RData')

#b.rows <- max(do.call(c, lapply(b, nrow)))
#b.cols <- do.call(c, lapply(b, ncol))

#sigmas <- rep(0.6 / sum(b.cols), b.rows)

set.seed(357)
chk1 <- propose_mvnorm_cube(1, vcov_prop_RE, sigmas)
set.seed(357)
chk2 <- propose_mvnorm_cube_2(1, vcov_prop_RE, sigmas)
all.equal(chk1, chk2)

set.seed(257)
#sigmas <- runif(dim(vcov_prop_RE)[3])
chk1 <- propose_mvnorm_cube(1, vcov_prop_RE, sigmas)
set.seed(257)
#sigmas <- runif(dim(vcov_prop_RE)[3])
chk2 <- propose_mvnorm_cube(1, vcov_prop_RE, sigmas)
all.equal(chk1, chk2)


S <- vcov_prop_RE
S <- chol_cube(vcov_prop_RE)
dim(S)

S2 <- vcov_prop_RE
S2 <- chol_cube(vcov_prop_RE)

S[,,1]
S2[,,2]

set.seed(1)
chk1 <- propose_mvnorm_mat(1, S, sigmas)
set.seed(1)
chk2 <- propose_mvnorm_mat(1, S2, sigmas)

all.equal(chk1, chk2)
