
SS1 <- transf_eta(eta_H[[1]][, 2], Funs_FunForms[[1]][[2]])
SS2 <- transf_etaC(eta_H[[1]][, 1, drop= F], Funs_FunForms[[1]][[1]])

all.equal(SS1, SS2)

TT1 <- create_Wlong(eta_H, FunForms_per_outcome, U_H, Funs_FunForms)
TT2 <- create_WlongC(eta_H, U_H, model_info)

all.equal(TT1, TT2, check.attributes = FALSE)

library("rbenchmark")
benchmark(R = create_Wlong(eta_H, FunForms_per_outcome, U_H, Funs_FunForms),
          Cpp = create_WlongC(eta_H, U_H, model_info),
          replications = 10000)
