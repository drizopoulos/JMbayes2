fm1 <- lme(log(serBilir) ~ year * sex + I(year^2) + age + prothrombin,
           data = pbc2, random = ~ year | id)
fm2 <- lme(serChol ~ ns(year, 3, B = c(0, 10)) + sex + age, data = pbc2,
           random = ~ year | id, na.action = na.exclude)
fm2. <- lme(I(serChol / 150) ~ ns(year, 2, B = c(0, 10)) + sex + age, data = pbc2,
            random = ~ ns(year, 2, B = c(0, 10)) | id, na.action = na.exclude)
fm3 <- mixed_model(hepatomegaly ~ sex + age, data = pbc2,
                   random = ~ 1 | id, family = binomial())
fm4 <- mixed_model(ascites ~ year + age, data = pbc2,
                   random = ~ year | id, family = binomial())
fm5 <- lme(prothrombin ~ ns(year, 2, B = c(0, 10)) * drug, data = pbc2,
           random = ~ ns(year, 2, B = c(0, 10)) | id,
           control = lmeControl(opt = "optim"))
###
# D1

Mixed_objects <- list(fm1, fm3, fm4)

D_lis <- lapply(Mixed_objects, extract_D)
D1 <- bdiag(D_lis)

###
# D2

Mixed_objects <- list(fm1, fm2, fm3, fm4)

D_lis <- lapply(Mixed_objects, extract_D)
D2 <- bdiag(D_lis)

###
# D3

Mixed_objects <- list(fm1, fm2, fm3, fm4, fm5, fm1)

D_lis <- lapply(Mixed_objects, extract_D)
D3 <- bdiag(D_lis)

#########################

fm_s29_pbc <- gls(log(serBilir) ~ year + year:drug + year * sex + age, data = pbc2,
                  correlation = corCAR1(form = ~ year | id))

D4 <- getVarCov(fm_s29_pbc, individual = 271)
class(D4) <- "matrix"

D5 <- getVarCov(fm_s29_pbc, individual = 32)
class(D5) <- "matrix"

#########################

