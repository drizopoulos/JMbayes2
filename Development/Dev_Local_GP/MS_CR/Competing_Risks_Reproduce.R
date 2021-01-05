library(JM)
library(JMbayes)
library(JMbayes2)

##############################
## Do things the JM, CR way ##
## Limited to 1 longitudinal##
## process                  ##
##############################

## slides: http://www.drizopoulos.com/courses/EMC/ESP72.pdf
## slide 163
pbc2.idCR <- crLong(pbc2.id, statusVar = "status",
                    censLevel = "alive", nameStrata = "CR")
pbc2.idCR[pbc2.idCR$id %in% c(1,2,5),
          c("id", "years", "status", "CR", "status2")]
## slide 164
coxFit.CR <- coxph(Surv(years, status2) ~ drug * strata(CR),
                   data = pbc2.idCR, x = TRUE)
## equivalently:
pbc2.idCR$drug.1 <- (pbc2.idCR$drug=="D-penicil" & pbc2.idCR$CR=="transplanted")+0L
pbc2.idCR$drug.2 <- (pbc2.idCR$drug=="D-penicil" & pbc2.idCR$CR=="dead"        )+0L
coxFit.CR2<- coxph(Surv(years, status2) ~ drug.1 + drug.2 + strata(CR),
                   data = pbc2.idCR, x = TRUE)


## slide 164 & 165
lmeFit.CR <- lme(log(serBilir) ~ drug * year, data = pbc2,
                 random = ~ year | id)

jointFit.CR <- jointModel(lmeFit.CR, coxFit.CR, timeVar = "year",
                          method = "spline-PH-aGH", CompRisk = TRUE,
                          interFact = list(value = ~ CR, data = pbc2.idCR))

## same as jointFit.CR -- different parameterization (use coxFit.CR2 and change interFact).
jointFit.CR2<- jointModel(lmeFit.CR, coxFit.CR2, timeVar = "year",
                          method = "spline-PH-aGH", CompRisk = TRUE,
                          interFact = list(value = ~ strata(CR) - 1, data = pbc2.idCR))



###################################
## Do things the JMbayes, MS way ##
## Allows competing risks and >1 ##
## longitudinal process          ##
###################################
library(mstate)


## follow section 5 of
## https://cran.r-project.org/web/packages/mstate/vignettes/Tutorial.pdf
pbc2.id[pbc2.id$id %in% c(1,2,5),
        c("id", "years", "status", "status2")]


## prepare tmat
tmat <- trans.comprisk(2, names = c("alive", "transplanted", "dead"))
tmat
## use msprep()
pbc2.id$statTX <- as.numeric(pbc2.id$status=="transplanted")
pbc2.id$statDD <- as.numeric(pbc2.id$status=="dead")
?msprep
pbc2.idMS <- msprep(time =   c(NA, "years", "years"),
                    status = c(NA, "statTX", "statDD"),
                    data = pbc2.id,
                    keep = "drug",
                    trans = tmat)


## Crucial step!  I omitted this in my previous comment.
## Expand covariates per transition. See before/after expansion.
pbc2.idMS[pbc2.idMS$id %in% c(1,2,5),]
pbc2.idMS <- expand.covs(pbc2.idMS, "drug", append = TRUE, longnames = FALSE)
pbc2.idMS[pbc2.idMS$id %in% c(1,2,5),]


## explicitly put in the expanded covariates and add strata() and add cluster()
coxFit.MS <- coxph(Surv(Tstart, Tstop, status) ~ drug.1 + drug.2 + strata(trans) + cluster(id),
                   data = pbc2.idMS, model=TRUE, x=TRUE)

## reproduce jointFit.CR with mvglmer/mvJointModelBayes (?)
mvglmerFit.MS <- mvglmer(list(log(serBilir) ~ drug * year + (year|id)),
                         data = pbc2,
                         families=list(gaussian)
)


## Crucial step:
## Specify interactions in order to allow different effect of the outcome for each risk
interacts <- list("log(serBilir)" = ~ strata(trans) - 1)

jointFit.MS <-
  mvJointModelBayes(mvglmerFit.MS,
                    coxFit.MS,
                    timeVar="year",
                    Interactions=interacts,
                    multiState=TRUE,
                    data_MultiState=pbc2.idMS,
                    idVar_MultiState = "id",
                    control = list(equal.strata.knots=TRUE,
                                   equal.strata.bound.knots=TRUE)
  )

###################################
## Do things the JMbayes2, MS way##
## Allows competing risks and >1 ##
## longitudinal process          ##
###################################
fforms <- list("log(serBilir)" = ~ value(log(serBilir)) + value(log(serBilir)):(strata(trans) - 1))

fitjm2 <- jm(Surv_object = coxFit.MS, Mixed_objects = lmeFit.CR, time_var = 'year', 
             data_Surv = pbc2.idMS, id_var = 'id', functional_forms = fforms)


## similar:
summary(jointFit.CR2)
summary(jointFit.MS)


## similar:
coxFit.CR
coxFit.MS

## similar:
summary(lmeFit.CR)
summary(mvglmerFit.MS)
