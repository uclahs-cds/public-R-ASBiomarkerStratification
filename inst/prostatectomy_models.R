library(caret)
library(rpart)
library(rpart.plot)
library(C50)
library(ProstateCancer.ASBiomarkerSynergy)
library(pROC)
library(here)
# Boosted Trees
library(gbm)
library(partykit)

set.seed(1333)

biodb <- default.load.data(onlyBiodb = TRUE);

train.control <- trainControl(
    method = "repeatedcv",
    number = 10,
    repeats = 5
)

prostatectomy.models <- AS.models(biodb, train.control, target = 'Prostatectomy');

saveRDS(prostatectomy.models, here('data/Prostatectomy_models.Rds'));

prostatectomy.resamps <- resamples(list(GBM = prostatectomy.models$gbm.fit,
                                  C5.0 = prostatectomy.models$c50.fit,
                                  rpart = prostatectomy.models$rpart.fit))

(sum.prostatectomy.resamps <- summary(prostatectomy.resamps))

prostatectomy.models.F <- AS.models(biodb, train.control, target = 'Prostatectomy', metric = 'F');

saveRDS(prostatectomy.models.F, here('data/Prostatectomy_models_F.Rds'));

prostatectomy.resamps.F <- resamples(list(GBM.F = prostatectomy.models.F$gbm.fit,
                                        C5.0.F = prostatectomy.models.F$c50.fit,
                                        rpart.F = prostatectomy.models.F$rpart.fit))

(sum.prostatectomy.resamps.F <- summary(prostatectomy.resamps.F))


