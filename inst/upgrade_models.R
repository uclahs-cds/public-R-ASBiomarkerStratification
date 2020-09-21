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

bx.models <- AS.models(biodb, train.control, target = 'BiopsyUpgraded');

saveRDS(bx.models, here('data/bx_upgrade_models.Rds'));

upgrade.resamps <- resamples(list(GBM = bx.models$gbm.fit,
                          C5.0 = bx.models$c50.fit,
                          rpart = bx.models$rpart.fit))

(sum.upgrade.resamps <- summary(upgrade.resamps))
