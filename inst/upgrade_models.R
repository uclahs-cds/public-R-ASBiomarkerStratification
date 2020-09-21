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

# https://topepo.github.io/caret/subsampling-for-class-imbalances.html
sampling = "up";
train.control.up <- train.control;
train.control.up$sampling <- "up";

bx.models.up <- AS.models(biodb, train.control.up);

saveRDS(bx.models.up, here('data/progressed_models_upsampled2.Rds'));

resamps <- resamples(list(GBM = bx.models$gbm.fit,
                          C5.0 = bx.models$c50.fit,
                          rpart = bx.models$rpart.fit))

resamps.up <- resamples(list(GBM.up = bx.models.up$gbm.fit,
                             C5.0.up = bx.models.up$c50.fit,
                             rpart.up = bx.models.up$rpart.fit))

resamps.joined <- resamples(list(GBM = bx.models$gbm.fit,
                                 C5.0 = bx.models$c50.fit,
                                 rpart = bx.models$rpart.fit,
                                 GBM.up = bx.models.up$gbm.fit,
                                 C5.0.up = bx.models.up$c50.fit,
                                 rpart.up = bx.models.up$rpart.fit))

(sum.resamps <- summary(resamps))

(sum.resamps.up <- summary(resamps.up))

(sum.resamps.joined <- summary(resamps.joined))
