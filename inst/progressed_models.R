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

prog.mod <- AS.models(biodb, train.control);

saveRDS(prog.mod, here('data/progressed_models.Rds'));

# https://topepo.github.io/caret/subsampling-for-class-imbalances.html
sampling = "up";
train.control.up <- train.control;
train.control.up$sampling <- "up";

prog.mod.up <- AS.models(biodb, train.control.up);

saveRDS(prog.mod.up, here('data/progressed_models_upsampled2.Rds'));

### Compare rpart, C5.0, and GBM
resamps.joined <- resamples(list(GBM = prog.mod$gbm.fit,
                          C5.0 = prog.mod$c50.fit,
                          rpart = prog.mod$rpart.fit,
                          GBM.up = prog.mod.up$gbm.fit,
                          C5.0.up = prog.mod.up$c50.fit,
                          rpart.up = prog.mod.up$rpart.fit))

resamps <- resamples(list(GBM = prog.mod$gbm.fit,
                             C5.0 = prog.mod$c50.fit,
                             rpart = prog.mod$rpart.fit))

resamps.up <- resamples(list(GBM.up = prog.mod.up$gbm.fit,
                          C5.0.up = prog.mod.up$c50.fit,
                          rpart.up = prog.mod.up$rpart.fit))

(sum.resamps <- summary(resamps))

(sum.resamps.up <- summary(resamps.up))

(sum.resamps.joined <- summary(resamps.joined))

trellis.par.set(caretTheme())
dotplot(resamps, metric = c("F", "Sens", "Spec", "Accuracy"))
dotplot(resamps.up, metric = c("F", "Sens", "Spec", "Accuracy"))
dotplot(resamps.joined, metric = c("F", "Sens", "Spec", "Accuracy"))

diff.resamps <- diff(resamps)
bwplot(diff.resamps, layout = c(3, 1))

rpart.plot(prog.mod$rpart.fit$finalModel);
plot(prog.mod$c50.bestfit)

# Variable importance
lapply(list(prog.mod$rpart.fit, prog.mod$c50.fit, prog.mod$gbm.fit), varImp)
# do.call(cbind, lapply(list(rpart.fit, c50.fit, gbm.fit.progress), function(x) varImp(x)$importance))


