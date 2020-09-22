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

saveRDS(prog.mod, here('data/progressed_models_acc.Rds'));

### Compare rpart, C5.0, and GBM
resamps <- resamples(list(GBM = prog.mod$gbm.fit,
                             C5.0 = prog.mod$c50.fit,
                             rpart = prog.mod$rpart.fit))

(sum.resamps <- summary(resamps))

prog.mod.F <- AS.models(biodb, train.control, metric = 'F');

saveRDS(prog.mod.F, here('data/progressed_models_F.Rds'));

### Compare rpart, C5.0, and GBM
resamps.F <- resamples(list(F.GBM = prog.mod.F$gbm.fit,
                          F.C5.0 = prog.mod.F$c50.fit,
                          F.rpart = prog.mod.F$rpart.fit))

(sum.resamps.F <- summary(resamps.F))

trellis.par.set(caretTheme())
dotplot(resamps, metric = c("F", "Sens", "Spec", "Accuracy"))

diff.resamps <- diff(resamps)
diff.resamps.F <- diff(resamps.F)
bwplot(resamps)
bwplot(resamps.F, layout = c(3, 1))

rpart.plot(prog.mod$rpart.fit$finalModel);
plot(prog.mod$c50.bestfit)

# Variable importance
lapply(list(prog.mod$rpart.fit, prog.mod$c50.fit, prog.mod$gbm.fit), varImp)
# do.call(cbind, lapply(list(rpart.fit, c50.fit, gbm.fit.progress), function(x) varImp(x)$importance))

summary(
    resamples(list(GBM = prog.mod$gbm.fit,
                   C5.0 = prog.mod$c50.fit,
                   rpart = prog.mod$rpart.fit,
                   F.GBM = prog.mod.F$gbm.fit,
                   F.C5.0 = prog.mod.F$c50.fit,
                   F.rpart = prog.mod.F$rpart.fit))
)
