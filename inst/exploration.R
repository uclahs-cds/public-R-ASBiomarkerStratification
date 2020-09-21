library(ProstateCancer.ASBiomarkerSynergy);

data <- default.load.data();

psa.cols <- c('p2PSA', 'freePSA', 'PSAHyb', 'PercentFreePSA', 'PHI');

# Missing decimal place
data$biodb[3, 'p2PSA'] <- data$biodb[3, 'p2PSA'] / 100;

# Confirm PHI formula
PHI.computed <- with(data$biodb, {
    (p2PSA / freePSA) * sqrt(PSAHyb)
})

cbind(data$biodb[, psa.cols], PHI.computed);

mean((data$biodb$PHI - PHI.computed)^2, na.rm = TRUE);

# Use `requireNamespace("table1", quietly = TRUE)` to test if package is installed
table1::table1(~ Age + Weight + Height + BMI | Race, data = data$biodb);

table1::table1(~ p2PSA + freePSA + PSAHyb + PercentFreePSA | BiopsyUpgraded, data = data$biodb);

table1::table1(~ GeneticAncestry | Race, data = data$biodb, overall = NULL, topclass="Rtable1-grid");

# roc.phi <- pROC::roc(BiopsyUpgraded ~ PHI, data = data$biodb)
# plot(roc.phi, legacy.axes = TRUE)
