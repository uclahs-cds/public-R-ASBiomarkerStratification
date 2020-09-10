library(ProstateCancer.ASBiomarkerSynergy);

data <- default.load.data();
attach(data);

psa.cols <- c('p2PSA', 'freePSA', 'PSAHyb', 'PercentFreePSA', 'PHI');

# Missing decimal place
biodb[3, 'p2PSA'] <- biodb[3, 'p2PSA'] / 100;

# Confirm PHI formula
PHI.computed <- with(biodb, {
    (p2PSA / freePSA) * sqrt(PSAHyb)
})

cbind(biodb[, psa.cols], PHI.computed);

mean((biodb$PHI - PHI.computed)^2, na.rm = TRUE);

# Use `requireNamespace("table1", quietly = TRUE)` to test if package is installed
table1::table1(~ Age + Weight + Height + BMI | Race, data = biodb);

table1::table1(~ p2PSA + freePSA + PSAHyb + PercentFreePSA | BiopsyUpgraded, data = biodb);

table1::table1(~ GeneticAncestry | Race, data = biodb, overall = NULL, topclass="Rtable1-grid");


roc.phi <- pROC::roc(BiopsyUpgraded ~ PHI, data = biodb)
plot(roc.phi, legacy.axes = TRUE, smooth = TRUE)
