library(ProstateCancer.ASBiomarkerSynergy);
library(gtsummary)

biodb <- default.load.data(onlyBiodb = TRUE);

biomarkers <- load.biomarker.categories();
demo.vars <- biomarkers$variables[biomarkers$clinically.useful]

reduced.vars <- c(
  'Age',
  'BMI',
  'Race',
  'Ethnicity',
  'ProstateVolume',
  'PCA3',
  'T2ERG',
  'MiPSCancerRisk',
  'MiPSHighGradeCancerRisk',
  'PercentFreePSA',
  'PHI',
  'GeneticRiskScore',
  'RSIlesionSignal',
  'RSIlesionPIRADS',
  'PSADensity',
  'PHIDensity'
);

demo.vars <- c("Age", "Race", "Ethnicity", "Weight", "Height",
  "BMI", "SOCPSA", "MRIResult", "MRILesions", "HighestPIRADS",
  "BiopsyResult", "ProstateVolume", "PreviousGleason", "PreviousISUP",
  "StudyHighestGleason", "StudyHighestISUP", "Observation", "BiopsyUpgraded",
  "PCA3", "T2ERG", "MiPSCancerRisk", "MiPSHighGradeCancerRisk",
  "PSAHyb", "freePSA", "p2PSA", "PercentFreePSA", "PHI", "GeneticAncestry",
  "GeneticRiskScore", "GeneticRiskCategory", "GlobalScreeningArray",
  "GSAPositives", "BRCAMutation", "TNFaAverage", "TNFaSTD", "RSInormalSignal",
  "RSIlesionSignal", "ADCnormalSignal", "ADClesionSignal", "RSIlesionPIRADS",
  "RSIlesionCancer", "RSIlesionGleason", "RSIlesionISUP", "RSIlesionUpgraded",
  "RSIlesionObservation", "ProgressedToTreatment", "Prostatectomy",
  "UpgradedAndProgressed", "NoUpgradeAndProgressed", "DaysBxToProgression",
  "DaysBxToLastClinicalAppt", "DaysBxToLastReview", "DaysDxToUpgrade",
  "DaysDxToProgression", "DaysDxToLastClinicalAppt", "DaysDxToLastReview",
  "Mutation_BRCA1", "Mutation_BRCA2", "Mutation_ATM", "Mutation_MLH1",
  "Mutation_PMS2", "PSADensity", "PHIDensity", "Hispanic", "PHI.computed"
)

levels(biodb$BiopsyUpgraded) <- c('No Upgrade', 'Biopsy Upgraded')

tbl_summary(biodb[, c(reduced.vars, 'BiopsyUpgraded')], by = BiopsyUpgraded) %>%
  add_p()  %>%
  as_flex_table() %>%
  flextable::save_as_docx(path = here::here('euro_urology/tables/demographics.docx'))
