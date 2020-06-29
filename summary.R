
### pair.testing.loop.R ###################################
# The following R code take a set tests with their
# respective thresholds and perfoms a serial test

# Cleaning Environment:
rm(list=ls());

# Loading Data:
source('load.data.R');

#############################################
# Create data frame.       
#############################################
data.set <- data.frame(
  ID = as.character(biodb$Record.ID),
  BiopsyUpgraded = biodb$BiopsyUpgraded,
  PCA3 = biodb$PCA3,
  T2ERG = biodb$T2ERG,
  MiPSCancerRisk = biodb$MiPSCancerRisk,
  MiPSHGCancerRisk = biodb$MiPSHighGradeCancerRisk,
  PSAHyb = biodb$PSAHyb,
  freePSA = biodb$freePSA,
  p2PSA = biodb$p2PSA,
  PercentFreePSA = biodb$PercentFreePSA,
  PHI = biodb$PHI,
  GeneticAncestry = biodb$GeneticAncestry,
  GeneticRiskScore = biodb$GeneticRiskScore,
  GeneticRiskCategory = biodb$GeneticRiskCategory,
  GlobalScreeningArray = biodb$GlobalScreeningArray,
  GSAPositives = biodb$GSAPositives,
  BRCAMutation = biodb$BRCAMutation,
  Mutation1 = biodb$Mutation1,
  Mutation2 = biodb$Mutation.2,
  TNFaAverage = biodb$TNFaAverage,
  ADCnormalSignal = biodb$ADCnormalSignal,
  ADClesionSignal = biodb$ADClesionSignal,
  RSIlesionSignal = biodb$RSIlesionSignal,
  RSInormalSignal = biodb$RSInormalSignal,
  SOCPSA = biodb$SOCPSA
  );

# Defining Response and variables set:
response <- 'BiopsyUpgraded';
variables.set <- c(
  'PCA3',
  'T2ERG',
  'MiPSCancerRisk',
  'MiPSHGCancerRisk',
  'PSAHyb',
  'freePSA',
  'p2PSA',
  'PercentFreePSA',
  'PHI',
  'GeneticAncestry',
  'GeneticRiskScore',
  'GeneticRiskCategory',
  'GlobalScreeningArray',
  'GSAPositives',
  'BRCAMutation',
  'Mutation1',
  'Mutation2',
  'TNFaAverage',
  'ADCnormalSignal',
  'ADClesionSignal',
  'RSIlesionSignal',
  'RSInormalSignal',
  'SOCPSA'
  );

# The previous order comes from the spreadsheet.
# The following is by type:
variables.by.type <- c(
  'ADCnormalSignal',
  'ADClesionSignal',
  'RSIlesionSignal',
  'RSInormalSignal',
  'PCA3',
  'T2ERG',
  'MiPSCancerRisk',
  'MiPSHGCancerRisk',
  'PSAHyb',
  'freePSA',
  'p2PSA',
  'PercentFreePSA',
  'PHI',
  'SOCPSA',
  'TNFaAverage',
  'GeneticAncestry',
  'GeneticRiskScore',
  'GeneticRiskCategory',
  'GlobalScreeningArray',
  'GSAPositives',
  'BRCAMutation',
  'Mutation1',
  'Mutation2' 
  );

real.names.set <- c(
  'ADC normal Signal',
  'ADC lesion Signal',
  'RSI lesion Signal',
  'RSI normal Signal',
  'PCA3',
  'T2ERG',
  'MiPS Cancer Risk',
  'MiPS High Grade Cancer Risk',
  'PSA Hybrid',
  'free PSA',
  'p2PSA',
  'Percent Free PSA',
  'PHI',
  'SOCPSA',
  'TNFa Average',
  'Genetic Ancestry',
  'Genetic Risk Score',
  'Genetic Risk Category',
  'Global Screening Array',
  'GSA Positives',
  'BRCA Mutation',
  'Mutation 1',
  'Mutation 2' 
  );

#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 

common.size <- rep(NA, length(variables.set));
summary.data <- data.frame(
  test.name = common.size,
  min.value = common.size,
  first.quarter = common.size,
  median = common.size,
  mean = common.size,
  third.quarter	= common.size,
  max.values = common.size,
  na.values = common.size
  );

for (i in 1:length(variables.set)){
  deliver <- paste('Processing Test ', i , ' : ', real.names.set[i], sep = ' ');
  print(deliver);
  summary.data$test.name[i] = real.names.set[i];
  summary.data$min.value[i] = round(summary(data.set[, variables.set[i]])[1],2); 
  summary.data$first.quarter[i] = round(summary(data.set[, variables.set[i]])[2],2);
  summary.data$median[i] = round(summary(data.set[, variables.set[i]])[3],2);
  summary.data$mean[i] = round(summary(data.set[, variables.set[i]])[4],2);
  summary.data$third.quarter[i] = round(summary(data.set[, variables.set[i]])[5],2);
  summary.data$max.values[i] = round(summary(data.set[, variables.set[i]])[6],2);
  summary.data$na.values[i] = round(summary(data.set[, variables.set[i]])[7],2); 
  };

### WRITE SESSION PROFILE TO FILE #####################
save.session.profile(
  BoutrosLab.utilities::generate.filename(
    Sys.Date(),
    'pair.testing.loop',
    'txt'
    )
  );

