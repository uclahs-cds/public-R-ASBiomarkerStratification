library(ProstateCancer.ASBiomarkerSynergy);

# Project : Prostate Cancer Active Surveillance

### biodb.AUPRC.R #######################################
# Perform AUPRC plotting

# Load Data:
data <- default.load.data();
attach(data);

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
  PSADensity = biodb$freePSA/biodb$ProstateVolume,
  p2PSA = biodb$p2PSA,
  PercentFreePSA = biodb$PercentFreePSA,
  PHI = biodb$PHI,
  PHIDensity = biodb$PHI/biodb$ProstateVolume,
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
  'PSADensity',
  'p2PSA',
  'PercentFreePSA',
  'PHI',
  'PHIDensity',
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

# real.names.set <c(
#   'PCA3',
#   'T2ERG',
#   'MiPS Cancer Risk',
#   'MiPS High Grade Cancer Risk',
#   'PSA Hybrid',
#   'free PSA',
#   'p2PSA',
#   'Percent Free PSA',
#   'PHI',
#   'Genetic Ancestry',
#   'Genetic Risk Score',
#   'Genetic Risk Category',
#   'Global Screening Array',
#   'GSA Positives',
#   'BRCA Mutation',
#   'Mutation 1',
#   'Mutation 2',
#   'TNFa Average',
#   'ADC normal Signal',
#   'ADC lesion Signal',
#   'RSI lesion Signal',
#   'RSI normal Signal',
#   'SOCPSA'
#   );

# Testing Purposes:
# variables.set <- 'PCA3';

# Saving Data for Operating Points:
common.size <- rep(NA, length(variables.set));
operating.data <- data.frame(
  test = common.size,
  tp = common.size,
  fp = common.size,
  tn = common.size,
  fn = common.size,
  recalls = common.size,
  fpr = common.size,
  accuracies = common.size,
  precisions = common.size,
  f1.scores = common.size,
  specificity = common.size
);

# Defining Steps for the thresholds:
step <- 1;

# Making the PRC plots:
prcplots.set <- list();
rocplots.set <- list();
auprc <- rep(NA, length(variables.set));
ci.auprc.upper <- rep(NA, length(variables.set));
ci.auprc.lower <- rep(NA, length(variables.set));

auroc <- auprc;
ci.auroc.upper <- ci.auprc.upper;
ci.auroc.lower <- ci.auprc.lower;

auprc.label.coord <- c(
  0.75,
  0.71,
  0.39,
  0.28,
  0.66,

  0.66,
  0.48,
  0.72,
  0.38,
  0.80,
  0.45,
  0.38,


  0.33,
  0.24,
  0.21,
  0.47,
  0.45,

  0.59,
  0.60,
  0.48,
  0.35,
  0.38,

  0.4,
  0.39,
  0.63
);

auroc.label.coord <- c(
  0.02,
  0.02,
  0.02,
  0.02,
  0.02,
  0.02,
  0.02,
  0.02,
  0.02,
  0.02,
  0.02,
  0.02,
  0.02,
  0.02,
  0.02,
  0.02,
  0.02,
  0.02,
  0.02,
  0.02,
  0.02,
  0.02,
  0.02,
  0.02,
  0.02
);

index <- 1;
#sink('AUPRC.log');
for (variable in variables.set){
  print(paste(
    'AUPRC and PR for ',
    variable,
    sep = ''
  )
  );

  select <- c('ID', 'BiopsyUpgraded',  variable);
  biomarker.test <- data.set[,c(select)];
  colnames(biomarker.test) <- c('ID', 'actual', 'predicted');
  biomarker.test <- biomarker.test[complete.cases(biomarker.test),];
  biomarker.test <- biomarker.test[order(biomarker.test$predicted),];
  biomarker.test <- biomarker.test[0 <= biomarker.test$predicted,];

  # Original scores
  original <- biomarker.test$predicted;

  # Rescaling the predicted signal:
  numerator <- biomarker.test$predicted - min(biomarker.test$predicted);
  denominator <- max(biomarker.test$predicted) - min(biomarker.test$predicted);
  biomarker.test$predicted <- 100*numerator/denominator;

  lower.threshold <- min(biomarker.test$predicted) - step;
  upper.threshold <- max(biomarker.test$predicted) + step;
  print(paste(
    'Test ',
    variable,
    ' Lower Threshold : ',
    lower.threshold,
    sep = ''
  )
  );
  print(paste(
    'Test ',
    variable,
    ' Upper Threshold : ',
    upper.threshold,
    sep = ''
  )
  );
  thresholds.set <- seq(round(lower.threshold), round(upper.threshold), step);

  # Defining data frame for merit data:
  common.size <- rep(NA, length(thresholds.set));
  merit.data <- data.frame(
    threshold = thresholds.set,
    scores = common.size,
    tp = common.size,
    fp = common.size,
    tn = common.size,
    fn = common.size,
    recalls = common.size,
    fpr = common.size,
    accuracies = common.size,
    precisions = common.size,
    f1.scores = common.size,
    specificity = common.size,
    positive.list = common.size,
    positive.patients = common.size
  );

  ## Computing the Confusion Matrix:
  tp <- rep(NA, length(thresholds.set));
  tn <-	rep(NA, length(thresholds.set));
  fp <- rep(NA, length(thresholds.set));
  fn <- rep(NA, length(thresholds.set));

  # NOTE:
  # Confusion Matrix is going to be built.
  # We are expecvariable,ting something like this:
  #         Predicted
  #            0     1
  # Actual  0  tn   fp
  #         1  fn   tp
  #
  # It is a 2x2 matrix
  # Actual currently has 0 and 1. Predicted MAY NOT
  # have both. Thus, the matrix can be different in
  # certain scenarios.

  #print(thresholds.set)

  merit.data$threshold <- thresholds.set;
  merit.data$scores <- thresholds.set*denominator/100 + min(original);

  for (j in 1:length(thresholds.set)) {
    threshold <- thresholds.set[j];
    data <- biomarker.test;
    data$predicted[threshold <= data$predicted] <- 1;
    data$predicted[1 != data$predicted] <- 0;

    if (0 < length(which(1 == data$predicted))){

      # Predicted data with values of 1 when compared to actual data, two outcomes are
      # possible match with the actual value (True Positive) or miss it (False Positive)
      positive.index <- which(1 == data$predicted)
      actual.positive <- which(1 == data$actual);

      indices <- actual.positive[which(actual.positive %in% positive.index)]
      merit.data$positive.list[j] <- paste(original[indices], collapse=',');

      merit.data$positive.patients[j] <- paste(as.character(biomarker.test$ID[indices]) , collapse=',');

      # If indices for positive in predicted data are present in actual data, they are
      # TRUE POSITIVE:

      tp[j] <- length(which(actual.positive %in% positive.index));
      # The remaining is FALSE POSITIVE:
      fp[j] <- length(positive.index) - tp[j];
    }
    else {
      # This is the case when NO POSITIVE is present in predicted data
      tp[j] <- 0;
      fp[j] <- 0;
    }

    if (0 < length(which(0 == data$predicted))){

      # Similar as above. Two outcomes are possible: Predicted value of 0 that hits the actual
      # data are TRUE NEGATIVE. Otherwise, FALSE NEGATIVE.
      negative.index <- which(0 == data$predicted)
      actual.negative <- which(0 == data$actual);
      tn[j] <- length(which(actual.negative %in% negative.index));
      fn[j] <- length(negative.index) - tn[j];
    }
    else {
      # This is the case when NO NEGATIVE is present in predicted data
      tn[j] <- 0;
      fn[j] <- 0;
    }

  }

  # Generating Merit Data and save it in a R object:
  merit.data$tp <- tp;
  merit.data$tn <- tn;
  merit.data$fp <- fp;
  merit.data$fn <- fn;
  merit.data$specificity <- tn/(tn + fp);
  merit.data$fpr <- fp/(fp + tn);
  merit.data$precisions <- tp/(tp + fp);
  merit.data$recalls <- tp/(tp + fn);
  merit.data$accuracies <- (tp + tn)/(tp + fp + tn + fn);
  merit.data$f1.scores <- 2*merit.data$recalls*merit.data$precisions/(merit.data$recalls + merit.data$precisions);
  merit.data[which(0 == merit.data$recalls),]$precisions <- 1;
  merit.data[which(0 == merit.data$recalls),]$f1.scores <- 0;

  # Finding the operating point:
  operating.point <- length(merit.data$specificity[merit.data$specificity<0.95]) + 1;
  operating.data$test[index] <- variable;
  operating.data$threshold[index] <- merit.data[operating.point, ]$threshold;
  operating.data$tp[index] <- merit.data[operating.point, ]$tp;
  operating.data$tn[index] <- merit.data[operating.point, ]$tn;
  operating.data$fp[index] <- merit.data[operating.point, ]$fp;
  operating.data$fn[index] <- merit.data[operating.point, ]$fn;
  operating.data$specificity[index] <- merit.data[operating.point, ]$specificity;
  operating.data$fpr[index] <- merit.data[operating.point, ]$fpr;
  operating.data$precisions[index] <- merit.data[operating.point, ]$precisions;
  operating.data$recalls[index] <- merit.data[operating.point, ]$recalls;
  operating.data$accuracies[index] <- merit.data[operating.point, ]$accuracies;
  operating.data$f1.scores[index] <- merit.data[operating.point, ]$f1.scores;

  # Saving the Merit data as backup:
  output.filename <- BoutrosLab.utilities::generate.filename(
    project.stem = 'AS',
    file.core = variable,
    extension = 'RData',
    file.date = Sys.Date()
  );
  save(
    merit.data,
    file = output.filename
  );

  # Plotting the PR curve for each test:
  prcplot.data <- data.frame(
    recall = merit.data$recalls,
    precision = merit.data$precisions
  );

  # Computing area under the curve using the trapezoid rule:
  #area.data <- unique(prcplot.data[complete.cases(prcplot.data),]);
  area.data <- prcplot.data;
  auprc[index] <- compute.auc(area.data, 'precision', 'recall');
  print(paste('AUPRC for ' , variable, ' ', auprc[index], sep = ''));

  # Confidence Interval for AUPRC:
  ci <- compute.auc.ci(auprc[index], biomarker.test$actual);
  ci.auprc.upper[index] <- ci[1];
  ci.auprc.lower[index] <- ci[2];

  roc1 <- pROC::roc(biomarker.test$actual, biomarker.test$predicted);
  print(paste('Hello ', pROC::ci.auc(roc1,  conf.level=0.95, boot.n=4000, boot.stratified=TRUE, reuse.auc=TRUE), sep=''));

  # Assigning feature for each panel:
  if (length(variables.set) - 4 > index){
    xaxislab <- NULL;
    xtck <- c(0, 0);
  }
  else {
    xaxislab <- seq(0, 1.1, 0.2);
    xtck <- c(1, 0);
  }

  if (index %in% c(1, 6, 11, 16, 21)){
    yaxislab <- seq(0, 1.1, 0.2);
    ytck <- c(1, 0);
  }
  else {
    yaxislab <- NULL;
    ytck <- c(0, 0);
  }

  prcplots.set[[index]] <- make.aucplot(
    aucplot.data = prcplot.data,
    #aucplot.main = variable,
    aucplot.xlabel = NULL,
    aucplot.ylabel = NULL,
    aucplot.xaxislab = xaxislab,
    aucplot.yaxislab = yaxislab,
    aucplot.xtck = xtck,
    aucplot.ytck = ytck,
    aucplot.ylimits = c(0, 1.2),
    aucplot.text = variable,
    aucplot.coord.text = c(auprc.label.coord[index], 0.95),
    aucplot.info = c(paste('AUPRC: ', round(auprc[index],3), sep = '')),
    aucplot.coord.info = c(0.74, 0.95),
    aucplot.output = 'NULL',
    aucplot.variable = precision~recall,
    aucplot.line = FALSE
  );

  # Plotting the PR curve for each test:
  rocplot.data <- data.frame(
    recall = merit.data$recalls,
    fpr = merit.data$fpr
  );

  # AUROC data
  ci.auroc.lower[index] <- pROC::ci(pROC::roc(biomarker.test$actual, biomarker.test$predicted))[1];
  auroc[index] <- pROC::ci(pROC::roc(biomarker.test$actual, biomarker.test$predicted))[2];
  ci.auroc.upper[index] <- pROC::ci(pROC::roc(biomarker.test$actual, biomarker.test$predicted))[3];

  rocplots.set[[index]] <- make.aucplot(
    aucplot.data = rocplot.data,
    #aucplot.main = variable,
    aucplot.xlabel = NULL,
    aucplot.ylabel = NULL,
    aucplot.xaxislab = xaxislab,
    aucplot.yaxislab = yaxislab,
    aucplot.xtck = xtck,
    aucplot.ytck = ytck,
    aucplot.ylimits = c(0, 1),
    aucplot.text = variable,
    aucplot.coord.text = c(auroc.label.coord[index], 0.97),
    aucplot.info = NULL,
    aucplot.coord.info = NULL,
    aucplot.output = 'NULL',
    aucplot.variable = recall~fpr,
    aucplot.line = TRUE
  );

  index = index + 1;

}
