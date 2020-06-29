

# Project : Prostate Cancer Active Surveillance

### AUPRC.R #################################
# Perform AUPRC and PR plotting and analysis
# using BiopsyUpgraded and Biomarkers tests

library(pROC);

run <- 2;
if (1 == run){
  rm(list=ls());
  source('install.packages.R');
  source('load.data.R');
  }

### make.confusion.matrix ################################
# Description: Given two datasets, predicted and actual,
# compute the confusion matrix. 
#
# Variables:
# predicted.data  data frame with values from test
# actual.data     data frame with binary data (0 and 1)
#
# Examples:
#
# make.confusion.matrix(
#   predicted.data = pca3.data,
#   actual.data = biopsy.upgraded.data
#   );
###########################################################
make.confusion.matrix <- function(
  predicted.data,
  actual.data
  ){
  # Get the array index for each positive (1) value in the actual data:        
  positives <- which(1 == actual.data, arr.ind = TRUE);
  # Number of Positives (1) in the actual data:
  num.positives <- length(positives);
  
  # Sort the predicted data in decreasing order and get the indices:
  ordered.predicted <- order(predicted.data, decreasing = TRUE);

  # NOTE: The reason we sort it that way is because the highest value
  # of the predicted dataset must be a POSITIVE CASE. 
  
  # Now, use these indices to find what values the actual data
  # is reporting (0 or 1). The cumulative sum is the number of
  # TRUE POSITIVES:
  tp <- cumsum(actual.data[ordered.predicted]);

  # Using the indices NOT present in the ordered list will give
  # the number of FALSE POSITIVES:
  fp <- cumsum(!actual.data[ordered.predicted]);
  
  # The number of positive minus the number of true positives
  # are the number of FALSE NEGATIVES:
  fn <- num.positives - tp;

  # Number of Negatives (0):
  num.negatives <- length(actual.data) - length(positives);
  # Some negatives will be detected (TRUE NEGATIVE)
  tn <- num.negatives -fp;
 
  # Building the Confusion Matrix:
  confusion.matrix.elements <- data.frame(
    tp = tp,
    fp = fp,
    tn = tn,
    fn = fn
    );
    biomarker.test
    
  return(confusion.matrix.elements);
  }

### compute.auc ######################################################
# Description: Compute area under the curve
#
# Variables:
# 
# input.data  data frame with 2 columns
# yaxis.label yaxis label,
# xaxis.label xaxis label,
#
# Example:
#
# compute.auc(input.data, 'recall', 'precision');
#################################################################
compute.auc <- function(
  input.data,
  yaxis.label,
  xaxis.label
  ){
  area <- c(0);
  for (k in 1:dim(input.data)[1]){
    if (1 < k){
      y.portion <- (input.data[, yaxis.label][k] + area.data[, yaxis.label][k-1]);
      x.portion <- (input.data[, xaxis.label][k-1] - area.data[,xaxis.label][k]);
      area <- area + y.portion*x.portion/2;
    }
    }
  return(area);
  }

### compute.auc.ci #############################################
# Description: Compute confidence interval of an area
#
# Variables:
#
# area.value  area value,
# actual.data actual data
#
# Example
#
# compute.auc.ci(auprc, actual);
#################################################################
compute.auc.ci <- function(
  area.value,
  actual.data
  ){
  auc <- area.value;
  q0 <- auc*(1 - auc);
  q1 <- auc/(2 - auc) - auc^2;
  q2 <- 2*auc^2/(1+auc) - auc^2;
  n1 <- length(which(0 == actual.data));
  n2 <- length(which(1 == actual.data));
  numerator <- q0 + (n1 - 1)*q1 + (n2 - 1)*q2;
  denominator <- n1*n2;
  limit <- sqrt(numerator/denominator);
  z.critical95 <- 1.96;
  limit <- limit*z.critical95;
  ci.range <- c(auprc[index] + limit, auprc[index] - limit);
  return(ci.range);
  }

## forest.plot #################################################
# Description: Make a forest plot
#
# Variables:
# forest.data         data frame
# forest.formula      formula to use for forest plot
# forest.xaxis.label  axis label
# forest.yaxis.label  yaxis label
# forest.yaxis.tck    yaxis tick
# forest.center       where we center the data
# forest.output       output filename
#
# Example:
# forest.plot(
#    forest.data = simple.data,
#    forest.formula = reorder(feature, AUPRC) ~ lower + upper,
#    forest.xaxis.label = 'AUPRC',
#    forest.yaxis.label = 'Test',
#    forest.center = simple.data$AUPRC,
#    forest.output = auprc_plot.file
#    );
################################################################
forest.plot <- function(
  forest.data,
  forest.formula,
  forest.xlabel,
  forest.ylabel,
  forest.yaxis.cex,
  forest.yaxis.tck,
  forest.center
  #forest.output
  ){
  BoutrosLab.plotting.general::create.segplot(
    formula = forest.formula,
    data = forest.data,
    xlab.label = forest.xlabel,
    ylab.label = forest.ylabel,
    xaxis.cex = 2,
    yaxis.cex = forest.yaxis.cex,
    xlab.cex = 3,
    ylab.cex = 2.,
    yaxis.tck = forest.yaxis.tck,
    xlimits = c(0, 1.0),
    xat = seq(0, 1.0, 0.25),
    draw.bands = FALSE,
    symbol.cex = 2,
    centers = forest.center,
    #filename = forest.output,
    description = 'Forest Plot for AUC data',
    resolution = 300
    );
  }

## make.aucplot #######################################################
# Description:       
# Function to plot PR Curves       
#
# aucplot.data  input data frame
# aucplot.xlabel title for the x-axis
# aucplot.ylabel title for the y-axis 
# aucplot.xaxislab labels for x-axis  
# aucplot.yaxislab labels for y-axis
# aucplot.xtck defines tick for x-axis c(a,b), a,b: 0 or 1
# aucplot.ytck defines tick for y-axis c(a,b), a,b: 0 or 1
# aucplot.ylimits defines yaxis limits 
# aucplot.text Text with some info
# aucplot.coord.text  coordinates of the text
# aucplot.info AUPR info
# aucplot.coord.info AUPR info coordinates
# aucplot.output output filename
# aucplot.formula  formula to use
# aucplot.line draw line (True or False)
#
# Example:     
#
# make.aucplot(
#   aucplot.data = prc.data,
#   aucplot.main = 'Main Title'
#   aucplot.xlabel = 'Recall',
#   aucplot.ylabel = 'Precision',
#   aucplot.xaxislab = seq(0, 1.1, 0.2),
#   aucplot.yaxislab = seq(0, 1.1, 0.2),
#   aucplot.xtck = c(1,0),
#   aucplot.ytck = c(1,0),
#   aucplot.ylimits = c(0,1),
#   aucplot.text = 'PCA3',
#   aucplot.coord.text = c(0.3, 0.15),
#   aucplot.info = c('AUC 123'),
#   aucplot.coords.info = c(0,0), 
#   aucplot.output = 'prc_plot.tiff',
#   aucplot.variable = 'precision ~ recall',
#   aucplot.line = FALSE
#   );
########################################################################
make.aucplot <- function(
  aucplot.data,
  #aucplot.main,
  aucplot.xlabel,
  aucplot.ylabel,
  aucplot.xaxislab,
  aucplot.yaxislab,
  aucplot.xtck,
  aucplot.ytck,
  aucplot.ylimits, 
  aucplot.text,
  aucplot.coord.text,
  aucplot.info,
  aucplot.coord.info,
  aucplot.output,
  aucplot.variable,
  aucplot.line
  ){
  BoutrosLab.plotting.general::create.scatterplot(
    formula = aucplot.variable,
    data = aucplot.data,
    #main = aucplot.main,
    #main.cex = 1,
    xlab.label = aucplot.xlabel,
    ylab.label = aucplot.ylabel,
    xaxis.lab = aucplot.xaxislab,
    xaxis.cex = 1.0,
    yaxis.cex =	1.0,
    yaxis.lab =	aucplot.yaxislab,
    xaxis.tck = aucplot.xtck,
    yaxis.tck = aucplot.ytck,
    xlimits = c(0,1),
    ylimits = aucplot.ylimits,
    xat = seq(0, 1.1, 0.2),
    yat = seq(0, 1.1, 0.2),
    add.xyline = aucplot.line,
    add.points = TRUE,
    xyline.col = 'black',
    type = c('p', 'l'),
    key = list(
      text = list(
        lab = aucplot.text,
        cex = 1,
        font = 'bold'
        ),
        x = aucplot.coord.text[1],
        y = aucplot.coord.text[2],
        padding.text = 3,
        corner = c(0,1)
      ),
    text.fontface = 'bold',
    add.text = TRUE,
    text.x = aucplot.coord.info[[1]],
    text.y = aucplot.coord.info[[2]],
    text.cex = c(1, 1),
    text.label = aucplot.info,
    resolution = 300
    );
  }

#biodb.bak <- biodb
#biodb <- biodb[which(biodb$Race==1),]

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

  roc1 <- roc(biomarker.test$actual, biomarker.test$predicted);
  print(paste('Hello ', ci.auc(roc1,  conf.level=0.95, boot.n=4000, boot.stratified=TRUE, reuse.auc=TRUE), sep=''));
  
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
  ci.auroc.lower[index] <- ci(roc(biomarker.test$actual, biomarker.test$predicted))[1];
  auroc[index] <- ci(roc(biomarker.test$actual, biomarker.test$predicted))[2];
  ci.auroc.upper[index] <- ci(roc(biomarker.test$actual, biomarker.test$predicted))[3];
  
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
#sink();

plot <- TRUE;
if (TRUE == plot){
  print('Plot the AUPRC data');
  # Plot the AUPRC values in a Forest Plot:
  auprc.data <- data.frame(
    feature = variables.set,
    AUPRC = auprc,
    upper = ci.auprc.upper,
    lower = ci.auprc.lower
    );
  
  auprc_plot.file <- BoutrosLab.utilities::generate.filename(
    Sys.Date(),
    'aupr_plot',
    'png'
    )
  
  # AUROC Plot:
  auroc.data <- data.frame(
    feature = variables.set,
    AUROC = auroc,
    upper = ci.auroc.upper,
    lower = ci.auroc.lower
    );

  auroc_plot.file <- BoutrosLab.utilities::generate.filename(
    Sys.Date(),
    'auroc_plot',
    'png'
    )

  simple.data <- as.data.frame(
    cbind(
      #seq(1, length(auprc.data$AUPRC),1),
      as.character(auprc.data$feature),
      auprc.data$AUPRC,
      auprc.data$upper,
      auprc.data$lower,
      auroc.data$AUROC,
      auroc.data$upper,
      auroc.data$lower
      )
      );
  
  colnames(simple.data) <- c(
    #'index',
    'feature',
    'AUPRC',
    'upper.auprc',
    'lower.auprc',
    'AUROC',
    'upper.auroc',
    'lower.auroc'
     );
  simple.data <- simple.data[with(simple.data, order(AUPRC, AUROC, decreasing = TRUE)),];
  rownames(simple.data) <- NULL;
  simple.data[,2:ncol(simple.data)] <- sapply( simple.data[,2:ncol(simple.data)] , as.character )
  simple.data[,2:ncol(simple.data)] <- sapply( simple.data[,2:ncol(simple.data)] , as.numeric )
  
  auprc.plot <- forest.plot(
  #forest.plot(                                                                                                                                                       
    forest.data = simple.data[,c('feature', 'AUPRC', 'lower.auprc', 'upper.auprc')],
    forest.formula = reorder(feature, AUPRC) ~ lower.auprc + upper.auprc,
    forest.xlabel = 'AUPRC',
    forest.ylabel = NULL, #'Test',  
    forest.yaxis.cex = 2,
    forest.yaxis.tck = c(1,0),
    forest.center = simple.data$AUPRC
    #forest.output = auprc_plot.file                                                                                                                                  
    );
  
  auroc.plot <- forest.plot(
  #forest.plot(
    forest.data = simple.data[,c('feature', 'AUPRC', 'AUROC', 'lower.auroc', 'upper.auroc')],
    forest.formula = reorder(feature, AUPRC) ~ lower.auroc + upper.auroc,
    forest.xlabel = 'AUROC',
    forest.ylabel = NULL, #'Test',
    forest.yaxis.cex = 0,
    forest.yaxis.tck = c(1,0),
    forest.center = simple.data$AUROC
    # forest.output = auroc_plot.file
    );
  
  # Plot the AUPRC values in a Forest Plot:           
  #simple.data <- data.frame(
  #  feature = variables.set,
  #  AUROC = auroc,
  #  upper = ci.auroc.upper,
  #  lower = ci.auroc.lower
  #  );

  auprc_plot.file <- BoutrosLab.utilities::generate.filename(
    Sys.Date(),
    'auroc_plot',
    'png'
    )

  # Plotting the Multipanel Plot for the AUPRC and the AUROC:
  BoutrosLab.plotting.general::create.multipanelplot(
    plot.objects = list(auprc.plot, auroc.plot),
    height = 15,
    width = 20,
    plot.objects.heights = c(0.3),
    plot.objects.widths = c(0.1, 0.1),
    layout.height = 1,
    layout.width = 2,
    x.spacing = 0.1,
    left.padding = 25,
    xlab.axis.padding = 5,
    filename = 'auprc_auroc.png',
    resolution = 300
    );
  
  # Plotting the Multipanel Plot for all the PRC:    
  multipanel.file <- BoutrosLab.utilities::generate.filename(
    Sys.Date(),
    'multipanelplot_prc',
    'png'
    )
  BoutrosLab.plotting.general::create.multipanelplot(
    plot.objects = prcplots.set,
    height = 15,
    width = 15,
    plot.objects.heights = c(0.75, 0.75, 0.75, 0.75, 0.75),
    plot.objects.widths = c(0.75, 0.75, 0.75, 0.75, 0.75),
    layout.height = 5,
    layout.width = 5,
    xlab.label = 'Recall',
    ylab.label = 'Precision',
    xlab.cex = 2,
    ylab.cex = 2,
    xlab.axis.padding = 0,
    ylab.axis.padding = 0,
    x.spacing = -0.5,
    y.spacing = -0.5,
    filename = multipanel.file,
    resolution = 300
    );
  
  # Plotting the Multipanel Plot for all the PRC: 
  multipanel.file <- BoutrosLab.utilities::generate.filename(
    Sys.Date(),
    'multipanelploprt_roc',
    'png'
    )
  BoutrosLab.plotting.general::create.multipanelplot(
    plot.objects = rocplots.set,
    height = 15,
    width = 15,
    plot.objects.heights = c(0.75, 0.75, 0.75, 0.75, 0.75),
    plot.objects.widths = c(0.75, 0.75, 0.75, 0.75, 0.75),
    layout.height = 5,
    layout.width = 5,
    xlab.label = '1 - Specificity',
    ylab.label = 'Sensitivity',
    xlab.cex = 2,
    ylab.cex = 2,
    xlab.axis.padding = 0,
    ylab.axis.padding = 0,
    x.spacing = -0.5,
    y.spacing = -0.5,
    filename = multipanel.file,
    resolution = 300
    );
  }

### WRITE SESSION PROFILE TO FILE #####################     
save.session.profile(
  BoutrosLab.utilities::generate.filename(
    Sys.Date(),
    'AUPRC.info',
    'txt'
    )
  );
