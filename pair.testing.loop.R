
### pair.testing.loop.R ###################################
# The following R code take a set tests with their
# respective thresholds and perfoms a serial test

# Cleaning Environment:
rm(list=ls());

# Loading Data:
source('load.data.R');

# Loading Confusion Matrix Function:
source('confusion.matrix.R');

### generate.legend ###################################
# Description:                                                 
# Function to generate scale used in hexbin plot
#
# Input variables:                                             
# legend.colours  color scale to be used               
# legend.labels   give label to each color
# legend.title    title of the scale                           
# legend.output   output filename                                   
# legend.boolean  continous (TRUE or FALSE)
#
# Example:       
# sample.names.legend  <- generate.legend( LinGray(11, beg=1, end=92),
#   rev(c('1', '16', '31', '46', '61', '76', '90', '105', '120', '135', '150')),
#   expression(bold('Counts')),
#   'FALSE'
#   );
#######################################################                                    
generate.legend <- function(
  legend.colours,
  legend.labels,
  legend.title,
  legend.boolean){
  list(
    colours = legend.colours,
    labels = legend.labels,
    title = legend.title,
    height = 15,
    width = 5,
    size = 4,
    between.row = 3,
    between.col = 3,
    #tck = 2,
    continuous = legend.boolean,
    title.cex = 5,
    axis.label = 5,
    #label.cex = 20.0,
    border = 'black'
    );
  }

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

# If we want to make a fair comparison between pairs,
# all of them must have the same set of subjects
# participating in all of the tests :
data.set <- data.set[complete.cases(data.set), ];

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
  'PSADensity',
  'p2PSA',
  'PercentFreePSA',
  'PHI',
  'PHIDensity',
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
  'PSA Density',
  'p2PSA',
  'Percent Free PSA',
  'PHI',
  'PHI Density',
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

# Optimal Set:
optimal.set.file <- '2020-06-28_AS_optimal.set.RData';
load(optimal.set.file);
biomarker.tests <- optimal.set$test;
biomarker.thresholds <- optimal.set$threshold;

# Making a data frame to pair the tests with their respective thresholds:  
common.size <- rep(NA, length(biomarker.tests)*(length(biomarker.tests)));
pair.tests.names <- data.frame(
  testA = common.size,
  testB = common.size,
  testA.threshold = common.size,
  testB.threshold = common.size
  );

index <- 1;
for (i in 1:length(biomarker.tests)){
  for (j in 1:(length(biomarker.tests))){
    #if (i != j){
      pair.tests.names$testA[index] <- biomarker.tests[i];
      pair.tests.names$testB[index] <- biomarker.tests[j];
      pair.tests.names$testA.threshold[index] <- biomarker.thresholds[i];
      pair.tests.names$testB.threshold[index] <- biomarker.thresholds[j];
      index = index + 1;
    #}
    }
  }

common.size = rep(NA, length(variables.set));
pair.tests.sensitivity <- data.frame(
  PCA3 = common.size,
  T2ERG = common.size,
  MiPSCancerRisk = common.size,
  MiPSHGCancerRisk = common.size,
  PSAHyb = common.size,
  freePSA = common.size,
  p2PSA = common.size,
  PercentFreePSA = common.size,
  PHI = common.size,
  GeneticAncestry = common.size,
  GeneticRiskScore = common.size,
  GeneticRiskCategory = common.size,
  GlobalScreeningArray = common.size,
  GSAPositives = common.size,
  BRCAMutation = common.size,
  Mutation1 = common.size,
  Mutation2 = common.size,
  TNFaAverage = common.size,
  ADCnormalSignal = common.size,
  ADClesionSignal = common.size,
  RSIlesionSignal = common.size,
  RSInormalSignal = common.size,
  SOCPSA = common.size
  );
rownames(pair.tests.sensitivity) <- variables.set;

pair.tests.specificity <- data.frame(
  PCA3 = common.size,
  T2ERG = common.size,
  MiPSCancerRisk = common.size,
  MiPSHGCancerRisk = common.size,
  PSAHyb = common.size,
  freePSA = common.size,
  p2PSA = common.size,
  PercentFreePSA = common.size,
  PHI = common.size,
  GeneticAncestry = common.size,
  GeneticRiskScore = common.size,
  GeneticRiskCategory = common.size,
  GlobalScreeningArray = common.size,
  GSAPositives = common.size,
  BRCAMutation = common.size,
  Mutation1 = common.size,
  Mutation2 = common.size,
  TNFaAverage = common.size,
  ADCnormalSignal = common.size,
  ADClesionSignal = common.size,
  RSIlesionSignal = common.size,
  RSInormalSignal = common.size,
  SOCPSA = common.size
  );
rownames(pair.tests.specificity) <- variables.set;

# NOTE: Since we are comparing two tests, we need to make sure
# subjects did both tests

# Testing Purposes:
# pair.tests.names <- pair.tests.names[23,];

for (i in 1:dim(pair.tests.names)[1]){
  testA <- pair.tests.names$testA[i];
  testB <- pair.tests.names$testB[i];
  
  print(paste('Processing ', testA, ' - ', testB, sep = ''));
  testA.threshold <- pair.tests.names$testA.threshold[i];
  testB.threshold <- pair.tests.names$testB.threshold[i];
  
  working.data <- data.set[,c('ID', 'BiopsyUpgraded', testA, testB)] 
  working.data <- working.data[complete.cases(working.data),];

  testA.set <- working.data[,c('ID', 'BiopsyUpgraded', testA)];
  colnames(testA.set) <- c('ID', 'actual', 'predicted');
  testA.set$predicted <- round(as.numeric(testA.set$predicted));

  testB.set <- working.data[,c('ID', 'BiopsyUpgraded', testB)];
  colnames(testB.set) <- c('ID', 'actual', 'predicted');
  testB.set$predicted <- round(as.numeric(testB.set$predicted));

  # Computing Confusion Matrix for Tests A and B:
  testA.matrix = confusion.matrix(testA.set, testA.threshold);
  testB.matrix = confusion.matrix(testB.set, testB.threshold);  

  if (!is.na(testA.matrix$tp.positive.patients)){
    tp.patients <- unlist(strsplit(testA.matrix$tp.positive.patients, ','));
    tp.indices <- which(testA.set$ID %in% tp.patients);
    testB.matrix.tpA = confusion.matrix(testA.set[tp.indices,], testB.threshold);
    fp.patients <- unlist(strsplit(testA.matrix$fp.positive.patients, ','));
    fp.indices  <- which(testA.set$ID %in% fp.patients);

    testB.matrix.fpA = confusion.matrix(testA.set[fp.indices,], testB.threshold);
    }

  if (is.na(testA.matrix$tp.positive.patients)){
    testA.matrix$tp <- 0;
    testA.matrix$fp <- 0;
    }
  
  # Computing Sensitivity, Specificity for test A and test B:
  sensitivityA <- testA.matrix$tp/(testA.matrix$tp + testA.matrix$fn);
  specificityA <- testA.matrix$tn/(testA.matrix$tn + testA.matrix$fp);

  sensitivityB <- testB.matrix$tp/(testB.matrix$tp + testB.matrix$fn);
  specificityB <- testB.matrix$tn/(testB.matrix$tn + testB.matrix$fp);

  # Computing Sensitivity, Specificity for tree:
  if (is.na(testB.matrix.tpA$tn) || is.na(testB.matrix.tpA$fp) ){
    sensitivityB.tree  <- 0;
    }
    else {
      sensitivityB.tree  <- testB.matrix.tpA$tp/(testB.matrix.tpA$tp + testB.matrix.tpA$fn);
    }

  if (is.na(testB.matrix.fpA$tn) || is.na(testB.matrix.fpA$fp) ){
    specificityB.tree  <- 0;
    }
    else {
      specificityB.tree  <- testB.matrix.fpA$tn/(testB.matrix.fpA$tn + testB.matrix.fpA$fp);
    }
   
  print(specificityA);
  print(specificityB.tree);
  
  overall.sensitivity <- sensitivityA*sensitivityB.tree;
  overall.specificity <- specificityA + (1-specificityA)*specificityB.tree;

  #overall.accuracy <-
  #(testB.matrix.tpA$tp + testB.matrix.tpA$tn)/(testB.matrix.tpA$tp + testB.matrix.tpA$tn+testB.matrix.tpA$fn+testA.matrix$fn)
    
  #print(sensitivityA)
  #print(specificityA)

  #print(sensitivityB)
  #print(specificityB)

  #print(overall.sensitivity)
  #print(overall.specificity)

  #pair.tests.sensitivity[];  
  #pair.tests.sensitivity[testA, testA] <- sensitivityA;
  pair.tests.sensitivity[testB, testA] <- overall.sensitivity
  pair.tests.specificity[testB, testA] <- overall.specificity

  }

# Scaling the colors:              
key.min = 0;
key.max = 1;
key.colour.interval.num = 100;
key.scale = c(
  seq(
    key.min,
    0,
    -key.min/key.colour.interval.num
    ),
  seq(
    0,
    key.max,
    key.max/key.colour.interval.num
    )
  );
key.scale = unique(key.scale);

# Plotting Heatmap:
chr.cov.colours <- c(
  rep('dodgerblue',4),
  rep('gold',2),
  rep('firebrick3', 9),
  rep('darkgreen', 8));

top.covariate <- list(
  rect = list(
    col = 'white',
    fill = chr.cov.colours,
    lwd = 1.5
    )
  );

sample.covariate <- list(
  rect = list(
    col = 'white',
    fill =  chr.cov.colours,
    lwd = 1.5
    )
  );

# Reordering by test type:
pair.tests.sensitivity <- pair.tests.sensitivity[,c(variables.by.type)];
pair.tests.sensitivity <- pair.tests.sensitivity[c(variables.by.type),];

# Reordering by test type: 
pair.tests.specificity <- pair.tests.specificity[,c(variables.by.type)];
pair.tests.specificity <- pair.tests.specificity[c(variables.by.type),];

specificity.plot <- BoutrosLab.plotting.general::create.heatmap(
  x = pair.tests.specificity,
  main = 'Specificity',
  xaxis.lab = real.names.set,
  xaxis.cex = 1.75,
  yaxis.lab = NULL,
  yaxis.cex = 1.75,
  #ylab.label = 'Test A',
  #ylab.axis.padding = 2,
  #filename = "specificity.png",
  colourkey.cex = 0.75,
  colourkey.labels.at = seq(0, 1, 0.2),
  at = key.scale,
  covariates = sample.covariate,
  covariates.top = top.covariate,  
  #xaxis.col = c(rep('dodgerblue',4), rep('gold',2), rep('firebrick3', 9), rep('darkgreen', 8)),
  #yaxis.col = c(rep('dodgerblue',4), rep('gold',2), rep('firebrick3', 9), rep('darkgreen', 8)),
  clustering.method = 'none',
  print.colour.key = FALSE,
  resolution = 300
  );

sensitivity.plot <- BoutrosLab.plotting.general::create.heatmap(
  x = pair.tests.sensitivity,
  main = 'Sensitivity',
  xaxis.lab = real.names.set,
  xaxis.cex = 1.75,
  yaxis.lab = real.names.set,
  yaxis.cex = 1.75,
  #filename = "sensitivity.png",
  colourkey.cex = 0.75,
  colourkey.labels.at = seq(0, 1, 0.2),
  at = key.scale,
  covariates.top = top.covariate,  
  #xaxis.col = c(rep('dodgerblue',4), rep('gold',2), rep('firebrick3', 9), rep('darkgreen', 8)),
  #yaxis.col = c(rep('dodgerblue',4), rep('gold',2), rep('firebrick3', 9), rep('darkgreen', 8)),
  clustering.method = 'none',
  print.colour.key = FALSE,
  resolution = 300
  );

test.legend  <- generate.legend(        
  legend.colours = c('dodgerblue','gold','firebrick3', 'darkgreen'),  
  legend.labels = c('Imaging', 'Urine', 'Blood/Urine', 'Genetics'),         
  legend.title = expression(bold('Test')),                    
  legend.boolean = 'FALSE'              
  );                

heatmap.legend  <- generate.legend(
  legend.colours = c('red', 'white', 'blue'),
  legend.labels = c('1', '0'),
  legend.title = expression(bold('Intensity')),
  legend.boolean = 'TRUE'
  );

legendG <- legend.grob(
  list(
    legend = test.legend,
    legend = heatmap.legend
    ),
    label.cex = 2,
    title.cex = 2
  );

BoutrosLab.plotting.general::create.multipanelplot(
  plot.objects = list(sensitivity.plot, specificity.plot),
  height = 15,
  width = 30,
  plot.objects.heights = c(0.3),
  plot.objects.widths = c(0.65, 0.5),
  layout.height = 1,
  layout.width = 2,
  x.spacing = 0.05,
  legend = list(right = list(fun = legendG)),
  left.padding = 5,
  xlab.axis.padding = 6,
  filename = 'pairs_sensitivity_specificity.png',
  resolution = 300
  );

### WRITE SESSION PROFILE TO FILE #####################
save.session.profile(
  BoutrosLab.utilities::generate.filename(
    Sys.Date(),
    'pair.testing.loop',
    'txt'
    )
  );

