library(ProstateCancer.ASBiomarkerSynergy);

# Project : Prostate Cancer Active Surveillance

### cv.test.R #######################################
# Perform AUPRC and PR plotting and analysis
# using BiopsyUpgraded and Biomarkers tests

# Load Data:
data <- default.load.data();
attach(data);

### loocv ############################################
# Description : Leave one out cross-validation
#
# Variables:
# tests.names tests names
# input.data  two-columns data-frame
# thresholds.set  set of values to be used to define
#                 true and false
#
# Example:
#
# outcome <- loocv(
#   tests.names = tests,
#   input.data = my.data,
#   thresholds.set = c(10, 20)
#   );
#
######################################################
loocv.pair <- function(
  tests.names,
  input.data,
  thresholds.set
  ){

  thresholdsA <- unique(thresholds.set$thresholdsA);
  thresholdsB <- unique(thresholds.set$thresholdsB);

  size <- rep(NA, length(thresholdsA)*length(thresholdsB));
  outcome <- data.frame(
    row.number = size,
    testA.name = size,
    testB.name = size,
    thresholdA = size,
    thresholdB = size,
    sensitivityAB = size,
    specificityAB = size,
    precisionAB = size,
    accuracyAB = size,
    f1.scoreAB = size,
    classification = size
    );
  r <- 1;

  print(thresholdsA)
  print(thresholdsB)

  # Training the model (test A, test B):
  size <- rep(NA, dim(input.data)[1]*length(thresholdsA)*length(thresholdsB));
  tests.paired <- data.frame(
    row.number.rm = size,
    testA.name = size,
    testB.name = size,
    thresholdA = size,
    thresholdB = size,
    sensitivityAB = size,
    specificityAB = size,
    precisionAB = size,
    accuracyAB = size,
    f1.scoreAB = size
    );

  m <- 1;
  # Loop through each score of the test:
  for (i in 1:nrow(input.data)){
    testA.message <- paste('Removing row ', i , ' from data', sep = '');
    print(testA.message);
    # Removing one point:
    trainA.data <- input.data[-i,];

    # Data for test A:
    trainA.data <- input.data[-i, c('ID', 'actual', 'predictedA')];
    colnames(trainA.data) <- c('ID', 'actual', 'predicted');
    #print(trainA.data)

    # Data for test B:
    trainB.data <- input.data[-i, c('ID', 'actual', 'predictedB')];
    colnames(trainB.data) <- c('ID', 'actual', 'predicted');
    #print(trainB.data)

    # Loop through Test A thrseholds:
    for (j in 1:length(thresholdsA)){
      # Confusion Matrix for test A:
      testA.df <- confusion.matrix(trainA.data, thresholdsA[j]);
      #print(testA.df)
      #print(paste('TN Test A ', testA.df$fp, sep = ''));
      #print(paste('FN Test A ', testA.df$tn, sep = ''));

      if (!is.na(testA.df$tp.positive.patients)){
        tp.patients <- unlist(strsplit(testA.df$tp.positive.patients, ','));
        tp.indices <- which(trainA.data$ID %in% tp.patients);
        }

      if (!is.na(testA.df$fp.positive.patients)){
        fp.patients <- unlist(strsplit(testA.df$fp.positive.patients, ','));
        fp.indices  <- which(trainA.data$ID %in% fp.patients);
        };

      # Computing Sensitivity for test A:
      if (0 == testA.df$tp && 0 == testA.df$fn){
        sensitivityA <- 0;
        }
        else if(is.na(testA.df$tp) && is.na(testA.df$fn)){
          sensitivityA <- 0;
        }
        else {
          sensitivityA <- testA.df$tp/(testA.df$tp + testA.df$fn);
        }

      #print(paste('TN ' , testA.df$tn, sep = ''));
      #print(paste('FP ', testA.df$fp, sep = ''));

      # Computing Specificity for test A:
      if (0 == testA.df$tn && 0 == testA.df$fp){
        specificityA <- 0;
        }
        else if (!is.na(testA.df$tn) && is.na(testA.df$fp)){
          specificityA <- 1;
        } else {
        specificityA <- testA.df$tn/(testA.df$tn + testA.df$fp);
        }

      if (!is.na(testA.df$tp.positive.patients)  && !is.na(testA.df$fp.positive.patients)){

        # Loop through Test B:
        for (k in 1:length(thresholdsB)){

          #print(paste('i ',i, ' j ', j,  ' k', k, sep = ''))

          tests.paired$row.number.rm[m] = i;
          tests.paired$testA.name[m] = tests.names[1];
          tests.paired$testB.name[m] = tests.names[2];
          tests.paired$thresholdA[m] <- thresholdsA[j];
          tests.paired$thresholdB[m] <- thresholdsB[k];
          #on.screen <- paste('Point ', i ,' ThresholdA ', thresholdsA[j], ' ThresholdB ', thresholdsB[k]);
          #print(on.screen);
          testB.df.tpA = confusion.matrix(trainB.data[tp.indices,], thresholdsB[k]);
          testB.df.fpA = confusion.matrix(trainB.data[fp.indices,], thresholdsB[k]);

          #print(paste('TP Tree ', testB.df.tpA$tp, sep = ''))
          #print(paste('FN Tree ', testB.df.tpA$fn, sep = ''))

          # Computing Sensitivity, Specificity for test B in each branch of the tree:
          if (is.na(testB.df.tpA$tn) || is.na(testB.df.tpA$fp) ){
            sensitivityB.tree  <- 0;
            }
            else {
              if (is.na(testB.df.tpA$fn) && !is.na(testB.df.tpA$tp)){
                sensitivityB.tree  <- 1;
                }
                else {
                  sensitivityB.tree  <- testB.df.tpA$tp/(testB.df.tpA$tp + testB.df.tpA$fn);
                }
            }

          #print(paste('FP Tree ', testB.df.fpA$fp, sep = ' '));
          #print(paste('TN Tree ', testB.df.fpA$tn, sep = ' '));

          if (is.na(testB.df.fpA$tn)){
            if (is.na(testB.df.fpA$fp)  || !is.na(testB.df.fpA$fp)){
              specificityB.tree  <- 0;
            }
            }
            else {
              if (is.na(testB.df.fpA$fp) && !is.na(testB.df.fpA$tn)){
                specificityB.tree <- 1;
                }
                else {
                  specificityB.tree  <- testB.df.fpA$tn/(testB.df.fpA$tn + testB.df.fpA$fp);
                }
            }

          # Overall Precision:
          if (is.na(testB.df.tpA$tp)){
            overall.precision <- 0;
            } else {
              if (is.na(testB.df.fpA$fp)){
                overall.precision <- 1;
                }
                else {
                  overall.precision <- testB.df.tpA$tp/(testB.df.tpA$tp + testB.df.fpA$fp)
                }
            }
          # Overall Accuracy:
          total.tn <- testA.df$tn + testB.df.fpA$tn;
          total.fp <- testA.df$fp + testB.df.fpA$fp;
          total.true <- testB.df.tpA$tp + testB.df.fpA$tn;
          total.false <- total.tn + total.fp;
          overall.accuracy <- (total.true)/(total.true + total.false)

          #print(specificityA);
          #print(specificityB.tree);

          overall.sensitivity <- sensitivityA*sensitivityB.tree;
          overall.specificity <- specificityA + (1-specificityA)*specificityB.tree;

          tests.paired$sensitivityAB[m] <- overall.sensitivity;
          tests.paired$specificityAB[m] <- overall.specificity;
          tests.paired$precisionAB[m] <- overall.precision;
          tests.paired$accuracyAB[m] <- overall.accuracy;

          overall.f1score	<- 2*overall.precision*overall.sensitivity/(overall.precision + overall.sensitivity);
          tests.paired$f1.scoreAB[m]  <- overall.f1score;

          m <- m + 1;
          #print(overall.sensitivity);
          #print(overall.specificity);

        } # The thresholds loop for test B ends here

        }

        temp <- tests.paired[which(tests.paired$row.number.rm == i), ];
        temp <- temp[which(temp$thresholdA == thresholdsA[j]), ];
        #print(temp)

        # Best Model testA-testB based on F1.score::
        best.f1.score <- max(temp$f1.scoreAB, na.rm = TRUE);
        index <- which(best.f1.score == temp$f1.scoreAB);
        # Since temp dataframe will increase, index also will increase
        # according to the thresholds in A.
        index <- index[1];
        #print(temp[index,])

        # The Optimal Thresholds for test A and B:
        optimalA <- temp[index,]$thresholdA;
        optimalB <- temp[index,]$thresholdB;

        #print(optimalA)
        #print(optimalB)


        #print(r)

        outcome$row.number[r] <- i;
        outcome$testA.name[r] <- temp[index,]$testA.name;
        outcome$testB.name[r] <- temp[index,]$testB.name;
        outcome$thresholdA[r] <- optimalA;
        outcome$thresholdB[r] <- optimalB;
        outcome$sensitivityAB[r] <- temp[index,]$sensitivityAB;
        outcome$specificityAB[r] <- temp[index,]$specificityAB;
        outcome$precisionAB[r] <- temp[index,]$precisionAB;
        outcome$accuracyAB[r] <- temp[index,]$accuracyAB;
        outcome$f1.scoreAB[r] <- temp[index,]$f1.scoreAB;

       #print(confusion.matrix(trainB.data[tp.indices,], optimalA));
       test.data <- input.data[ i, c('ID', 'actual', 'predictedA')];
       colnames(test.data) <- c('ID', 'actual', 'predicted');
       #print(test.data)

       # Test A evaluates test data (1 data point):
       testA.one <- confusion.matrix(test.data, optimalA);
       #print('test A one')
       #print(testA.one[,c('tn', 'tp', 'fn', 'fp')]);

       if (0 < length(testA.one) ) {

         if (!is.na(testA.one$tp) &&  1 == testA.one$tp){
           print(paste('TP Point ', i, ' Processing Test  Optimal A : ', optimalA , ' Optimal B ', optimalB, sep = ''));
           tp.test <- unlist(strsplit(testA.one$tp.positive.patients, ','));
           tp.test.indices <- which(trainA.data$ID %in% tp.test);
           testB.one.tpA = confusion.matrix(test.data, optimalB);
           print('Test B evaluates TP from Test A :  testB.one.tpA ');
           print(testB.one.tpA)
           if (0 == testB.one.tpA$tp || is.na(testB.one.tpA$tp)){
             classification  <- 'tp.misclassified';
             }
             else {
               classification <- 'tp.classified';
             }
           outcome$classification[r] <- classification;
           r <- r + 1;
           }

         if (1 == testA.one$tn){
           print(paste('TN Point ', i, ' Processing Test  Optimal A : ', optimalA , ' Optimal B ', optimalB, sep = ''));
           outcome$classification[r] <- 'tp.classifiedA';
           r <- r + 1;
           }

         if (!is.na(testA.one$fp) && (1 == testA.one$fp)){
           print(paste('FP Point ', i, ' Processing Test  Optimal A : ', optimalA , ' Optimal B ', optimalB, sep = ''));
           fp.test <- unlist(strsplit(testA.one$fp.positive.patients, ','));
           fp.test.indices  <- which(trainA.data$ID %in% fp.test);
           testB.one.fpA = confusion.matrix(test.data, optimalB);
           print('Test B evaluates FP from Test A');
           #print(testB.one.fpA$tn);
           if (0 == testB.one.fpA$tn || is.na(testB.one.fpA$tn) ){
             classification  <- 'tn.misclassified';
             }
             else {
               classification <- 'tn.classified';
             }
           outcome$classification[r] <- classification;
           r <- r + 1;
           }
         }












    } # The thresholds loop for test A ends here

   }  # The loop to remove one point ends here

   return(outcome)

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

# Testing Purposes:
variables.set <- c('PCA3', 'T2ERG');
#variables.set <- c('PCA3');

optimal.set <- data.frame(
  test = rep(NA, length(variables.set)),
  threshold = rep(NA, length(variables.set))
  );

# Making a data frame to pair the tests with their respective thresholds:
common.size <- rep(NA, length(variables.set)*(length(variables.set) -1));
pair.tests.names <- data.frame(
  testA = common.size,
  testB = common.size
  );

index <- 1;
for (i in 1:length(variables.set)){
  for (j in 1:(length(variables.set))){
    if (j!=i) {
      pair.tests.names$testA[index] <- variables.set[i];
      pair.tests.names$testB[index] <- variables.set[j];
      index = index + 1;
      }
    }
  }

count <- 1;
for (i in 1:1){ #dim(pair.tests.names)[1]){
  testA <- pair.tests.names$testA[i];
  testB <- pair.tests.names$testB[i];

  print(paste('Processing Pair ', testA , '-' , testB,  sep = ''));

  # Getting the data we need:
  select <- c('ID', 'BiopsyUpgraded',  testA, testB);
  biomarker.test <- data.set[,c(select)];
  colnames(biomarker.test) <- c('ID', 'actual', 'predictedA', 'predictedB');
  # Some rows has NA or incomplete data:
  biomarker.test <- biomarker.test[complete.cases(biomarker.test),];

  # Test A Thresholds:
  thresholdsA.set <- sort(biomarker.test$predictedA);
  # Test A Thresholds:
  thresholdsB.set <- sort(biomarker.test$predictedB);

  # Put all Thresholds together:
  thresholdsAB <- data.frame(
    thresholdsA = thresholdsA.set,
    thresholdsB = thresholdsB.set
    );

  # Leave-One-Out Cross-Validation for pair (TestA, TestB)
  outcome.cv <- loocv.pair(
    tests.names = c(testA, testB),
    input.data = biomarker.test,
    thresholds.set = thresholdsAB
    );

  outcome.cv.filename <- BoutrosLab.utilities::generate.filename(
    project.stem = 'AS',
    file.core = paste('outcome.cv' , testA, testB, sep = '_'),
    extension = 'RData',
    file.date = Sys.Date()
    );

  save(
    outcome.cv,
    file = outcome.cv.filename
    );
  # Save the thresholds. In case there are more than 1 threshold,
  # choose the lowest:
  #optimal.set$test[count] <- variable;
  #optimal.set$threshold[count] <- min(unique(outcome$threshold), na.rm = TRUE);
  #count <- count + 1;

  }



#
# # Saving the object:
# output.filename <- BoutrosLab.utilities::generate.filename(
#  project.stem = 'AS',
#  file.core = 'optimal.set',
#  extension = 'RData',
#  file.date = Sys.Date()
#  );
# save(
#  optimal.set,
#  file = output.filename
#  );
#

detach(data);
