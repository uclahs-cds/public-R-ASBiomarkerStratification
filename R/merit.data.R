
# Project : Prostate Cancer Active Surveillance

#' @export
compute.merit.data <- function(
  predicted.data,
  target,
  biomarker.test,
  original.data
  ){

  # Rescaling the predicted signal:
  numerator <- predicted.data[,target] - min(predicted.data[,target]);
  denominator <- max(predicted.data[,target]) - min(predicted.data[,target]);
  predicted.data <- 100*numerator/denominator;

  print(numerator)
  print(denominator)

  lower.threshold <- min(predicted.data[,target]) - step;
  upper.threshold <- max(predicted.data[,target]) + step;
  print(paste(
    ' Lower Threshold : ',
    lower.threshold,
    sep = ''
    )
    );
  print(paste(
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
    positive.list = common.size
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

  merit.data$threshold <- thresholds.set;
  merit.data$scores <- thresholds.set*denominator/100 + min(original.data);
  print(min(original.data))

  print(merit.data$scores)
  print(original.data)


  for (j in 1:length(thresholds.set)) {
    threshold <- thresholds.set[j];
    data <- biomarker.test;
    data[,target][threshold <= data[,target]] <- 1;
    data[,target][1 != data[,target]] <- 0;

    #print('data')
    #print(data);

    if (0 < length(which(1 == data[,target]))){

      # Predicted data with values of 1 when compared to actual data, two outcomes are
      # possible match with the actual value (True Positive) or miss it (False Positive)
      positive.index <- which(1 == data[,target])
      actual.positive <- which(1 == data[,'actual']);

      #print(positive.index);
      #print(actual.positive);

      indices <- actual.positive[which(actual.positive %in% positive.index)]

      #print(paste('Hello ',  indices));

      #print(original.data)
      merit.data$positive.list[j] <- paste(original.data[indices], collapse=',');

      #print('HERE')

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
#
#    if (0 < length(which(0 == data$predicted))){
#
#      # Similar as above. Two outcomes are possible: Predicted value of 0 that hits the actual
#      # data are TRUE NEGATIVE. Otherwise, FALSE NEGATIVE.
#      negative.index <- which(0 == data[,target])
#      actual.negative <- which(0 == data[,actual]);
#      tn[j] <- length(which(actual.negative %in% negative.index));
#      fn[j] <- length(negative.index) - tn[j];
#      }
#      else {
#        # This is the case when NO NEGATIVE is present in predicted data
#        tn[j] <- 0;
#        fn[j] <- 0;
#      }
#
    }
#
#  # Generating Merit Data and save it in a R object:
#  merit.data$tp <- tp;
#  merit.data$tn <- tn;
#  merit.data$fp <- fp;
#  merit.data$fn <- fn;
#  merit.data$specificity <- tn/(tn + fp);
#  merit.data$fpr <- fp/(fp + tn);
#  merit.data$precisions <- tp/(tp + fp);
#  merit.data$recalls <- tp/(tp + fn);
#  merit.data$accuracies <- (tp + tn)/(tp + fp + tn + fn);
#  merit.data$f1.scores <- 2*merit.data$recalls*merit.data$precisions/(merit.data$recalls + merit.data$precisions);
#  merit.data[which(0 == merit.data$recalls),]$precisions <- 1;
#  merit.data[which(0 == merit.data$recalls),]$f1.scores <- 0;

  return(merit.data);
  }
