
# Project : Prostate Cancer Active Surveillance

### confusion.matrix #################################
# Description: Computes the confusion matrix
#
# Variables
# input.data  two columns data frame (actual and predicted) 
# thresholds.set Thresholds values to compute the confusion matrix
#
# Example:
#
#  result <- confusion.matrix(
#    input.data = my.data,
#    thresholds.set = c(10, 20)
#    );
#
######################################################
confusion.matrix <- function(
  input.data,
  thresholds.set
  ){    

  common.size <- rep(NA, length(thresholds.set));
  merit.data <- data.frame(
    tp = common.size,
    fp = common.size,
    fn = common.size,
    tn = common.size,
    f1.score = common.size,
    threshold = common.size,
    tp.positive.list = common.size,
    tp.positive.patients = common.size,
    fp.positive.list = common.size,
    fp.positive.patients = common.size
  );
    
  index <- 1;
  for (threshold in thresholds.set){
    data <- input.data;
    data$predicted[threshold <= data$predicted] <- 1;
    data$predicted[1 != data$predicted] <- 0;
    
    if (0 < length(which(1 == data$predicted))){
    
      # Predicted data with values of 1 when compared to actual data, two outcomes are
      # possible match with the actual value (True Positive) or miss it (False Positive)
      positive.index  <- which(1 == data$predicted);
      #print(positive.index)
      actual.positive <- which(1 == data$actual);
      indices <- actual.positive[which(actual.positive %in% positive.index)];
      #print(indices)
        
      # We sort them by scores. It will be helpful later:
      temp <- data.frame(
        positive.list = round(input.data$predicted[c(indices)]),
        positive.patients = as.character(input.data$ID[indices])
        );  
      temp <- temp[order(temp$positive.list),];

      merit.data$tp.positive.list[index] <- paste(temp$positive.list, collapse = ',');
      merit.data$tp.positive.patients[index] <- paste(temp$positive.patients , collapse = ',');  

      # If indices for positive in predicted data are present in actual data, they are
      # TRUE POSITIVE:
      merit.data$tp[index] <- length(which(actual.positive %in% positive.index));
      # The remaining is FALSE POSITIVE:
      merit.data$fp[index] <- length(positive.index) - merit.data$tp[index];

      # Getting the False Positive Members:
      fp.indices <- positive.index[!positive.index %in% indices];
      temp <- data.frame(
        positive.list = round(input.data$predicted[c(fp.indices)]),
        positive.patients = as.character(input.data$ID[fp.indices])
        );
      temp <- temp[order(temp$positive.list),];
      merit.data$fp.positive.list[index] <- paste(temp$positive.list, collapse = ',');
      merit.data$fp.positive.patients[index] <- paste(temp$positive.patients , collapse = ',');
      }
      else {
        # This is the case when NO POSITIVE is present in predicted data
        tp <- 0;
        fp <- 0;
      }
    
    if (0 < length(which(0 == data$predicted))){
    
      # Similar as above. Two outcomes are possible: Predicted value of 0 that hits the actual
      # data are TRUE NEGATIVE. Otherwise, FALSE NEGATIVE.
      negative.index <- which(0 == data$predicted);
      actual.negative <- which(0 == data$actual);    
      merit.data$tn[index] <- length(which(actual.negative %in% negative.index));    
      merit.data$fn[index] <- length(negative.index) - merit.data$tn[index]; 
      }
      else {
        # This is the case when NO NEGATIVE is present in predicted data
        merit.data$tn[index] <- 0;
        merit.data$fn[index] <- 0; 
      }

      merit.data$threshold[index] <- threshold;
    
    index <- index + 1; 
    }
    
    merit.data$f1.score <- 2*merit.data$tp/(2*merit.data$tp + merit.data$fp + merit.data$fn);
    return(merit.data);    
  }

### WRITE SESSION PROFILE TO FILE #####################     
save.session.profile(
  BoutrosLab.utilities::generate.filename(
    Sys.Date(),
    'confusion.matrix',
    'txt'
    )
  );
