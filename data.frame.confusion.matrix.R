
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
data.frame.confusion.matrix <- function(
  input.data,
  threshold
  ){    

  common.size <- rep(NA, length(threshold));
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
    
  
    data <- input.data;
    data$predicted[threshold <= data$predicted] <- 1;
    data$predicted[1 != data$predicted] <- 0;

    return(data)
}

### WRITE SESSION PROFILE TO FILE #####################     
save.session.profile(
  BoutrosLab.utilities::generate.filename(
    Sys.Date(),
    'data.frame.confusion.matrix',
    'txt'
    )
  );
