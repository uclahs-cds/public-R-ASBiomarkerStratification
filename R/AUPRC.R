### AUPRC.R #################################
# Perform AUPRC and PR plotting and analysis
# using BiopsyUpgraded and Biomarkers tests

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
#' @export
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
#' @export
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
#' @export
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
#' @export
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
#' @export
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

#' @export
plot.AUPRC <- function() {
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
