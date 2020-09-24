#' Models for Progressed to Treatment
#'
#' @param biodb
#' @param train.control
#' @param target
#' @param metric
#'
#' @return
#' @export
#'
#' @examples
AS.models <- function(
    biodb,
    train.control = NULL,
    target = c('ProgressedToTreatment', 'BiopsyUpgraded', 'Prostatectomy'),
    metric = c('Accuracy', 'AUC', 'F', 'Kappa', 'Precision', 'Recall', 'ROC',
               'Sens', 'Spec'),
    exclude.vars = NULL) {
    target <- match.arg(target);
    metric <- match.arg(metric);

    biokey.variables <- c('Age',
                          'Race',
                          'Ethnicity',
                          'Weight',
                          'Height',
                          'BMI', 'MRIResult', 'MRILesions',
                          'BiopsyResult',
                          'ProstateVolume',
                          # 'Observation',
                          'PCA3',
                          'T2ERG', 'MiPSCancerRisk', 'MiPSHighGradeCancerRisk', 'PSAHyb',
                          'freePSA', 'p2PSA', 'PercentFreePSA', 'PHI', 'GeneticAncestry',
                          # 'PreviousGleason', 'StudyHighestGleason', 'RSIlesionGleason', # Only use the ISUP grade group over Gleason
                          # 'StudyHighestISUP', 'HighestPIRADS',
                          'PreviousISUP',
                          'GeneticRiskScore',
                          'TNFaAverage', 'TNFaSTD', # Tumor necrosis factor?
                          #'GeneticRiskCategory', Don't need since genetic risk score is continuous version of this
                          'GlobalScreeningArray', # This is just an indicator if any of the follow are > 0
                          'GSAPositives', 'BRCAMutation',
                          'Mutation.BRCA1', 'Mutation.BRCA2', 'Mutation.ATM', # 'Mutation.MLH1', 'Mutation.PMS2', Not enough data
                          'RSInormalSignal',
                          'RSIlesionSignal',
                          'ADCnormalSignal', 'ADClesionSignal',
                          # 'RSIlesionPIRADS', 'RSIlesionCancer',  'RSIlesionUpgraded', 'RSIlesionISUP',
                          'PSADensity', 'PHIDensity');

    if('BiopsyUpgraded' == target) {
        exclude.vars <- c(exclude.vars, 'BiopsyResult');
    }

    biokey.variables <- setdiff(biokey.variables, exclude.vars)

    missing.target <- is.na(biodb[, target])
    # Variable we want to predict
    y.target <- biodb[!missing.target, target];

    # Drop the missing target values
    X <- biodb[!missing.target, biokey.variables]

    gbmGrid <-  expand.grid(interaction.depth = c(1, 5, 9),
                            n.trees = (1:30)*50,
                            shrinkage = c(0.001, 0.01, 0.1),
                            n.minobsinnode = 20)

    custom.summary <- function (data, lev = NULL, model = NULL) {
        c(
            defaultSummary(data, lev, model),
            twoClassSummary(data, lev, model),
            prSummary(data, lev, model)
        )
    }

    if(is.null(train.control)) {
        train.control <- trainControl(
            method = 'repeatedcv',
            number = 10,
            repeats = 5,
            classProbs = TRUE,
            summaryFunction = custom.summary
            # https://topepo.github.io/caret/subsampling-for-class-imbalances.html
            # Don't seem to be enough of a class imbalance to make a difference
            # sampling = 'up'
            )
        } else {
            train.control$classProbs <- TRUE
            train.control$summaryFunction <- custom.summary
        }

    y <- y.target
    levels(y) <- c('no', 'yes');

    rpart.fit <- train(
        X,
        y,
        method = 'rpart',
        metric = metric,
        trControl = train.control,
        tuneLength = 30
    )

    c50.grid <- expand.grid(
        .winnow = c(TRUE, FALSE),
        .trials = c(1, 5, 10, 15, 20),
        .model = 'tree'
    )

    c50.fit <- train(
        X,
        y,
        method = 'C5.0',
        metric = metric,
        trControl = train.control,
        tuneGrid = c50.grid
    )

    c50.params <- as.list(c50.fit$bestTune)
    c50.params$x <- X
    c50.params$y <- y.target

    c50.bestfit <- do.call(C5.0, c50.params)

    gbm.fit <- train(
        X,
        y,
        method = 'gbm',
        metric = metric,
        trControl = train.control,
        tuneGrid = gbmGrid,
        verbose = FALSE
        )

    list(
        rpart.fit = rpart.fit,
        c50.fit = c50.fit,
        c50.bestfit = c50.bestfit,
        gbm.fit = gbm.fit
        )
    }
