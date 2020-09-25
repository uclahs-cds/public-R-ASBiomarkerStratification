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
    metric = c('F', 'Accuracy', 'AUC', 'Kappa', 'Precision', 'Recall', 'ROC',
               'Sens', 'Spec'),
    exclude.vars = NULL,
    predict.missing = TRUE,
    seed = NULL) {
    target <- match.arg(target);
    metric <- match.arg(metric);

    if(!is.null(seed)) set.seed(seed);

    biokey.variables <- c('Age',
                          'Race',
                          'Hispanic',
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
                          'TNFaAverage', 'TNFaSTD', # Tumor necrosis factor
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

        # Change name to multiclass name
        if(predict.missing && 'F' == metric) {
            metric <- 'Mean_F1'
        }
    }

    biokey.variables <- setdiff(biokey.variables, exclude.vars)

    missing.target <- is.na(biodb[, target])

    if(predict.missing) {
        X <- biodb[, biokey.variables];
        y.target <- as.character(biodb[, target]);
        y.target[missing.target] <- "Missing";
        y.target <- as.factor(y.target);

        y <- y.target;
        levels(y) <- c('no', 'yes', 'Missing');
    } else {
        # Variable we want to predict
        y.target <- biodb[!missing.target, target];
        y <- y.target
        levels(y) <- c('no', 'yes');

        # Drop the missing target values
        X <- biodb[!missing.target, biokey.variables]
    }

    gbm.grid <-  expand.grid(interaction.depth = c(1, 5, 9),
                            n.trees = (1:30)*50,
                            shrinkage = c(0.001, 0.01, 0.1),
                            n.minobsinnode = 10 # 20
                            )

    xgb.grid <- expand.grid(
        nrounds = seq(from = 200, to = 1000, by = 50),
        eta = c(0.025, 0.05, 0.1, 0.3),
        max_depth = c(2, 3, 4, 5, 6),
        gamma = 0,
        colsample_bytree = 1,
        min_child_weight = 1,
        subsample = 1
        );

    custom.summary <- function (data, lev = NULL, model = NULL) {
        if(length(levels(y)) == 2) {
            c(
                defaultSummary(data, lev, model),
                twoClassSummary(data, lev, model),
                prSummary(data, lev, model)
                )
        } else if(length(levels(y)) > 2) {
            c(
                defaultSummary(data, lev, model),
                multiClassSummary(data, lev, model)
            )
        }

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

    # XGB needs numeric input.
    # Convert ordered factors to numeric
    # ordered.cols <- unlist(lapply(biodb, is.ordered))
    X.ints <- X
    # biodb.ints[, ordered.cols] <- lapply(biodb.ints[, ordered.cols], function(x) as.numeric(x))
    # Convert ordered variable to numeric
    X.ints$PreviousISUP <- as.numeric(as.character(X$PreviousISUP))

    # Convert to dummy variables
    dummy.formula <- paste0("~ ", paste0(biokey.variables, collapse = " + "))
    X.dummy.ints <- predict(dummyVars(dummy.formula, data = X.ints, fullRank = TRUE), newdata = X)

    print("Fitting XGB model...");
    xgb.fit <- train(
        X.dummy.ints,
        y,
        trControl = train.control,
        tuneGrid = xgb.grid,
        metric = metric,
        method = "xgbTree"
    )

    print("Completed fitting XGB model...");

    model.file <- paste('xgb', target, metric, seed, 'model.RDS', sep = "_");
    saveRDS(object = xgb.fit, file = here::here(paste0('models/', model.file)));

    # Add up-sampling for rpart
    train.control.up <- train.control;
    train.control.up$sampling <- "up";


    print("Fitting rpart Model...");
    rpart.fit <- train(
        X,
        y,
        method = 'rpart',
        metric = metric,
        trControl = train.control.up,
        tuneLength = 30
    )
    print("Completed fitting rpart model...");

    model.file <- paste('rpart', target, metric, seed, 'model.RDS', sep = "_");
    saveRDS(object = rpart.fit, file = here::here(paste0('models/', model.file)));

    print("Fitting gbm Model...");
    gbm.fit <- train(
        X,
        y,
        method = 'gbm',
        metric = metric,
        trControl = train.control,
        tuneGrid = gbm.grid,
        verbose = FALSE
        )
    print("Completed fitting gbm model...");

    model.file <- paste('gbm', target, metric, seed, 'model.RDS', sep = "_");
    saveRDS(object = gbm.fit, file = here::here(paste0('models/', model.file)));

    list(
        rpart.fit = rpart.fit,
        xgb.fit = xgb.fit,
        gbm.fit = gbm.fit
        )
    }
