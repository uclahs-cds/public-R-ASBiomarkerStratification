#' Models for AS endpoints
#'
#' @param biodb
#' @param train.control
#' @param target
#' @param metric
#' @param models
#' @param exclude.vars
#' @param include.vars Include additional variables for the models
#' @param predict.missing
#' @param seed
#' @param rpart.cost
#' @param rm.NoUpgradeAndProgressed
#' @param reduced.model Use the reduced set of features for the model
#' @param baseline.model Use the clinically relevant demographic features for the model
#'
#' @return
#' @export
#'
#' @examples
AS.models <- function(
    biodb,
    train.control = NULL,
    target = c('ProgressedToTreatment', 'BiopsyUpgraded', 'Prostatectomy'),
    metric = c('F', 'Accuracy', 'PR-AUC', 'Kappa', 'Precision', 'Recall', 'ROC-AUC',
               'Sens', 'Spec'),
    models = c('xgb', 'gbm', 'rpart'),
    exclude.vars = NULL,
    include.vars = NULL,
    predict.missing = FALSE,
    seed = NULL,
    rpart.cost = NULL,
    rm.NoUpgradeAndProgressed = FALSE,
    reduced.model = FALSE,
    baseline.model = FALSE,
    suffix = ''
    ) {
    target <- match.arg(target);
    metric <- match.arg(metric);
    models <- match.arg(models, several.ok = TRUE);
    biomarkers <- load.biomarker.categories();

    if (baseline.model) {
        biokey.variables <- biomarkers$variable[biomarkers$clinically.useful == 1 & biomarkers$category == 'Demographics'];
        }
    else if (reduced.model) {
        biokey.variables <- biomarkers$variable[biomarkers$clinically.useful == 1];
        }
    else {
        biokey.variables <- biomarkers$variable;
        }

    if ('BiopsyUpgraded' == target) {
        exclude.vars <- c(exclude.vars, 'BiopsyResult');

        # Change name to multiclass name
        if (predict.missing && 'F' == metric) {
            metric <- 'Mean_F1';
            }
        }

    # Part of file name
    model.id <- paste(target, metric, seed, sep = '_');
    if (!is.null(suffix)) {
        model.id <- paste(model.id, suffix, sep = '_');
        }

    biokey.variables <- setdiff(biokey.variables, exclude.vars);

    if (!is.null(include.vars)) {
        biokey.variables <- union(biokey.variables, include.vars)
        }

    missing.target <- is.na(biodb[, target]);

    # Only keep the patient that did not leave AS voluntarily if we have the rm.NoUpgradeAndProgressed flag
    valid.patients <- ! rm.NoUpgradeAndProgressed |
        (! (biodb$NoUpgradeAndProgressed == 1) |
        is.na(biodb$NoUpgradeAndProgressed));

    if (predict.missing) {
        X <- biodb[, biokey.variables];
        y.target <- as.character(biodb[, target]);
        y.target[missing.target] <- 'Missing';
        y.target <- as.factor(y.target);

        y <- y.target;
        levels(y) <- c('no', 'yes', 'Missing');
        }
    else {
        # Variable we want to predict
        y.target <- biodb[!missing.target & valid.patients, target];
        y <- y.target
        levels(y) <- c('no', 'yes');
        # Most functions expect the positive case to be the first factor
        # factor(y, levels = c('yes', 'no'));

        # Drop the missing target values
        X <- biodb[!missing.target & valid.patients, biokey.variables]
        }

    gbm.grid <- gbm.hyper.grid();

    xgb.grid <- expand.grid(
        nrounds = seq(from = 200, to = 1000, by = 50),
        eta = c(0.025, 0.05, 0.1, 0.3),
        max_depth = c(4, 5, 6),
        gamma = 0,
        colsample_bytree = 1,
        min_child_weight = 1,
        subsample = 1
        );

    if (is.null(train.control)) {
        train.control <- trainControl(
            method = 'repeatedcv',
            number = 10,
            repeats = 5,
            classProbs = TRUE,
            savePredictions = T,
            summaryFunction = custom.summary
            # https://topepo.github.io/caret/subsampling-for-class-imbalances.html
            # Don't seem to be enough of a class imbalance to make a difference
            # sampling = 'up'
            )
        }
    else {
        train.control$savePredictions <- TRUE;
        train.control$classProbs <- TRUE;
        train.control$summaryFunction <- custom.summary;
        }

    # XGB needs numeric input.
    # Convert ordered factors to numeric
    X.ints <- X;
    # Convert ordered variable to numeric
    if ('PreviousISUP' %in% colnames(X.ints)) {
        X.ints$PreviousISUP <- as.numeric(as.character(X$PreviousISUP));
        }

    # Convert to dummy variables
    dummy.formula <- paste0('~ ', paste0(biokey.variables, collapse = ' + '));
    X.dummy.ints <- predict(caret::dummyVars(dummy.formula, data = X.ints, fullRank = TRUE), newdata = X.ints);

    print(paste('Fitting models for:', target, 'optimizing', metric));

    if ('xgb' %in% models) {
        print('Fitting XGB model...');
        if (!is.null(seed)) set.seed(seed);
        xgb.fit <- caret::train(
            X.dummy.ints,
            y,
            trControl = train.control,
            tuneGrid = xgb.grid,
            metric = metric,
            method = 'xgbTree'
        );

        print('Completed fitting XGB model...');

        model.file <- paste('xgb', model.id, 'model.RDS', sep = '_');
        print(paste0('Saving file to: ',  here::here(paste0('models/', model.file))));
        saveRDS(object = xgb.fit, file = here::here(paste0('models/', model.file)));
    }

    if ('rpart' %in% models) {
        if (!is.null(rpart.cost)) {
            print('Fitting rpart model with cost matrix...');
            print(rpart.cost)
            }
        else {
            print('Fitting rpart model ...');
            }


        if (!is.null(seed)) set.seed(seed);
        rpart.fit <- caret::train(
            X,
            y,
            method = 'rpart',
            metric = metric,
            trControl = train.control,
            tuneLength = 30,
            parms = list(loss = rpart.cost)
            );
        print('Completed fitting rpart model...');

        model.file <- paste('rpart', model.id, 'model.RDS', sep = '_');
        saveRDS(object = rpart.fit, file = here::here(paste0('models/', model.file)));
    }

    if ('gbm' %in% models) {
        if (!is.null(seed)) set.seed(seed);

        print('Fitting gbm model');
        gbm.fit <- caret::train(
            X,
            y,
            method = 'gbm',
            metric = metric,
            trControl = train.control,
            tuneGrid = gbm.grid,
            verbose = FALSE
            );
        print('Completed fitting gbm model...');

        model.file <- paste('gbm', model.id, 'model.RDS', sep = '_');
        # print(paste0('Saving file to: ',  here::here(paste0('models/', model.file))));
        saveRDS(object = gbm.fit, file = here::here(paste0('models/', model.file)));
    }
}


#' Custom summary for caret model
#'
#' @param data
#' @param lev
#' @param model
#'
#' @return
#' @export
#'
#' @examples
custom.summary <- function(data, lev = NULL, model = NULL) {
    if (length(lev) == 2) {
        two.class.sum <- caret::twoClassSummary(data, lev, model);
        names(two.class.sum) <- c('ROC-AUC', 'Sens', 'Spec');
        pr.summary <- caret::prSummary(data, lev, model);
        names(pr.summary) <- c('PR-AUC', 'Precision', 'Recall', 'F');
        c(
            caret::defaultSummary(data, lev, model),
            two.class.sum,
            pr.summary
        )
        }
    else if (length(lev) > 2) {
        c(
            caret::defaultSummary(data, lev, model),
            caret::multiClassSummary(data, lev, model)
        )
        }
    }


#' GBM hyper parameter grid
#'
#' @return
#' @export
#'
#' @examples
gbm.hyper.grid <- function() {
    expand.grid(interaction.depth = c(1, 3, 5),
                n.trees = seq(200, 1300, by = 50),
                shrinkage = c(0.001, 0.01, 0.1),
                n.minobsinnode = 10)
    }

#' Caret train control settings
#' @export
AS.train.control <- caret::trainControl(
    method = 'repeatedcv',
    number = 10,
    repeats = 5,
    classProbs = TRUE,
    savePredictions = T,
    summaryFunction = custom.summary
    );
