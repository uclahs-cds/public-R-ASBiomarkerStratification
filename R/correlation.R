# variables to be used for heatmap
cor.variables <- c(
    'Age',
    'Weight',
    'Height',
    'BMI',
    'ProstateVolume',
    'HighestPIRADS',
    'ADClesionSignal',
    'RSIlesionSignal',
    'PCA3',
    'T2ERG',
    'MiPSCancerRisk',
    'MiPSHighGradeCancerRisk',
    'PSAHyb',
    'freePSA',
    'PSADensity',
    'p2PSA',
    'PercentFreePSA',
    'PHI',
    'PHIDensity',
    'SOCPSA',
    'TNFaAverage',
    'GeneticRiskScore',
    'GlobalScreeningArray',
    'GSAPositives',
    'BRCAMutation',
    "Mutation_BRCA1", "Mutation_BRCA2", "Mutation_ATM", "Mutation_MLH1", "Mutation_PMS2"
);

create.forestplot <- function(biodb, ...) {
    # Compute the univariate effect-sizes (AUROC)
    uni.auc.ci <- lapply(biodb[, cor.variables], function(predictor) {
        pROC::ci.auc(
            roc(response = biodb$BiopsyUpgraded,
                predictor= predictor,
                plot = FALSE,
                direction = '<',
                levels = c('no', 'yes')))
    })

    labels <- label.or.name(biodb[, cor.variables])

    segplot.data <- as.data.frame(do.call(rbind, uni.auc.ci))
    colnames(segplot.data) <- c('min', 'point', 'max')
    segplot.data$variable <- as.factor(labels)
    segplot.data$significant <- segplot.data$min > 0.5

    sig.colours <- c("grey", "darkorange1")

    legend <- list(
        inside = list(
            fun = draw.key,
            args = list(
                key = list(
                    points = list(
                        col = rev(sig.colours),
                        pch = 19,
                        lty = 1,
                        cex = 1.25
                    ),
                    text = list(
                        lab = c("p < 0.05", "p >= 0.05")
                    ),
                    padding.text = 1,
                    cex = 1.25
                )
            ),
            x = 0.8,
            y = 0.15,
            corner = c(0,1)
        )
    )

    create.segplot(
        formula = reorder(variable, point) ~ min + max,
        data = segplot.data,
        xlimits = c(0, 1),
        level = segplot.data$significant,
        col.regions = c("grey", "darkorange1"),
        centers = segplot.data$point,
        xat = seq(0, 1, by = 0.2),
        abline.col = "lightgrey",
        abline.lty = 2,
        abline.v = 0.5,
        ylab.label = '',
        xlab.label = 'AUROC',
        legend = legend,
        height = 8,
        width = 12,
        left.padding = 10,
        ...
        # This doesn't work right. How are the y-axis values ordered?
        # yaxis.col = ifelse(segplot.data$significant, 'darkorange1', 'black')
    )
}

#' Creates a correlation heatmap for the AS cohort
#'
#' @param biodb the data frame of the patients
#' @param ... Additional arguments to pass to `BoutrosLab.plotting.general::create.heatmap`
#'
#' @return Output from `BoutrosLab.plotting.general::create.heatmap`
#' @export
#'
#' @examples
#' biodb <- default.load.data(onlyBiodb = TRUE)
#' create.heatmap.AS(biodb)
#' create.heatmap.AS(biodb,
#'   filename = here('figures/corr_heatmap.tiff'))
create.heatmap.AS <- function(biodb, ...) {
    labels <- label.or.name(biodb[, cor.variables])

    # Computing the Correlations for heatmap:
    heatmap.data <- vector(
        mode = "list",
        length = length(cor.variables)
    );

    numeric.biodb.heatmap.vars <- lapply(biodb[, cor.variables], as.numeric);
    target.corr.data <- unlist(lapply(numeric.biodb.heatmap.vars, cor, y = as.numeric(biodb$BiopsyUpgraded), method = "spearman", use = "complete.obs"))

    simple.data <- cor(as.data.frame(numeric.biodb.heatmap.vars), method = "spearman", use = "complete.obs")

    # Scaling the colors:
    key.min = -1;
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

    sample.cov.legend <- list(
        legend = list(
            colours = c(
                'black',
                'dodgerblue',
                'gold',
                'firebrick3',
                'darkgreen'
            ),
            labels = c(
                'Patient Features',
                'Imaging',
                'Urine',
                'Blood/Urine',
                'Genetics'
            ),
            title = 'Test Methodology'
        )
    );



    chr.cov.colours <- c(
        rep('black',7),
        rep('dodgerblue',8),
        rep('gold',2),
        rep('firebrick3', 11),
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

    # Plotting Heatmap:
    corr.heatmap <- BoutrosLab.plotting.general::create.heatmap(
        x = simple.data,
        xaxis.lab = labels,
        xaxis.cex = 2,
        yaxis.lab = labels,
        yaxis.cex = 2,
        colourkey.cex = 2,
        # covariates = sample.covariate,
        #covariates.top = top.covariate,
        #covariate.legend = sample.cov.legend,
        legend.side = 'right',
        legend.title.cex = 2.2,
        legend.cex = 2,
        legend.title.just = 'left',
        legend.between.row = 0.1,
        legend.border.padding = 3,
        colourkey.labels.at = seq(-1, 1, 0.2),
        at = key.scale,
        axis.xlab.padding = 13,
        plot.dendrograms = 'none',
        height = 18,
        width = 23,
        left.padding = 20,
        bottom.padding = 3,
        resolution = 200
    );

    bx.upgrade.corr <- create.heatmap(
        x = t(target.corr.data[corr.heatmap$x.limits]),
        clustering.method = 'none',
        at = key.scale,
        print.colour.key = FALSE,
        grid.row = TRUE,
        grid.col = TRUE,
        yaxis.tck = 0,
        yaxis.lab = 'Biopsy Upgraded'
    )

    bx.upgrade.barplot.data <- data.frame(
        x = corr.heatmap$x.limits,
        y = target.corr.data[corr.heatmap$x.limits]
    )

    max_cor <- max(abs(target.corr.data))
    bx.ylim <- round(max_cor, digits = 1)

    bx.upgrade.corr.barplot <- create.barplot(
        formula = y ~ x,
        data = bx.upgrade.barplot.data,
#        yat = seq(-bx.ylim, bxy, by = 0.4),
        ylimits = c(-bx.ylim, bx.ylim)
    )

    create.multiplot(
        plot.objects = list(corr.heatmap, bx.upgrade.corr.barplot),
        # panel.heights = c(0.05, 1),
        panel.heights = c(0.5, 1),
        y.relation = 'free',
        height = 18,
        width = 23,
        left.padding = 30,
        right.padding = 2,
        bottom.padding = 5,
        xlab.to.xaxis.padding = 18,
        legend = list(
            right = list(
                x = 0.10,
                y = 0.50,
                fun = legend.grob(sample.cov.legend)
            )
        ),
        ylab.label = list(
            'Biopsy Upgraded\nCorrelations',
            '',
            ''
        ),
        # yat = c(
        #    seq_along(corr.heatmap$y.scales$labels)),
        xaxis.rot = 90,
        # yaxis.lab = list(
        #      corr.heatmap$y.scales$labels,
        #      'Biopsy Upgraded'
        #  ),
        ...
    )
}
