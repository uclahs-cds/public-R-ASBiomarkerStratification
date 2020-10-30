#' Creates a correlation heatmap for the AS cohort
#'
#' @param biodb
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
create.heatmap.AS <- function(biodb, ...) {
    # variables to be used for heatmap:
    variables <- c(
        'Age',
        'Weight',
        'Height',
        'BMI',
        'ProstateVolume',
        'MRIResult',
        'MRILesions',
        'HighestPIRADS',
        'BiopsyResult',
        'ADCnormalSignal',
        'ADClesionSignal',
        'RSIlesionSignal',
        'RSInormalSignal',
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
        'GeneticRiskCategory',
        'GlobalScreeningArray',
        'GSAPositives',
        'BRCAMutation',
        "Mutation_BRCA1", "Mutation_BRCA2", "Mutation_ATM", "Mutation_MLH1", "Mutation_PMS2"
    );

    labels <- label.or.name(biodb[, variables])

    # Computing the Correlations for heatmap:
    heatmap.data <- vector(
        mode = "list",
        length = length(variables)
    );

    for (i in 1:length(variables) ){
        for (j in 1:length(variables)  ){
            biodbA <- as.numeric(biodb[, variables[i] ]);
            biodbB <- as.numeric(biodb[, variables[j] ]);
            heatmap.data[[i]] <- c(
                heatmap.data[[i]],
                cor(
                    biodbA,
                    biodbB,
                    method = "spearman",
                    use = "complete.obs"
                )
            );
        }
    }

    # Data for Heatmap:
    simple.data <- as.data.frame(
        do.call(cbind, heatmap.data)
    );
    colnames(simple.data) <- variables;
    rownames(simple.data) <- variables;

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
    BoutrosLab.plotting.general::create.heatmap(
        x = simple.data,
        xaxis.lab = labels,
        xaxis.cex = 1,
        yaxis.lab = labels,
        yaxis.cex = 1,
        colourkey.cex = 1,
        covariates = sample.covariate,
        covariates.top = top.covariate,
        covariate.legend = sample.cov.legend,
        legend.side = 'right',
        legend.title.cex = 2,
        legend.cex = 1.75,
        legend.title.just = 'left',
        legend.between.row = 0.1,
        legend.border.padding = 3,
        colourkey.labels.at = seq(-1, 1, 0.2),
        at = key.scale,
        axis.xlab.padding = 13,
        #  clustering.method = 'none',
        plot.dendrograms = 'right',
        right.dendrogram.size = 5,
        height = 18,
        width = 18,
        left.padding = 5,
        bottom.padding = 3,
        resolution = 1000,
        ...
    );
}
