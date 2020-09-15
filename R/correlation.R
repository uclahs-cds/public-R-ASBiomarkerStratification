#' Title
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
        'Mutation1',
        'Mutation.2'
    );


    real.names.set <- c(
        'Age',
        'Weight',
        'Height',
        'BMI',
        'Prostate Volume',
        'MRI Result',
        'MRI Lesions',
        'Highest PIRADS',
        'Biopsy Result',
        'ADC normal Signal',
        'ADC lesion Signal',
        'RSI lesion Signal',
        'RSI normal Signal',
        'PCA3',
        'T2ERG',
        'MiPS Cancer Risk',
        'MiPS High Grade Cancer Risk',
        'PSA Hybrid',
        'free PSA',
        'PSA Density',
        'p2PSA',
        'Percent Free PSA',
        'PHI',
        'PHI Density',
        'SOCPSA',
        'TNFa Average',
        'Genetic Risk Score',
        'Genetic Risk Category',
        'Global Screening Array',
        'GSA Positives',
        'BRCA Mutation',
        'Mutation 1',
        'Mutation 2'
    );


    # Adding PSA Density and PHI Density:
    biodb.added <- data.frame(
        biodb,
        PSADensity = biodb$freePSA/biodb$ProstateVolume,
        PHIDensity = biodb$PHI/biodb$ProstateVolume
    );

    # Computing the Correlations for heatmap:
    heatmap.data <- vector(
        mode = "list",
        length = length(variables)
    );

    for (i in 1:length(variables) ){
        for (j in 1:length(variables)  ){
            biodbA <- as.numeric(biodb.added[, variables[i] ]);
            biodbB <- as.numeric(biodb.added[, variables[j] ]);
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
            colours = c('black'),
            labels = c(''),
            title = 'Patients Features'
        ),
        legend = list(
            colours = c(
                'dodgerblue',
                'gold',
                'firebrick3',
                'darkgreen'
            ),
            labels = c(
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
        xaxis.lab = real.names.set,
        xaxis.cex = 0.75,
        yaxis.lab = real.names.set,
        yaxis.cex = 0.75,
        colourkey.cex = 0.75,
        covariates = sample.covariate,
        covariates.top = top.covariate,
        covariate.legend = sample.cov.legend,
        legend.side = 'right',
        legend.title.cex = 1,
        legend.cex = 1,
        legend.title.just = 'left',
        legend.between.row = 0.1,
        legend.border.padding = 0.5,
        colourkey.labels.at = seq(-1, 1, 0.2),
        at = key.scale,
        axis.xlab.padding = 4,
        #  clustering.method = 'none',
        plot.dendrograms = 'right',
        height = 9,
        width = 9,
        left.padding = 4,
        resolution = 300,
        ...
    );
}
