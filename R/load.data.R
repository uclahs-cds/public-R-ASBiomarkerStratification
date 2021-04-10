factor.ISUP <- function(x) {
  x <- as.factor(x);
  levels(x) <- c('No cancer', 'Low', 'Intermediate favorable', 'Intermediate unfavorable', 'High', 'Very High');
  x;
  }

factor.Gleason <- function(x) {
  factor(x, levels = c('0+0', '3+3', '3+4', '4+3', '4+4', '5+5'), ordered = TRUE);
  }

factor.PIRADS <- function(x) {
  factor(x, levels = 0:5, ordered = TRUE);
  }

mutation.dummy <- function(x) {
  mutation1 <- as.integer(x['Mutation1']);
  mutation2 <- as.integer(x['Mutation.2']);
  if (is.na(mutation1) && is.na(mutation2)) {
    # Missing value, make all dummy variables missing as well
    res <- list(Mutation_BRCA1 = NA, Mutation_BRCA2 = NA, Mutation_ATM = NA, Mutation_MLH1 = NA, Mutation_PMS2 = NA, Germline.variants = NA);
    }
  else {
    res <- list(BRCA1 = 0, BRCA2 = 0, ATM = 0, MLH1 = 0, PMS2 = 0);
    if (mutation1 > 0) res[mutation1] <- 1
    if (mutation2 > 0) res[mutation2] <- 1
    germline.variants <- paste0(names(res)[res > 0], collapse = ' & ')
    if (germline.variants == '') germline.variants <- 'None'
    names(res) <- paste0('Mutation_', names(res))
    res$Germline.variants <- germline.variants
    }
  res
  }

#' Loads the data for the Active Surveillance project
#'
#' @param biomark.pathBiomarkers xlsx file
#' @param biomark.key.path path to the Biomarkers key xlsx file
#' @param genetics path to the genetics xlsx file
#' @param onlyBiodb boolean indicating if we should only return the biodb (main data)
#' @export
load.data.AS <- function(biomark.path,
                         biomark.key.path,
                         genetics.path,
                         biomark.categories.path,
                         onlyBiodb = FALSE) {
  factor.cols <- c('Race', 'MRIResult',
                   'Ethnicity',
                   'BiopsyResult','Observation',
                   'GeneticAncestry', 'GeneticRiskCategory',
                   # 'Mutation1', 'Mutation.2', 'GlobalScreeningArray', 'BRCAMutation',
                   'RSIlesionCancer',  'RSIlesionUpgraded',
                   'RSIlesionObservation',
                   'ProgressedToTreatment', 'BiopsyUpgraded',
                   'Prostatectomy',
                   'UpgradedAndProgressed', 'NoUpgradeAndProgressed'
                   );

  ISUP.cols <- c('PreviousISUP',  'RSIlesionISUP', 'StudyHighestISUP');
  Gleason.cols <- c('PreviousGleason', 'StudyHighestGleason', 'RSIlesionGleason');
  PIRADS.cols <- c('RSIlesionPIRADS', 'HighestPIRADS');

  mutation.cols <- c('Mutation1', 'Mutation.2');

  targets <- c('ProgressedToTreatment', 'BiopsyUpgraded', 'Prostatectomy');

  # Open Raw data:
  biodb <- xlsx::read.xlsx(
    biomark.path,
    sheetIndex = 1,
    header = TRUE,
    stringsAsFactors = FALSE
  );

  # Create dummy variables from the 2 mutation columns
  mutation.dummy.vars <- do.call(rbind.data.frame,
                                 apply(biodb[, mutation.cols], 1, mutation.dummy)
                                 );


  biodb <- cbind.data.frame(biodb, mutation.dummy.vars);

  # Remove the old mutation columns
  biodb[, mutation.cols] <- NULL

  biodb$Germline.variants <- relevel(
    x = as.factor(biodb$Germline.variants),
    ref = 'None')

  #mutation.combinations <- unique(biodb[, c('Mutation_BRCA1', 'Mutation_BRCA2', 'Mutation_ATM', 'Mutation_MLH1', 'Mutation_PMS2')])
  #rownames(mutation.combinations) <- NULL

  # Convert to factor
  biodb[,factor.cols] <- lapply(biodb[,factor.cols], as.factor);

  # Convert ISUP and Gleason to factors with pre-specified levels
  biodb[,ISUP.cols] <-  lapply(biodb[, ISUP.cols], factor.ISUP);
  biodb[,Gleason.cols] <-  lapply(biodb[, Gleason.cols], factor.Gleason);
  biodb[,PIRADS.cols] <- lapply(biodb[, PIRADS.cols], factor.PIRADS);

  # Add PSA and PHI Density
  biodb$PSADensity <- biodb$SOCPSA / biodb$ProstateVolume;
  biodb$PHIDensity <- biodb$PHI / biodb$ProstateVolume;

  # Add MyProstateScore = MiPSHighGradeCancerRisk
  biodb$MyProstateScore <- biodb$MiPSHighGradeCancerRisk;
  biodb$MPSDensity <- biodb$MyProstateScore / biodb$ProstateVolume;

  # Update levels
  levels(biodb$Race) <- c('White', 'African-American', 'Asian');
  levels(biodb$Ethnicity) <- c('Non-Hispanic', 'Hispanic');
  levels(biodb$GeneticAncestry) <- c('European', 'African', 'East Asian', 'Native American');
  levels(biodb$GeneticRiskCategory) <- c('Low', 'Normal', 'High');
  levels(biodb$MRIResult) <- c('No Legion', 'Legion Found');
  levels(biodb$HighestPIRADS) <- c('No lesion', 'Very low', 'Low', 'Intermediate', 'High', 'Very high');
  levels(biodb$BiopsyResult) <- c('Negative', 'Positive');
  levels(biodb$Observation) <- c('MRI Positive/Biopsy Positive', 'MRI Positive/Biopsy Negative',
                                 'MRI Negative/Biopsy Positive', 'MRI Negative/Biopsy Negative');

  # Note: When we remove the NoUpgradeAndProgressed then BiopsyUpgrade represents any aggressive disease progression.
  # New target that combines all of the other targets
  # biodb$AggressiveDisease <- as.factor(apply(biodb[, targets], 1, function(r) any(r == 1, na.rm = TRUE)))
  # targets <- c(targets, 'AggressiveDisease')

  # Re-levels the targets to no/yes levels
  biodb[, targets] <- lapply(biodb[, targets], `levels<-`, value = c('no', 'yes'))

  # Rename Ethnicity to Hispanic
  biodb$Hispanic <- as.factor(biodb$Ethnicity);

  attr(biodb$Weight, 'label') <- 'Weight (kg)';
  attr(biodb$Height, 'label') <- 'Height (cm)';
  attr(biodb$ProstateVolume, 'label') <- 'Prostate Volume (cm^3)';
  attr(biodb$p2PSA, 'label') <- '[-2]proPSA';
  attr(biodb$freePSA, 'label') <- 'Free PSA';
  attr(biodb$MRIResult, 'label') <- 'MRI Result';
  attr(biodb$MRILesions, 'label') <- 'MRI Lesions';
  attr(biodb$HighestPIRADS, 'label') <- 'PI-RADS';
  attr(biodb$BiopsyResult, 'label') <- 'Biopsy Result';
  attr(biodb$ADCnormalSignal, 'label') <- 'ADC normal Signal';
  attr(biodb$ADClesionSignal, 'label') <- 'ADC lesion Signal';
  attr(biodb$RSIlesionSignal, 'label') <- 'RSI lesion Signal';
  attr(biodb$RSInormalSignal, 'label') <- 'RSI normal Signal';
  attr(biodb$RSIlesionPIRADS, 'label') <- 'PI-RADS'
  attr(biodb$MiPSCancerRisk, 'label') <- 'MiPS Cancer Risk';
  attr(biodb$MiPSHighGradeCancerRisk, 'label') <- 'MyProstateScore';
  attr(biodb$MyProstateScore, 'label') <- 'MyProstateScore';
  attr(biodb$MPSDensity, 'label') <- 'MPS Density';
  attr(biodb$SOCPSA, 'label') <- 'PSA';
  attr(biodb$PSAHyb, 'label') <- 'PSA Hybrid';
  attr(biodb$PSADensity, 'label') <- 'PSA Density';
  attr(biodb$PercentFreePSA, 'label') <- '% Free PSA';
  attr(biodb$PHIDensity, 'label') <- 'PHI Density';
  attr(biodb$TNFaAverage, 'label') <- 'TNFa';
  attr(biodb$GeneticRiskScore, 'label') <- 'Genetic Risk Score';
  attr(biodb$GeneticRiskCategory, 'label') <- 'Genetic Risk Category';
  attr(biodb$GlobalScreeningArray, 'label') <- 'Global Screening Array';
  attr(biodb$GSAPositives, 'label') <- 'GSA Positives';
  attr(biodb$BRCAMutation, 'label') <- 'BRCA Mutation';
  attr(biodb$Mutation_BRCA1, 'label') <- 'BRCA1';
  attr(biodb$Mutation_BRCA2, 'label') <- 'BRCA2';
  attr(biodb$Mutation_ATM, 'label') <- 'ATM';
  attr(biodb$Mutation_MLH1, 'label') <- 'MLH1';
  attr(biodb$Mutation_PMS2, 'label') <- 'PMS2';
  attr(biodb$Germline.variants, 'label') <- 'Deleterious germline variants';
  attr(biodb$PreviousISUP, 'label') <- 'Previous ISUP';
  attr(biodb$StudyHighestGleason, 'label') <- 'Study Highest Gleason';
  attr(biodb$StudyHighestISUP, 'label') <- 'Study Highest ISUP';
  attr(biodb$ProgressedToTreatment, 'label') <- 'Progressed to Treatment';

  biodb$PHI.computed <- with(biodb, {
    (p2PSA / freePSA) * sqrt(PSAHyb)
  });

  biodb$FollowUpTime <- biodb$DaysDxToLastReview;
  biodb$FollowUpTime[!is.na(biodb$DaysDxToUpgrade)] <- biodb$DaysDxToUpgrade[!is.na(biodb$DaysDxToUpgrade)];

  biokey <- xlsx::read.xlsx(
    biomark.key.path,
    sheetIndex = 1,
    header = TRUE,
    stringsAsFactors = FALSE
    );

  biokey <- biokey[!is.na(biokey$Column.ID), ];

  biokey.gen <- xlsx::read.xlsx(
    genetics.path,
    sheetIndex = 1,
    header = TRUE,
    stringsAsFactors = FALSE
    );

  bio.categories <- xlsx::read.xlsx(
    biomark.categories.path,
    sheetIndex = 1,
    header = TRUE,
    stringsAsFactors = FALSE
    );

  bio.categories$Category <- as.factor(bio.categories$Category);

  if (onlyBiodb) {
    biodb
    }
  else {
    list(
      biodb = biodb,
      biokey = biokey,
      biokey.gen = biokey.gen,
      bio.categories = bio.categories
      );
    }
  }

#' Loads the data for analysis with default file names
#' @param onlyBiodb boolean indicating if we should only return the biodb (main data)
#' @export
default.load.data <- function(onlyBiodb = FALSE) {
  file.names <- c(
    'MRI DOD Biomarkers Database_Boutros - 2020.11.2.xlsx',
    'MRI DOD Biomarkers Database_Boutros - Key.xlsx',
    'MRI DOD Biomarkers Genetics_Boutros - 2019.11.12.xlsx',
    'MRI DOD Biomarkers Category_Boutros.xlsx'
    );
  file.paths <- here::here(paste0('data/', file.names));
  do.call('load.data.AS', c(as.list(file.paths), onlyBiodb));
  }

#' Read in the biomarker categories
#'
#' @param file.name
#'
#' @return
#' @export
#'
#' @examples
load.biomarker.categories <- function(file.name = 'biomarkers_categories.xlsx') {
  res <- xlsx::read.xlsx(
    here::here(paste0('data/', file.name)),
    sheetIndex = 1,
    header = TRUE,
    stringsAsFactors = FALSE
    );
  rownames(res) <- res$variable;
  res[order(res$order), ];
  }
