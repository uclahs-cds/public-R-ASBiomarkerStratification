factor.ISUP <- function(x) {
  factor(x, levels = 0:5, ordered = TRUE)
  }

factor.Gleason <- function(x) {
  factor(x, levels = c("0+0", "3+3", "3+4", "4+3", "4+4", "5+5"), ordered = TRUE)
}

factor.PIRADS <- function(x) {
  factor(x, levels = 0:5, ordered = TRUE)
}

mutation.dummy <- function(x) {
  mutation1 <- as.integer(x["Mutation1"]);
  mutation2 <- as.integer(x["Mutation.2"]);
  if(is.na(mutation1) && is.na(mutation2)) {
    # Missing value, make all dummy variables missing as well
    res <- list(Mutation_BRCA1 = NA, Mutation_BRCA2 = NA, Mutation_ATM = NA, Mutation_MLH1 = NA, Mutation_PMS2 = NA);
  }
  else {
    res <- list(Mutation_BRCA1 = 0, Mutation_BRCA2 = 0, Mutation_ATM = 0, Mutation_MLH1 = 0, Mutation_PMS2 = 0);
    # What happens when we do res[0] <- 1? Doesn't seem to do anything
    if(mutation1 > 0) res[mutation1] <- 1
    if(mutation2 > 0) res[mutation2] <- 1
    }
  res
  }

# Project: Prostate Cancer Active Surveillance

#######################################################
# Put input data into objects:
#######################################################

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

  mutation.cols <- c("Mutation1", "Mutation.2");

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

  # Convert to factor
  biodb[,factor.cols] <- lapply(biodb[,factor.cols], as.factor);

  # Convert ISUP and Gleason to factors with pre-specified levels
  biodb[,ISUP.cols] <-  lapply(biodb[, ISUP.cols], factor.ISUP);
  biodb[,Gleason.cols] <-  lapply(biodb[, Gleason.cols], factor.Gleason);
  biodb[,PIRADS.cols] <- lapply(biodb[, PIRADS.cols], factor.PIRADS);

  # Add PSA and PHI Density
  biodb$PSADensity <- biodb$PSAHyb / biodb$ProstateVolume;
  biodb$PHIDensity = biodb$PHI / biodb$ProstateVolume;

  # TODO: Add all other binary variables as logical
  # biodb$ProgressedToTreatment <- as.logical(biodb$ProgressedToTreatment)
  # biodb$BiopsyUpgraded <- as.logical(biodb$BiopsyUpgraded)

  # Update levels
  levels(biodb$Race) <- c('White', 'African-American', 'Asian');
  levels(biodb$Ethnicity) <- c('Non-Hispanic', 'Hispanic');
  levels(biodb$GeneticAncestry) <- c('European', 'African', 'East Asian', 'Native American');
  levels(biodb$GeneticRiskCategory) <- c('Low', 'Normal', 'High');
  # levels(biodb$MRIResult) <- c('No Legion', 'Legion Found');
  levels(biodb$HighestPIRADS) <- c('No lesion', 'Very low', 'Low', 'Intermediate', 'High', 'Very high');
  # levels(biodb$BiopsyResult) <- c('Negative', 'Positive');
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

  attr(biodb$Weight, 'label') <- "Weight (kg)";
  attr(biodb$Height, 'label') <- "Height (cm)";
  attr(biodb$ProstateVolume, 'label') <- "Prostate Volume (cm^3)";
  attr(biodb$p2PSA, 'label') <- "[-2]proPSA";
  attr(biodb$freePSA, 'label') <- "free PSA";

  biodb$PHI.computed <- with(biodb, {
    (p2PSA / freePSA) * sqrt(PSAHyb)
  });

  biokey <- xlsx::read.xlsx(
    biomark.key.path,
    sheetIndex = 1,
    header = TRUE,
    stringsAsFactors = FALSE
    );

  biokey <- biokey[!is.na(biokey$Column.ID), ]

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

  if(onlyBiodb) {
    biodb
    } else {
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
    'MRI DOD Biomarkers Database_Boutros - 2020.10.28.xlsx',
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
  xlsx::read.xlsx(
    here::here(paste0('data/', file.name)),
    sheetIndex = 1,
    header = TRUE,
    stringsAsFactors = FALSE
  )
}
