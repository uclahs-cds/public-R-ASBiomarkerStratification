factor.ISUP <- function(x) {
  factor(x, levels = 0:5, ordered = TRUE)
  }

factor.Gleason <- function(x) {
  factor(x, levels = c("0+0", "3+3", "3+4", "4+3", "4+4", "5+5"), ordered = TRUE)
  }

mutation.dummy <- function(x) {
  mutation1 <- as.integer(x["Mutation1"]);
  mutation2 <- as.integer(x["Mutation.2"]);
  if(is.na(mutation1) && is.na(mutation2)) {
    # Missing value, make all dummy variables missing as well
    res <- list(Mutation.BRCA1 = NA, Mutation.BRCA2 = NA, Mutation.ATM = NA, Mutation.MLH1 = NA, Mutation.PMS2 = NA);
  }
  else {
    res <- list(Mutation.BRCA1 = 0, Mutation.BRCA2 = 0, Mutation.ATM = 0, Mutation.MLH1 = 0, Mutation.PMS2 = 0);
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
#' @export
load.data.AS <- function(biomark.path,
                         biomark.key.path,
                         genetics.path,
                         biomark.categories.path) {
  factor.cols <- c('Race', 'Ethnicity', 'MRIResult', 'HighestPIRADS',
                   'BiopsyResult','Observation',
                   'BiopsyUpgraded', 'GeneticAncestry', 'GeneticRiskCategory',
                   # 'Mutation1', 'Mutation.2', 'GlobalScreeningArray', 'BRCAMutation',
                   'RSIlesionPIRADS',
                   'RSIlesionCancer',  'RSIlesionUpgraded',
                   'RSIlesionObservation', 'ProgressedToTreatment', 'Prostatectomy',
                   'UpgradedAndProgressed', 'NoUpgradeAndProgressed'
                   );

  ISUP.cols <- c('PreviousISUP',  'RSIlesionISUP', 'StudyHighestISUP');
  Gleason.cols <- c('PreviousGleason', 'StudyHighestGleason', 'RSIlesionGleason');

  mutation.cols <- c("Mutation1", "Mutation.2");

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

  # Add PSA and PHI Density
  biodb$PSADensity <- biodb$freePSA / biodb$ProstateVolume;
  biodb$PHIDensity = biodb$PHI / biodb$ProstateVolume;

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

  attr(biodb$ProstateVolume, 'label') <- "Prostate Volume (cm^3)"

  biokey <- xlsx::read.xlsx(
    biomark.key.path,
    sheetIndex = 1,
    header = TRUE,
    stringsAsFactors = FALSE
    );

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

  # Return the loaded objects
  list(
    biodb = biodb,
    biokey = biokey,
    biokey.gen = biokey.gen,
    bio.categories = bio.categories
    );
  }

#' Loads the data for analysis with default file names
#' @export
default.load.data <- function() {
  file.names <- c(
    'MRI DOD Biomarkers Database_Boutros - 2020.04.20.xlsx',
    'MRI DOD Biomarkers Database_Boutros - Key.xlsx',
    'MRI DOD Biomarkers Genetics_Boutros - 2019.11.12.xlsx',
    'MRI DOD Biomarkers Category_Boutros.xlsx'
    );
  file.paths <- here::here(paste0('data/', file.names));
  do.call('load.data.AS', as.list(file.paths));
  }


