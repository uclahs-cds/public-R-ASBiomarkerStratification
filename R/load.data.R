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
  factor.cols <- c('Race', 'Ethnicity', 'MRIResult', 'MRILesions', 'HighestPIRADS',
                   'BiopsyResult', 'PreviousGleason', 'PreviousISUP', 'Observation',
                   'BiopsyUpgraded', 'GeneticAncestry', 'GeneticRiskCategory', 'GlobalScreeningArray',
                   'GSAPositives', 'BRCAMutation', 'Mutation1', 'Mutation.2', 'RSIlesionPIRADS',
                   'RSIlesionCancer', 'RSIlesionGleason', 'RSIlesionISUP', 'RSIlesionUpgraded',
                   'RSIlesionObservation', 'ProgressedToTreatment', 'Prostatectomy',
                   'UpgradedAndProgressed', 'NoUpgradeAndProgressed')
  # Open Raw data:
  biodb <- xlsx::read.xlsx(
    biomark.path,
    sheetIndex = 1,
    header = TRUE
  );

  # Convert to factor
  biodb[,factor.cols] <- lapply(biodb[,factor.cols], as.factor)
  biodb <- biodb[!apply(is.na(biodb), 2, all)];

  # Update levels
  levels(biodb$Race) <- c('White', 'African-American', 'Asian');
  levels(biodb$Ethnicity) <- c('Non-Hispanic', 'Hispanic');
  levels(biodb$GeneticAncestry) <- c('European', 'African', 'East Asian', 'Native American');
  levels(biodb$GeneticRiskCategory) <- c('Low', 'Normal', 'High');

  biokey <- xlsx::read.xlsx(
    biomark.key.path,
    sheetIndex = 1,
    header = TRUE
    );

  biokey.gen <- xlsx::read.xlsx(
    genetics.path,
    sheetIndex = 1,
    header = TRUE
    );

  bio.categories <- xlsx::read.xlsx(
    biomark.categories.path,
    sheetIndex = 1,
    header = TRUE
    );

  bio.categories$Category <- as.factor(bio.categories$Category);

  #######################################################
  # Sort data by Age, Race and Ethnicity:
  #######################################################

  sorted.biodb <- biodb[order(biodb$Age, biodb$Race, biodb$Ethnicity),]

  # Return the loaded objects
  list(
    biodb = biodb,
    biokey = biokey,
    biokey.gen = biokey.gen,
    sorted.biodb = sorted.biodb,
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


