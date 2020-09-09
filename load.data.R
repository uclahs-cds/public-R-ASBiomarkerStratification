
# Project: Prostate Cancer Active Surveillance

#######################################################
# Loading Library:
#######################################################
library(BoutrosLab.utilities);
library(BoutrosLab.plotting.general);
library(xlsx);
library(here);

#######################################################
# Put input data into objects:
#######################################################

path <- here('data');

# Latest dataset 04/20/2020:
bio.file <- paste(
  path,
  'MRI DOD Biomarkers Database_Boutros - 2020.04.20.xlsx',
  sep = '/'
  );

biokey.file <- paste(
  path,
  'MRI DOD Biomarkers Database_Boutros - Key.xlsx',
  sep = '/'
  );

biokey.gen <- paste(
  path,
  'MRI DOD Biomarkers Genetics_Boutros - 2019.11.12.xlsx',
  sep = '/'
  );

# Open Raw data:
biodb <- read.xlsx(
  bio.file,
  1,
  header = TRUE
  );
biodb <- biodb[!apply(is.na(biodb), 2, all)];

biokey <- read.xlsx(
  biokey.file,
  1,
  header = TRUE
  );
biokey <- biokey[! apply(is.na(biokey), 2, all)];

biokey.gen <- read.xlsx(
  biokey.gen,
  1,
  header = TRUE
  );
biokey.gen <- biokey.gen[! apply(is.na(biokey.gen), 2, all)];

#######################################################
# Sort data by Age, Race and Ethnicity:
#######################################################

sorted.biodb <- biodb[order(biodb$Age, biodb$Race, biodb$Ethnicity),]

### WRITE SESSION PROFILE TO FILE #####################
save.session.profile(
  BoutrosLab.utilities::generate.filename(
    Sys.Date(),
    'grab.data.info',
    'txt'
    )
  );
