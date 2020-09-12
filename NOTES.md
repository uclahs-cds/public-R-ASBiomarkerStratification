# Private Notes

- Look at the load data file here:
/data/users/alfonsolam/projects/ProstateCancer-ASBiomarkerSynergy/analysis

## Questions
  - Predicting `BiopsyUpgraded` and `ProgressedToTreatment`?
  - What are we optimizing? Prediction? Cost?
  - What is the difference between SOCPSA and PSAHyb? Looks like PHI is computed from PSAHyb
    - Generally correlated but a few outliers (one extreme one that could be input error)
  - Is there a reason `DaysBxToUpgrade` is missing?

## References
  - The Prostate Health Index (PHI): a new test for the detection of prostate cancer https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3943368
    - In multiple prospective international trials, this composite measurement has been shown to outperform conventional PSA and free PSA measurements.
    - Unlike PCA3 and TMPRSS2:ERG, PHI is also consistently associated with Gleason score and upgrading during active surveillance.

## Data
  - PHI = ([-2]proPSA/free PSA) × √PSA.
    - PHI = (p2PSA / freePSA) * sqrt(PSAHyb)
    - One subject M0003 has an incorrect value of 1561
    
## Terms
  - Magnetic resonance imaging (MRI)
    - Sensitivity of MRI for the detection of extracapsular extension has been reported to range from 13% to 95% (!).
  - Restriction spectrum imaging (RSI)
  - Prostate Health Index (PSI)
    - The Prostate Health Index (PHI) is a new formula that combines all three forms (total PSA, free PSA and p2PSA) into a single score that can be used to aid in clinical decision-making

## Methods
  - From Alfono's work the way of creating the data for survival analysis is to:
    1. Use the `DaysDxToUpgrade` if the subject was upgraded (event = 1)
    2. Otherwise, use `DaysBxToLastReview` (event = 0)
  - Seems like it would make more sense to use both `Dx` (Diagnosis) or `Bx` (Biopsy) for each instead
    
