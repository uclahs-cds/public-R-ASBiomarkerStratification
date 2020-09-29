# Notes

- Look at the load data file here:
/data/users/alfonsolam/projects/ProstateCancer-ASBiomarkerSynergy/analysis

## Project

There is a clinical need to predict *before* surgery if a man has aggressive disease so that we can decide if they need surgery at all.  The quality of life benefits are thus huge (avoiding therapy entirely!).  The problem is that because it’s pre-surgery, we do not have the full cancer to study.  Instead we use biopsies (spatially-restricted samples of the cancer), radiology (imaging like MRI) and minimally-invasive biomarkers (e.g. blood or urine tests).  It's unclear which of those different strategies is best, and how those strategies should be sequenced or ordered.  A collaborator at UTHSCSA (University of Texas Health Sciences Center San Antonio, I think) Dr. Michael Liss is a surgeon who is thinking hard about these problems.  He's put together a really nice cohort of ~100 patients where basically every possible biomarker has been generated and we want to figure out 'what is the optimal biomarker we can make using all tests'.  That will put an upper-bound to accuracy which we can go investigate in a prospective clinical trial.  We can then go backwards to start figuring out if there are ways to simplify/cheapen that test.

  - To revisit the initial descriptive analysis (feel free to keep any/all code you find useful, but no requirement for that)
  - Do some careful CV to identify the optimal model and operating point (probably F1 for the latter)
  - Take a look at full and truncated time-to-event for that final operating point
  

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
    
## Background

  - The Prostate Health Index (PHI): a new test for the detection of prostate cancer https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3943368
    - In multiple prospective international trials, this composite measurement has been shown to outperform conventional PSA and free PSA measurements.
    - Unlike PCA3 and TMPRSS2:ERG, PHI is also consistently associated with Gleason score and upgrading during active surveillance.

The Prostate Health Index (PHI) is computed as
$$
\text{PHI} = ([-2]\text{proPSA}/\text{free PSA}) \times \sqrt{\text{PSA}}
$$

In the data set this corresponds to `(p2PSA / freePSA) * sqrt(PSAHyb)`

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
    
