# Notes

- Look at the load data file here:
/data/users/alfonsolam/projects/ProstateCancer-ASBiomarkerSynergy/analysis

## Project

There is a clinical need to predict *before* surgery if a man has aggressive disease so that we can decide if they need surgery at all.  The quality of life benefits are thus huge (avoiding therapy entirely!).  The problem is that because it’s pre-surgery, we do not have the full cancer to study.  Instead we use biopsies (spatially-restricted samples of the cancer), radiology (imaging like MRI) and minimally-invasive biomarkers (e.g. blood or urine tests).  It's unclear which of those different strategies is best, and how those strategies should be sequenced or ordered.  A collaborator at UTHSCSA (University of Texas Health Sciences Center San Antonio, I think) Dr. Michael Liss is a surgeon who is thinking hard about these problems.  He's put together a really nice cohort of ~100 patients where basically every possible biomarker has been generated and we want to figure out 'what is the optimal biomarker we can make using all tests'.  That will put an upper-bound to accuracy which we can go investigate in a prospective clinical trial.  We can then go backwards to start figuring out if there are ways to simplify/cheapen that test.

  - To revisit the initial descriptive analysis (feel free to keep any/all code you find useful, but no requirement for that)
  - Do some careful CV to identify the optimal model and operating point (probably F1 for the latter)
  - Take a look at full and truncated time-to-event for that final operating point
  
In general a patient with GS 3+3 (now increasingly called ISUP Grade Group 1 or ISUP GG1 or just GG1) will go on AS while a 3+4 (GG2) may but probably will not and a 4+3 (gg3) definitely will not.  Thus you can also look at the patients who did get surgery and see if a GG1 patient on biopsy (a tiny fraction of the tumor) was actually GG2 or GG3 (or worse) after surgery when the entire tumor is available for checking.  Thus looking at the surgical pathology can remove the spatial variability component (called undersampling) of the biopsy procedure out.

Clinically, AS means *not* treating a patient.  This is of course statistically superior for a bulk population, but still leads to problems when there is a FN (somebody has aggressive disease, but that isn’t recognized and thus they are inappropriately on AS).  Many patients will therefore voluntarily elect to exit AS prematurely, seeking a treatment that they may not derive benefit from.  Thus the key clinical problem is unrecognized aggressive disease, which both hinders those patients directly and more broadly reduces confidence in AS.  The goal of our work is thus to give better identification of which patients should exit AS, so that those who should not have more confidence in their decision.  Surgery/Radiotherapy (Sx/Rt) is extremely expensive ($20-50k) and has long-term morbidities (quality-of-life impacts).  We sometimes talk about this in terms of “health preference values”, which reflects how much a change in life-quality relates to a full-year of health.  Very roughly, every year on AS relative to definitive local therapy (Sx or Rt) saves the patient ~0.2 years of fully-healthy life, so prolonging AS as long as possible has big advantages.  For example, my dad has been on AS for almost 7 years now, and it’s very possible he will be on it for the rest of his life.  But, his initial immediate instinct was “cut this thing out of me”, and having a son who could walk him through the decision-making made a big difference there.  So in the long-run, the question is really “can these biomarkers improve confidence men have to remain on AS” as much as anything else -- ~15% each year *voluntarily* leave AS for treatment.
  

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
    
