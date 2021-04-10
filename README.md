# Comparing Non-Invasive Biomarkers for Prediction of Prostate Cancer Progression on Active Surveillance

## Introduction

Widespread use of serum prostate-specific antigen (PSA) screening results in 50% of new cases of prostate cancer being diagnosed with low-grade localized disease.1-4 Standard of care for these cases is to defer immediate treatment in favor of active surveillance (AS), a low-toxicity management which entails close monitoring with PSA tests, repeat prostate biopsies and multi-parametric MRI (mpMRI).5-8 Around 30% of men on AS progress or electively undergo definitive treatment within two years of diagnosis.8,9 Efficacious low-cost, low-toxicity strategies to limit exit from AS are needed to reduce over-treatment and improve health-related quality of life (QOL).

## Building Package
Using devtools at the root of the directory
```{r}
devtools::install()
```
or from the command line
```{r}
R CMD INSTALL --no-multiarch --with-keep.source project-ProstateCancer-ASBiomarkerSynergy
```

## Models and results
The models comparing clinico-epidemiologic, radiologic, and molecular features is found in [inst/category_models.R](inst/category_models.R)
The sequential models were built using [inst/sequential_models.R](inst/sequential_models.R).
The RSI and PI-RADS comparison models are found here: [inst/rsi_baseline_models.R](inst/rsi_baseline_models.R)
.
