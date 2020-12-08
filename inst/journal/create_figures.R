# Correlation heatmap
library(BoutrosLab.ASBiomarkerSynergy);
library(here)

biodb <- default.load.data(onlyBiodb = TRUE)

create.heatmap.AS(biodb,
                  filename = here('euro_urology/figures/corr_heatmap.tiff')
                  )
