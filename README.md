# Bryophytes and vascular plants biodiversity patterns along a latitudinal gradient in peatlands
This repository includes the data and R scripts to reproduce analyses and figures found in the article Bryophytes and vascular plants show diverging taxonomic and functional biodiversity patterns along a latitudinal gradient in bogs and fens by Élise Deschênes, Monique Poulin, Marie-Hélène Brice, Pierre Legendre & Stéphanie Pellerin for publication in Journal of Biogeography. 

The analyses were carried out with R version 4.2.2  (a free software environment for statistical computing and graphics) and require the installation of a recent version of it.

**The following packages must be installed to run the scripts:**

- ggeffects
- nlme
- emmeans
- ggplot2
- dplyr
- car
- gridExtra
- svglite
- FD
- adespatial
- vegan
- eulerr
- ggcorrplot
- GGally

**Below are the R commands to install them all:**

```
install.packages(
  c("ggeffects", "nlme", "emmeans", "ggplot2", "dplyr", "car", "gridExtra", "svglite",
  "FD", "adespatial", "vegan", "eulerr", "ggcorrplot", "GGally")
)
```

To reproduce the entire analysis including data cleaning, analyses and figures, run:
```
source("scripts/Figure1_Alpha.R")
source("scripts/Figure2_Beta.R")
source("scripts/Figures3_4_Composition.R")
source("scripts/FigureS1_S2_S3.R")
```
All data used for the analyses can be found in the data folder.

Scripts Figure1_Alpha.R, Figure2_Beta.R, Figures3_4_Composition.R and FigureS1_S2_S3.R performed all the analyses and produced the figures.
