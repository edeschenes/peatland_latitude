![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7094810.svg)

# Mosses and vascular plants biodiversity patterns along a latitudinal gradient in peatlands
This repository includes the data and R scripts to reproduce analyses and figures found in the article Mosses and vascular plants show diverging taxonomic and functional biodiversity patterns along a latitudinal gradient in bogs and fens by Élise Deschênes, Monique Poulin, Marie-Hélène Brice, Pierre Legendre & Stéphanie Pellerin for publication in Journal of Biogeography. 

The analyses were carried out with R version 4.2.2  (a free software environment for statistical computing and graphics) and require the installation of a recent version of it.

**The following packages must be installed to run the scripts:**

- ggeffects
- effects
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
- Reshape2
- ggpmisc


**Below are the R commands to install them all:**

```
install.packages(
  c("ggeffects","effects", "nlme", "emmeans", "ggplot2", "dplyr", "car", "gridExtra", "svglite",
  "FD", "adespatial", "vegan", "eulerr", "ggcorrplot", "GGally", "reshape2", "ggpmisc")
)
```

To reproduce the entire analysis including data cleaning, analyses and figures, run:
```
source("scripts/Figure2_Alpha.R")
source("scripts/Figure3_LCBD.R")
source("scripts/Figures_4_5_Composition.R")
source("scripts/Figures_supp.R")
```
All data used for the analyses can be found in the data folder.

Scripts Figure2_Alpha.R, Figure3_LCBD.R, Figures_4_5_Composition.R and Figures_supp.R performed all the analyses and produced the figures.


