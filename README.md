# ProducDiv
## Reef fish biodiversity promotes productivity

This repository contains the data and code used for the analysis and figure generation the ProducDiv project. The project takes a macroecological approach to assess biodiversity-productivity links in global reef fishes across space and time. 

Disclaimer: This repository may still undergo revisions and changes.

**Corresponding author**: Cher Chow, cherchow (at) proton.me

#### *In review*

---

### Requirements

The analysis pipeline in this repo runs in an R environment (version 4.3.3). The code's been run in both Linux and Mac operating systems, but not trialed in PCs. Here are the following package dependencies:

|      Package         |      Version     |
|----------------------|------------------|
|     brms             |     2.22.0       |
|     dggridR          |     3.1.0        |
|     dplyr            |     1.1.4        |
|     ggdist           |     3.3.2        |
|     ggplot2          |     3.5.2        |
|     iNEXT            |     3.0.1        |
|     lme4             |     1.1.37       |
|     lubridate        |     1.9.4        |
|     patchwork        |     1.3.1        |
|     quantreg         |     6.1          |
|     readr            |     2.1.5        |
|     rfishbase        |     5.0.0        |
|     rfishprod        |     0.0.3        |
|     rnaturalearth    |     1.0.1        |
|     sf               |     1.0.16       |
|     stringr          |     1.5.1        |
|     tidybayes        |     3.0.7        |
|     tidync           |     0.3.0        |
|     tidytable        |     0.11.2       |
|     xgboost          |     1.7.7.1      |

### Repository structure

- **`analysis`**: analysis pipeline R scripts
- **`src`**: source and derived data
- `load_data.r`: useful script to load data required across the pipeline. It calls model output R objects, but these aren't included in the repo here due to file size restrictions. Please contact authors if needed.

### Data

  - `biogeo.shp`: shapefile of global EEZs assigned biogeographic region sensu Kulbicki et al. 2014
  - `dt_completeness.csv`: data on sample completeness
  - `dt_protections.csv`: data on area-based no-take restrictions
  - `dt_tchange.csv`: data on S:Prod time-series links
  - `gaspar_traits.csv`: traits compiled by [Floeter et al. 2017](https://doi.org/10.1111/brv.12336) that we use to partially compile water column position traits
  - `parravicini_trophic_guilds.csv`: species diet traits from model extrapolations by [Parravicini et al. 2020](https://doi.org/10.1371/journal.pbio.3000702)
  - `traits_final.csv`: final compiled file of fish traits
  - `vars_rare.csv`: compiled model variables for rarefied time-series (used to estimate temporal trends)
  - `vars_single.csv`: compiled model variables for spatial single-year data
  - `vars_time_flux.csv`: compiled model variables for adjacent year temporal trends
  - `vars_time.csv`: compiled model variables for temporal trends

### Analysis directory

  - `01_gridding.r`: spatial gridding and sample-based rarefaction to split data into single-year and time-series subsets
  - `02_diversity.r`: calculate species richness and sample completeness
  - `03_prodtrait_compile.r`: compile traits for calculating productivity
  - `04_sst_compile.r`: spatiotemporal match of SST to surveys 
  - `05_fishprod.r`: estimate Kmax -> productivity for all species populations (with traits + SST)
  - `06_predictors.r`: compile all model predictors together
  - `07_temporal.r`: model temporal trends for time-series surveys
  - `08_models.r`: Bayesian modelling of richness with productivity
  - `09_nullmodel.r`: Agnostic assemblage simulation with modelling
  - `10_sensitivity.r`: Length and sample completeness sensitivity tests
  - `dag_spatial_N-S.txt`: Dagitty code for spatial DAG (N -> S)
  - `dag_spatial_S-N.txt`: Dagitty code for spatial DAG (S -> N)
  - `dag_temporal.txt`: Dagitty code for temporal DAG
  - `figure_supp.r`: Supplementary figure generation
  - `figures.R`: Main text figure generation