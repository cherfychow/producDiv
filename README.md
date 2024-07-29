# ProducDiv
## Quantifying global reef fish productivity and its links with assemblage biodiversity

---

#### Repository structure

- *`assembly_wrangling`*: individual data standardisation scripts for fish assemblage data
- *`dataset_outputs`*: individual standardised data
- *`analysis`*: analysis pipeline R scripts
- *`src`*: source data for analysis pipeline


### Analysis source data (`src`)

- *`fish_prodtraits.csv`*: compiled traits on water column position and reef association, maximum length (TL), diet (sensu [Parravicini et al. 2020](https://doi.org/10.1371/journal.pbio.3000702)), and Bayesian length-weight parameters for all observed fishes in the dataset.
- *`gaspar_traits.csv`*: traits compiled by [Floeter et al. 2017](https://doi.org/10.1111/brv.12336) that we use to partially compile water column position traits
- *`missing_traits_genus.csv`*: genus-level traits compiled manually from FishBase model estimates and most common trait values.
- *`missing_troph_sp.csv`*: diet traits for species that don't have machine readable diet information on FishBase
- *`oisst_landmask.nc`*: Spatial data to filter out land areas from the NOAA OISST v2 rasters
- *`parravicini_trophic_guilds.csv`*: species diet traits from model extrapolations by [Parravicini et al. 2020](https://doi.org/10.1371/journal.pbio.3000702)