
# ProducDiv
# Fitting Bayesian GLMMs

require(tidytable)
require(brms)

# Model 2: Temporal, full time series ---------------------------------------------------------
# does proportional diversity change result in proportional productivity change?

tlong_priors <- default_prior(prod_slope ~ S_slope + hgravity + (1|DataID) + (1|gridcell),
                              data = vars_time_long, family = gaussian())
tlong_priors
tlong_priors$prior[2] <- 'normal(-0.1, 3)' # human gravity, slight negative
tlong_priors$prior[3] <- 'normal(0.1, 3)' # richness effect, slight positive
set.seed(2)
mod_tlong <- brm(data = vars_time_long,
                 formula = prod_slope ~ S_slope + hgravity + (1|DataID) + (1|gridcell),
                 family = gaussian(),
                 prior = tlong_priors, iter = 3000)
save(mod_tlong, file = 'mod_tlong.RData')

# Model 1: Spatial ----------------------------------------------------------------------------
# what is the diversity-productivity relationship like across space

vars_single$hgravity <- vars_single$hgravity/500

# weakly informative priors
priors <- default_prior(prod_kg_yr ~ log(S) + hgravity + (1|DataID) + (1|gridcell),
              data = vars_single, family = Gamma(link = 'log'))
priors
priors$prior[3] <- 'normal(0.1, 3)' # richness effect, slight positive
priors$prior[2] <- 'normal(-0.1, 3)' # human gravity, slight negative

set.seed(1)
mod_spatial <- brm(data = vars_single,
                   formula = prod_kg_yr ~ log(S) + hgravity + (1|DataID) + (1|gridcell),
                   family = Gamma(link = 'log'),
                   prior = priors, iter = 3000)
save(mod_spatial, file = 'mod_spatial.RData')

