
# ProducDiv
# Fitting Bayesian GLMMs

require(tidytable)
require(brms)


# Model 1: Spatial ----------------------------------------------------------------------------

# weakly informative priors?
default_prior(prod_kg_yr ~ S + biogeo + hgravity + (1|DataID) + (1|gridcell), 
              data = vars_single, family = Gamma(link = 'log'))

mod_spatial <- brm(data = vars_single, 
                   formula = prod_kg_yr ~ S + (S|biogeo) + hgravity + (1|DataID) + (1|gridcell),
                   family = Gamma(link = 'log'), 
                   prior = c(
                     prior(normal(0.5,3), class = b, coef = S), # richness effect
                     prior(normal(-0.1, 1), class = 'b', coef = 'hgravity') # human gravity
                   ))
