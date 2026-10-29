

# ProducDiv
# Calculate asymptotic diversity and final prep of diversity predictors

require(iNEXT.3D)
require(tidytable)
require(dplyr)
require(ggplot2)

# requires 01_gridding for the data source

# CALCULATE OBSERVED RICHNESS ======================================================================

dt_div <- dt_sp %>% select(taxaID, worms_valid) %>% 
  left_join(dt_rareass, ., by = 'taxaID') %>% 
  group_by(DataID, seriesID, SamplingEventID, Year) %>% 
  summarise(S = n_distinct(worms_valid)) %>% ungroup


# CALCULATE SAMPLE COMPLETENESS ====================================================================

# from iNEXT by Hsieh et al. 2016
# @param x a vector of species abundances
# @param m a integer vector of rarefaction/extrapolation sample size
# @return a vector of estimated sample coverage function
inext_compl <- function(x, m){
  x <- x[x>0]
  n <- sum(x)
  f1 <- sum(x == 1)
  f2 <- sum(x == 2)
  f0.hat <- ifelse(f2 == 0, (n - 1) / n * f1 * (f1 - 1) / 2, (n - 1) / n * f1 ^ 2/ 2 / f2)  #estimation of unseen species via Chao1
  A <- ifelse(f1>0, n*f0.hat/(n*f0.hat+f1), 1)
  Sub <- function(m){
    #if(m < n) out <- 1-sum(x / n * exp(lchoose(n - x, m)-lchoose(n - 1, m)))
    if(m < n) {
      xx <- x[(n-x)>=m]
      out <- 1-sum(xx / n * exp(lgamma(n-xx+1)-lgamma(n-xx-m+1)-lgamma(n)+lgamma(n-m)))
    }
    if(m == n) out <- 1-f1/n*A
    if(m > n) out <- 1-f1/n*A^(m-n+1)
    out
  }
  sapply(m, Sub)		
}

# key index to start off the output table
dt_complete <- dt_single %>% select(DataID, seriesID, SamplingEventID, Abundance) %>% 
  group_by(DataID, seriesID, SamplingEventID) %>% summarise(totalN = sum(Abundance))
dt_complete$compness <- NA

for (i in 1:nrow(dt_complete)) {
  # turn densities back into counts
 tempdt <- dt_single %>% 
    filter(DataID == dt_complete$DataID[i],
           seriesID == dt_complete$seriesID[i],
           SamplingEventID == dt_complete$SamplingEventID[i]) %>% 
    group_by(worms_valid) %>% summarise(Abundance = sum(Abundance))
 scale <- 1/min(tempdt$Abundance)
 tempdt <- tempdt$Abundance * scale
 
  dt_complete$compness[i] <- inext_compl(x = tempdt, m = dt_complete$totalN[i])
}
rm(tempdt, scale)