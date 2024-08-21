
# ProducDiv
# Size spectra slopes

require(dplyr)
require(stringr)
require(beepr)

# SPECIFY FUNCTIONS ================================================================================

# adapted from Edwards et al. 2017: https://doi.org/10.1111/2041-210X.12641

# import neg-log-likelihood function from Edwards et al. 2017 for size spectra slopes
# starting values for MLE
  # b: value of b for which to calculate the negative log-likelihood
  # x: values of data x
  # xmin: minimum value of x, have as an input to avoid repeatedly calculating
  # xmax: maximum value of x, have as an input to avoid repeatedly calculating
spectra_negL = function(b, x, xmin, xmax) {
  n = length(x)
  if(xmin <= 0 | xmin >= xmax) stop("Parameters out of bounds for size spectra estimation")
  if(b != -1)
    { neglogLL = -n * log( ( b + 1) / (xmax^(b + 1) - xmin^(b + 1)) ) -
          b * sum(log(x))
    } else
    { neglogLL = n * log( log(xmax) - log(xmin) ) + sum(log(x))
    }
  return(neglogLL)
}

# function that gives lower and upper 95% confidence intervals
spectra_conf = function(x, xmin = 1, xmax = 100, mle) {
  # x = length data as a numeric vector
  # mle = object name of estimate output from MLE optimisation

  bvec = seq(mle$estimate - 0.5, mle$estimate + 0.5, 0.00001) 
      
  PLB.LLvals = vector(length=length(bvec))  # negative log-likelihood for bvec
  for(i in 1:length(bvec))
      {
          PLB.LLvals[i] = spectra_negL(bvec[i], x = x, xmin = xmin, xmax = xmax)   
      }
      
  critVal = mle$minimum + qchisq(0.95,1)/2
                      # 1 degree of freedom, Hilborn and Mangel (1997) p162.
  bIn95 = bvec[ PLB.LLvals < critVal ]
                      # b values in 95% confidence interval
  PLB.MLE.bConf = c(min(bIn95), max(bIn95))
  return(PLB.MLE.bConf)
}

# all inclusive
# returns a vector of estimate, lwr, upper
size_spectra = function(x) {

  xmin = min(x)
  xmax= max(x)

  # calculate starting b value
  b_start  <-  1/( log(min(x)) - sum(log(x)) / length(x)) - 1

  # instead of specifying b, we give the starting value as p to get the opimisation started
  spectra_MLE <-  nlm(spectra_negL, p = b_start, x = x, xmin = xmin, xmax = xmax)
  
  # calculate confidence intervals
  conf <- spectra_conf(x = x, mle = spectra_MLE, xmin = xmin, xmax = xmax)

  output  <- c(spectra_MLE$estimate, conf)
  names(output) <- c('estimate', 'lwr', 'upr')
  return(output)
}

# LOOP ALL SAMPLING EVENTS =========================================================================

# make a table object that we can loop starting values in
dt_size <- dt_rareass %>% select(rareID, SamplingEventID, Length, Abundance)
# repeat each length by its abundance
dt_size <- dt_size[rep(1:nrow(dt_size), times = dt_size$Abundance),] %>%
  select(!Abundance) %>% # collapse lengths
  group_by_at(vars(-Length)) %>% summarise(lengths = paste(Length, collapse = ','))

# make sure this holds up
n_distinct(dt_size$SamplingEventID) == nrow(dt_size)

lengths <- as.list(NA)
for (i in 1:nrow(dt_size)) {
  lengths[[i]] <- str_split(dt_size$lengths[i], ',', simplify = T) %>% as.numeric
}

dt_size$lengths <- NULL
dt_size$b <- NA
dt_size$b_lwr <- NA
dt_size$b_upr <- NA

for (i in 1:nrow(dt_size)) {

  spec <- tryCatch(size_spectra(x = lengths[[i]]), error = function(e) e)
  if(inherits(spec, "error")) next # error handling
  
  dt_size$b[i] <- spec[1]
  dt_size$b_lwr[i] <- spec[2]
  dt_size$b_upr[i] <- spec[3]

}
beep()

# make a unified dataframe that contains all the variables for each sampling event
dt_variables <- dt_div %>% filter(metric == 'S') %>% 
tidyr::pivot_wider(names_from = 'metric', values_from = 'value') %>%
select(rareID, SamplingEventID, S) %>%
left_join(., dt_prod, by = c('rareID', 'SamplingEventID')) %>%
left_join(., dt_size, by = c('rareID', 'SamplingEventID')) %>%
select(!Area) %>%
relocate(DataID, rareID, SamplingEventID, Year)
