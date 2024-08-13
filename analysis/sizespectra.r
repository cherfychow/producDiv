
# ProducDiv
# Size spectra slopes

require(dplyr)
require(purrr)

# import neg-log-likelihood function from Edwards et al. 2017 for size spectra slopes
negLL.PLB = function(b, x, n, xmin, xmax, sumlogx)
  {
    if(xmin <= 0 | xmin >= xmax) stop("Parameters out of bounds in negLL.PLB")
    if(b != -1)
      { neglogLL = -n * log( ( b + 1) / (xmax^(b + 1) - xmin^(b + 1)) ) -
            b * sumlogx
      } else
      { neglogLL = n * log( log(xmax) - log(xmin) ) + sumlogx
      }
    return(neglogLL)
  }

# starting values for MLE
# b: value of b for which to calculate the negative log-likelihood
#   x: values of data x
#   n: length(x), have as an input to avoid repeatedly calculating it 
#   xmin: minimum value of x, have as an input to avoid repeatedly calculating
#   xmax: maximum value of x, have as an input to avoid repeatedly calculating
#   sumlogx: sum(log(x)) as an input, to avoid repeatedly calculating

b = 1/( log(min(x)) - sum.log.x/length(x)) - 1

size_start <- dt_rareass %>% select(rareID, SamplingEventID, Length, Abundance)
# repeat each length by its abundance
size_start <- size_start[rep(1:nrow(size_start), each = size_start$Abundance),] %>%
  select(!Abundance) %>% # collapse lengths
  group_by_at(vars(-Length)) %>% summarise(lengths = paste(Length, collapse = ','))

PLB.minLL =  nlm(negLL.PLB, p = PL.bMLE, x = x, n = length(x),
    xmin=xmin, xmax=xmax, sumlogx=sum.log.x) #, print.level=2 )

PLB.bMLE = PLB.minLL$estimate

# 95% confidence intervals for MLE method.

PLB.minNegLL = PLB.minLL$minimum