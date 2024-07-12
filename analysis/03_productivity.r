
# ProducDiv
# Calculate productivity

require(dplyr)
require(rfishprod)
require(rfishbase)

iter = 10

# Check dataset repdata #
(repdata <- rfishprod:::repdata)

# Getting levels ready #
repdata <- tidytrait(repdata, db)

# Formula from Morais and Bellwood (2018) #
fmod <- formula(~ sstmean + MaxSizeTL + Diet + Position + Method)
# estimate VGBF Kmax
dt_kmax <- predKmax(repdata, dataset = db,
                   fmod = fmod, niter = iter,
                   return = 'pred')

