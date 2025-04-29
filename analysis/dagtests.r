
# ProducDiv
# DAG tests

require(dagitty)
require(tidytable)
require(stringr)
require(comets)
require(dHSIC)

# Define DAGs ======================================================================================

dag <- dagitty('dag {
bb="-3.522,-4.515,3.23,3.63"
"bioavailable energy" [latent,pos="0.443,-3.061"]
"complementarity/competition" [latent,pos="0.973,-1.548"]
"fishing pressure" [adjusted,pos="-0.553,-0.007"]
"size structure" [latent,pos="-0.016,1.560"]
"species diversity" [exposure,pos="0.979,1.547"]
"structural complexity" [latent,pos="1.880,1.540"]
"total abundance" [pos="0.473,-0.027"]
"total biomass" [pos="-1.052,-1.521"]
"unknown abiotic" [latent,pos="0.176,-2.258"]
SST [pos="-0.569,-3.061"]
biogeography [pos="1.506,-0.027"]
productivity [outcome,pos="-0.052,-1.521"]
"bioavailable energy" -> "complementarity/competition"
"bioavailable energy" -> "total abundance"
"bioavailable energy" -> "total biomass"
"complementarity/competition" -> productivity
"fishing pressure" -> "size structure"
"fishing pressure" -> "species diversity"
"fishing pressure" -> productivity
"size structure" -> productivity
"species diversity" -> "complementarity/competition"
"structural complexity" -> "species diversity"
"total abundance" -> "complementarity/competition"
"total abundance" -> "size structure"
"total abundance" -> "total biomass"
"total abundance" -> productivity
"unknown abiotic" -> "bioavailable energy"
"unknown abiotic" -> "complementarity/competition"
"unknown abiotic" -> SST
"unknown abiotic" -> productivity
SST -> productivity
biogeography -> "species diversity"
productivity -> "total biomass"
}
')

dag_t <- dagitty('dag {
bb="-4.873,-5.187,5.286,5.164"
"climate change" [latent,pos="3.383,-1.525"]
"coral cover" [latent,pos="-2.918,-1.502"]
"fishing pressure" [adjusted,pos="-2.089,0.954"]
"local extinction" [latent,pos="-1.050,-2.574"]
"range shifts" [latent,pos="0.774,-4.729"]
"size structure" [latent,pos="-0.928,2.439"]
"∆N" [adjusted,pos="1.105,3.109"]
"∆Prod" [outcome,pos="1.594,-0.843"]
"∆S" [exposure,pos="-1.045,-0.837"]
"∆SST_alpha" [adjusted,pos="2.816,1.301"]
"∆SST_gamma" [latent,pos="2.231,-3.289"]
emigration [latent,pos="0.425,-2.775"]
immigration [latent,pos="-2.054,-3.177"]
"climate change" -> "∆SST_alpha"
"climate change" -> "∆SST_gamma"
"coral cover" -> "∆S"
"fishing pressure" -> "size structure"
"fishing pressure" -> "∆S"
"local extinction" -> "∆N"
"local extinction" -> "∆S"
"range shifts" -> "local extinction"
"range shifts" -> emigration
"range shifts" -> immigration
"size structure" -> "∆Prod"
"∆N" -> "∆Prod"
"∆S" -> "∆Prod"
"∆SST_alpha" -> "∆Prod"
"∆SST_gamma" -> "range shifts"
emigration -> "∆S"
immigration -> "∆S"
}
')

# Some basic tests

isAcyclic(dag)
isAcyclic(dag_t)

adjustmentSets(dag, exposure = 'species diversity', outcome = 'productivity',
               type = 'minimal') # yep just fishing pressure
adjustmentSets(dag_t, exposure = '∆S', outcome = '∆Prod',
               type = 'minimal') # yep

impliedConditionalIndependencies(dag)
impliedConditionalIndependencies(dag_t)


# Prep variables ===================================================================================

# Some variables won't be in the model but need to be calculated for these data consistency checks

vars_single <- dt_single %>% 
  # total biomass
  left_join(., dt_traits[c(1, 4, 5)], by = 'taxaID') %>%
  mutate(biomass = Abundance * a * Length ^ b/1000) %>% 
  group_by(seriesID, SamplingEventID, Year, Month, date_res) %>% 
  summarise(totalbiomass = sum(biomass, na.rm = T), totalN = sum(Abundance)) %>% 
  select(seriesID, SamplingEventID, Month, Year, date_res, starts_with('total')) %>% 
  # total biomass + abundance
  left_join(vars_single, ., by = c('seriesID', 'SamplingEventID', 'Year')) %>%
  # SST
  left_join(., dt_sst, by = c('gridcell', 'Year', 'Month', 'date_res'))
summary(vars_single)
vars_single <- vars_single %>% select(!c(Month, date_res))

temp <- dt_rareass %>% group_by(DataID, seriesID, SamplingEventID, Year) %>% 
  summarise(totalN = sum(Abundance)) %>% ungroup

# full time series ∆N
mod_Nt <- as.list(rep(NA, nrow(vars_time_long)))
for (i in 1:nrow(vars_time_long)) {
  mod_Nt[[i]] <- temp %>% filter(DataID == vars_time_long$DataID[i],
                                 seriesID == vars_time_long$seriesID[i]) %>% 
    glm(data = .,
        totalN ~ Year, family = Gamma(link = 'log'))
}

vars_time_long$N_slope <- NA
for (i in 1:nrow(vars_time_long)) {
  vars_time_long$N_slope[i] <- coefficients(mod_Nt[[i]])[2]
}

# Data consistency =================================================================================

# Conditional independencies to check for

## og DAG ==========================================================================================

# fishing pressure ⊥ total abundance
# fishing pressure ⊥ SST
# fishing pressure ⊥ biogeography
# species diversity ⊥ total abundance
# species diversity ⊥ SST
# total abundance ⊥ biogeography
# total biomass ⊥ biogeography | fishing pressure, species diversity
# SST ⊥ biogeography
# biogeography ⊥ productivity | fishing pressure, species diversity

# Y ⊥ X | Z: y independent of x given z
# beware of high rmseas

vars_single %>% mutate(hgrav2 = log(hgrav2 + 1), totalN = log(totalN)) %>% 
  ciTest(Y = 'hgrav2', X = 'totalN', data = .)
vars_single %>% mutate(hgrav2 = log(hgrav2 + 1), biogeo = as.factor(biogeo)) %>% 
  ciTest(Y = 'hgrav2', X = 'biogeo', data = ., type = 'cis.pillai')
vars_single %>% mutate(hgrav2 = log(hgrav2 + 1), S = log(S)) %>% 
  ciTest(Y = 'S', X = 'totalN', data = .)
vars_single %>% mutate(S = log(S)) %>% 
  ciTest(Y = 'S', X = 'sst', data = .)
ciTest(Y = 'totalN', X = 'biogeo', data = vars_single, type = 'cis.chisq')
vars_single$biogeo <- as.factor(vars_single$biogeo)
ciTest(Y = 'biogeo', X = 'prod_kg', Z = c('hgrav2','S'), data = vars_single, type ='cis.pillai')
comet(biogeo ~ prod_kg | hgrav2 + S, data = vars_single, test = 'gcm')
ciTest(Y = 'sst', X = 'biogeo', data = vars_single, type = 'cis.chisq')
comet(totalbiomass ~ biogeo | hgrav2 + S, data = vars_single, test = 'gcm')

## Temporal DAG ====================================================================================

# fishing pressure ⊥ ∆SST_alpha
# fishing pressure ⊥ ∆N

ciTest(Y = 'hgrav2', X = 'sst_slope', data = vars_time_long)
ciTest(Y = 'hgrav2', X = 'N_slope', data = vars_time_long)

# tests passed
# final adjustments:
# OG DAG: fishing pressure
# temporal DAG: ∆SST, ∆N, fishing pressure