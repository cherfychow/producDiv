
# ProducDiv
# Calculate productivity from fish traits + SST with observed lengths

require(dplyr)
require(stringr)
require(rfishprod)
require(beepr)

# PREP =========================================================================

# give predKmax data on all fishes observed
# then calculate productivity per rareID

# need to keep some sampling info since SST are specific to location + time
# dt_assemblage, because it's full data
# traits
prodtrait <- dt_assemblage %>% 
  distinct(DataID, gridcell, SamplingEventID, Year, Month, genus_sp, Length) %>%
  # all family level = just family
  mutate(genus_sp = case_when(str_detect(genus_sp, '(?<=idae)\\ssp\\d*$') ~ str_remove(genus_sp, '\\ssp\\d*$'),
                              .default = str_remove(genus_sp, '\\d+'))) %>%
  # known uncertains irrelevant for productivity
  left_join(., dt_traits, by = c('genus_sp' = 'Species')) # traits

# SST
prodtrait <- left_join(prodtrait, dt_sst[-c(1,5,6)], 
                       by = c('gridcell', 'Year', 'Month')) %>% # SST
  rename(sstmean = sst, Size = Length, Species = genus_sp)

# our own version of tidytrait
sort(unique(prodtrait$Diet)) == sort(unique(dt_prodDB$Diet))
sort(unique(prodtrait$Position)) == sort(unique(dt_prodDB$Position))
# have to fit the db's factor ordering
prodtrait <- prodtrait %>% 
  mutate(Method = 'Otolith',
         Diet = factor(Diet, levels = sort(unique(prodtrait$Diet))),
         Position = factor(Position, levels = levels(dt_prodDB$Position)))
summary(prodtrait)
prodtrait <- prodtrait %>% filter(!is.na(b)) # only 2 very uncertain records
# make db match
dt_prodDB$Diet <- factor(dt_prodDB$Diet, levels = sort(unique(prodtrait$Diet)))


# check if any fish are larger than their max size
prodtrait %>% filter(Size >= MaxSizeTL) %>% nrow

# limit observed size to be below max size
prodtrait <- mutate(prodtrait, 
                    Size = case_when(Size > MaxSizeTL ~ MaxSizeTL - 0.1, 
                                     .default = Size))

# Formula from Morais and Bellwood (2018) #
growthf <- formula(~ sstmean + MaxSizeTL + Diet + Position)

# Kmax, growth coef ============================================================

# Kmax estimate per species-gridcell-year-dataID / species-samplingeventid-dataID
# recommended boostrapts on 100s to 1000s scale
set.seed(24)
vbgf_par <- prodtrait %>% distinct %>% 
  predKmax(., dataset = dt_prodDB,
           fmod = growthf, niter = 100, return = 'pred')

vbgf_par <- vbgf_par$pred # just keep predictions
rm(prodtrait, growthf)

summary(vbgf_par)
# negative Kmax values
vbgf_par %>% filter(Kmax <= 0) %>% nrow() # work on this later... see if 1000 runs solves it

# Mortality and gross somatic gain =============================================

# mortality
vbgf_par <- vbgf_par %>% mutate(
  M = predM(Lmeas = Size, Lmax = MaxSizeTL, Kmax = Kmax, method = 'Lorenzen'),
  M_l = predM(Lmeas = Size, Lmax = MaxSizeTL, Kmax = Kmax_lowq, method = 'Lorenzen'),
  M_u = predM(Lmeas = Size, Lmax = MaxSizeTL, Kmax = Kmax_uppq, method = 'Lorenzen'))

vbgf_par <- vbgf_par %>% filter(Kmax > 0) %>% # leave out the nonsensicals for now
  mutate(
    Size_t1 = applyVBGF(Lmeas = Size, Lmax = MaxSizeTL, Kmax = Kmax), # rate of change day^-1
    size_delta = Size_t1 - Size
  ) # length change

# lowers and uppers refer to the Kmax value used!
# not necessarily the lower or upper of that metric
vbgf_par <- vbgf_par %>% 
  mutate(biomass_gain = somaGain(a, b, Lmeas = Size, Lmax = MaxSizeTL, Kmax = Kmax),
         biomass_loss = somaLoss(M, Lmeas = Size, a, b),
         delta_g_day = biomass_gain - biomass_loss,
         biomass_gain_l = somaGain(a, b, Lmeas = Size, Lmax = MaxSizeTL, Kmax = Kmax_lowq),
         biomass_loss_l = somaLoss(M_l, Lmeas = Size, a, b),
         delta_g_day_l = biomass_gain_l - biomass_loss_l,
         biomass_gain_u = somaGain(a, b, Lmeas = Size, Lmax = MaxSizeTL, Kmax = Kmax_uppq),
         biomass_loss_u = somaLoss(M_u, Lmeas = Size, a, b),
         delta_g_day_u = biomass_gain_u - biomass_loss)

summary(vbgf_par)
# don't need to keep traits
vbgf_par <- vbgf_par %>% mutate(rareID = paste(DataID, gridcell, sep = '_')) %>% 
  select(rareID, SamplingEventID, Species, Size, Kmax:delta_g_day_u)

# Apply to assemblage data =====================================================

colnames(dt_rareass)
dt_rareass <- dt_rareass %>% 
  # known uncertain blind species names for productivity calculations
  mutate(prod_genus_sp = case_when(str_detect(genus_sp, '(?<=idae)\\ssp\\d*$') ~ str_remove(genus_sp, '\\ssp\\d*$'),
                                   .default = str_remove(genus_sp, '\\d+')))

# match up
dt_prod <- vbgf_par %>% select(rareID:Size, starts_with('delta')) %>% 
  left_join(dt_rareass, ., by = c('rareID', 'SamplingEventID', 
                                  'prod_genus_sp' = 'Species',
                                  'Length' = 'Size'))
# aggregate to species level
dt_prod <- dt_prod %>% 
  group_by(pick(!c(Abundance, starts_with('delta')))) %>%
  # multiply by count
  summarise(delta_g_day = Abundance * delta_g_day,
            delta_g_day_l = Abundance * delta_g_day_l,
            delta_g_day_u = Abundance * delta_g_day_u) %>% 
  ungroup %>% group_by(pick(!c(Length, delta_g_day, delta_g_day_l, delta_g_day_u))) %>% 
  # aggregate total length and biomass change for the species
  summarise(delta_g_day = sum(delta_g_day, na.rm = T),
            delta_g_day_l = sum(delta_g_day_l, na.rm = T),
            delta_g_day_u = sum(delta_g_day_u, na.rm = T)) %>% ungroup

# aggregate to site level

dt_prod <- dt_prod %>% ungroup %>% 
  group_by(rareID, Year, SamplingEventID) %>%
  summarise(delta_g_day = sum(delta_g_day, na.rm = T),
            delta_g_day_l = sum(delta_g_day_l, na.rm = T),
            delta_g_day_u = sum(delta_g_day_u, na.rm = T)) %>% 
  # convert to kg/yr
  mutate(delta_g_day = delta_g_day * 365 / 1000,
         delta_g_day_l = delta_g_day_l * 365 / 1000,
         delta_g_day_u = delta_g_day_u * 365 / 1000) %>% ungroup

colnames(dt_prod) <- str_replace(colnames(dt_prod), 'delta\\_g\\_day', 'prod_kg_yr')

# kg/ha/yr
dt_prod <- dt_rareass %>% distinct(rareID, DataID) %>% 
  left_join(dt_meta[c(1,7)], ., by = "DataID") %>% 
  mutate(Area = str_extract(`Grain size`, '\\d{2,}')) %>% 
  select(DataID, rareID, Area) %>% 
  left_join(dt_prod, ., by = "rareID")
dt_prod$Area <- as.numeric(dt_prod$Area)
# unit conversion
dt_prod <- dt_prod %>% mutate(across(starts_with('prod'), ~ .x * 10000 / Area))

## Visualise ===================================================================

require(ggplot2)
require(patchwork)

dt_div <- dt_rareass %>% distinct(rareID, DataID) %>% 
  left_join(dt_div, ., by = "rareID")

pretty <- theme_bw(base_size = 13, base_family = "Arial") +
  theme(panel.grid = element_blank())

# diversity
f_asyS <- ggplot(data = dt_div %>% filter(metric == 'S_asy')) +
  # geom_ribbon(aes(x = Year, ymin = lwr, ymax = upr, fill = rareID),
  #             alpha = 0.2) +
  geom_point(aes(x = Year, y = value, color = rareID), 
             shape = 21, size = 2, alpha = 0.6) +
  geom_smooth(aes(x = Year, y = value, color = rareID),
              method = 'loess', se = F) +
  guides(fill = 'none', color = 'none') +
  scale_color_cherulean(palette = 'cheridis', discrete = T) +
  scale_fill_cherulean(palette = 'cheridis', discrete = T) +
  labs(y = 'Extrapolated species richness')

f_S <- ggplot(data = dt_div %>% filter(metric == 'S')) +
  # geom_ribbon(aes(x = Year, ymin = lwr, ymax = upr, fill = rareID),
  #             alpha = 0.2) +
  geom_point(aes(x = Year, y = value, color = rareID), 
             shape = 21, size = 2, alpha = 0.6) +
  geom_smooth(aes(x = Year, y = value, color = rareID),
              method = 'loess', se = F) +
  guides(fill = 'none', color = 'none') +
  scale_color_cherulean(palette = 'cheridis', discrete = T) +
  scale_fill_cherulean(palette = 'cheridis', discrete = T) +
  labs(y = 'Species richness')

f_prod <- ggplot(data = dt_prod) +
  geom_point(aes(x = Year, y = prod_kg_yr, color = rareID),
             shape = 21, size = 2, alpha = 0.6) +
  geom_smooth(aes(x = Year, y = prod_kg_yr, color = rareID), method = 'loess', se = F) +
  guides(fill = 'none', color = 'none') +
  scale_color_cherulean(palette = 'cheridis', discrete = T) +
  labs(y = bquote("Biomass production (kg⋅ha"^-1~"⋅yr"^-1~")"))

(f_S + f_asyS + f_prod) * pretty

# together!

f_proddiv <- left_join(dt_div, dt_prod, by = c('rareID', 'Year', 'SamplingEventID')) %>% 
  filter(metric == 'S') %>% 
  ggplot(data = .) +
  geom_point(aes(x = value, y = prod_kg_yr, color = rareID)) +
  guides(color = 'none') + pretty +
  scale_color_cherulean(palette = 'cheridis', discrete = T) +
  labs(x = 'Species richness', y = bquote("Biomass production (kg⋅ha"^-1~"⋅yr"^-1~")"))

f_proddiv
f_proddiv * scale_y_log10()


left_join(dt_div, dt_prod, by = c('rareID', 'Year', 'SamplingEventID')) %>% 
  filter(metric == 'S_asy') %>% 
  ggplot(data = .) +
  geom_point(aes(x = value, y = prod_kg_yr, color = rareID)) +
  guides(color = 'none') + pretty +
  scale_color_cherulean(palette = 'cheridis', discrete = T) +
  labs(x = 'Estimated species richness', 
       y = bquote("Biomass production (kg⋅ha"^-1~"⋅yr"^-1~")")) +
  scale_x_log10() + scale_y_log10()

ggplot() +
  # geom_ribbon(aes(x = Year, ymin = lwr, ymax = upr, fill = rareID),
  #             alpha = 0.2) +
  geom_point(data = dt_div %>% filter(metric == 'S'), 
             aes(x = Year, y = value * 200, group = rareID), 
             shape = 21, size = 2, alpha = 0.6, color = 'tomato3') +
  geom_smooth(data = dt_div %>% filter(metric == 'S'),
              aes(x = Year, y = value * 200, group = rareID),
              method = 'loess', se = F, size = 0.75, color = 'tomato3') +
  geom_point(data = dt_prod, 
             aes(x = Year, y = prod_kg_yr, group = rareID), 
             shape = 21, size = 2, alpha = 0.6, color = 'mediumorchid3') +
  geom_smooth(data = dt_prod,
              aes(x = Year, y = prod_kg_yr, group = rareID),
              method = 'loess', se = F, size = 0.75, color = 'mediumorchid3') +
  scale_y_continuous(name = "Productivity", 
                     sec.axis = sec_axis(name = "Species richness", trans = ~. / 200)) +
  labs(x = 'Year') + facet_wrap(vars(DataID), scales = 'free_x') + pretty

