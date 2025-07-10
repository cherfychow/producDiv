

# ProducDiv
# Length binning sensitivity test

require(tidytable)
require(rfishprod)
require(ggplot2)
require(patchwork)

dt_length <- dt_rareass %>% ungroup %>%
  # match Kmax
  left_join(., dt_vbgf, 
            by = c('taxaID', 'DataID', 'seriesID', 'SamplingEventID')) %>% 
  filter(is.na(Length_lwr), Length <= MaxSizeTL) %>% sample_n(1000) %>%  
  select(Length, Species, a, b, MaxSizeTL, starts_with('Kmax'))

# get single ones
dt_length <- dt_single %>% ungroup %>%
  # match Kmax
  left_join(., prod_single, 
            by = c('taxaID', 'DataID', 'seriesID', 'SamplingEventID')) %>% 
  filter(is.na(Length_lwr), Length <= MaxSizeTL) %>% sample_n(1000) %>%  
  select(Length, Species, a, b, MaxSizeTL, starts_with('Kmax')) %>% 
  bind_rows(dt_length, .)

dt_length <- dt_length %>% 
  mutate(Length_bin = case_when(Length < 5 ~ 2.5,
                                between(Length, 5, 10) ~ 7.5,
                                Length >= 100 ~ floor(Length/50) * 50 + 25,
                                .default = floor(Length/10) * 10 + 5)) %>%  # artificially bin
  # constrain lengths like we did before
  mutate(Length_bin = case_when(Length_bin >= MaxSizeTL ~ MaxSizeTL - 0.1, .default = Length_bin),
         Length = case_when(Length >= MaxSizeTL ~ MaxSizeTL - 0.1, .default = Length))

# calculate gross productivity
dt_length <- dt_length %>% 
  mutate(bin_prod = somaGain(a, b, Lmeas = Length_bin, Lmax = MaxSizeTL, Kmax = Kmax),
         bin_prod_l = somaGain(a, b, Lmeas = Length_bin, Lmax = MaxSizeTL, Kmax = Kmax_lowq),
         bin_prod_u = somaGain(a, b, Lmeas = Length_bin, Lmax = MaxSizeTL, Kmax = Kmax_uppq),
         act_prod = somaGain(a, b, Lmeas = Length, Lmax = MaxSizeTL, Kmax = Kmax),
         act_prod_l = somaGain(a, b, Lmeas = Length, Lmax = MaxSizeTL, Kmax = Kmax_lowq),
         act_prod_u = somaGain(a, b, Lmeas = Length, Lmax = MaxSizeTL, Kmax = Kmax_uppq))

dt_length$delta_prod <- with(dt_length, abs(bin_prod - act_prod))
dt_length$delta_length <- with(dt_length, abs(Length_bin - Length))
dt_length$delta_length_rel <- with(dt_length, Length/MaxSizeTL)

dt_length <- dt_length %>% mutate(across(contains('prod'), ~.x / 1000)) # to kg/day


f_Lsens_pair <- ggplot(data = dt_length) +
  geom_linerange(aes(xmin = act_prod_l*1000, xmax = act_prod_u*1000, y = bin_prod*1000), 
                 linewidth = 0.25, colour = '#666', alpha = 0.1) +
  geom_linerange(aes(ymin = bin_prod_l*1000, ymax = bin_prod_u*1000, x = act_prod*1000), 
                 linewidth = 0.25, colour = '#666', alpha = 0.1) +
  geom_point(aes(x = act_prod*1000, y = bin_prod*1000), 
             shape = 21, size = 2, fill = '#FFF', alpha = 0.6) +
  labs(x = bquote('Productivity'['actual']~' (g/day)'),
       y = bquote('Productivity'['bin']~' (g/day)')) +
  style + scale_x_log10() + scale_y_log10()

f_Lsens_delta <- ggplot(data = dt_length) +
  geom_hline(aes(yintercept = 1), linetype = 'dashed', colour = '#AAA') +
  geom_point(aes(x = delta_length, y = delta_prod*1000 + 1),
             shape = 21, size = 2, fill = '#FFF', alpha = 0.6) +
  labs(x = '∆Length (cm)',
       y = '∆Productivity (g/day)') +
  style + scale_y_log10()

f_Lsens_rel <- ggplot(data = dt_length) +
  geom_point(aes(x = delta_length_rel, y = delta_prod*1000),
             shape = 21, size = 2, fill = '#FFF', alpha = 0.6) +
  labs(x = bquote('Length/L'['max']),
       y = '∆Productivity (g/day)') +
  style + scale_y_log10()

f_Lsens_max <-  ggplot(data = dt_length) +
  geom_point(aes(x = MaxSizeTL, y = delta_prod*1000),
             shape = 21, size = 2, fill = '#FFF', alpha = 0.6) +
  labs(x = bquote('L'['max']),
       y = '∆Productivity (g/day)') +
  style + scale_y_log10() + scale_x_log10()

f_Lsens_pair + f_Lsens_delta + f_Lsens_rel + f_Lsens_max


# Sample completeness ==============================================================================

require(lme4)
require(dplyr)
load('model_vars.RData')

m_spatial_comp <- glmer(data = vars_single %>% filter(compness > 0.8),
                       prod_kg ~ log(S) + log(hgravity + .0000390505) + (1 + log(S)|biogeo) + (1|DataID) + (1|gridcell:DataID) + (1|grav_ID),
                       family = Gamma(link = 'log'))
save(m_spatial_comp, file = 'm_spatial_comp.RData')

pred_spatial_comp <- vars_single %>% filter(compness > 0.8) %>% 
  mutate(hgravity = mean(hgravity)) %>% select(S, hgravity) %>% 
  distinct
temp <- predict(m_spatial_comp, newdata = pred_spatial_comp, re.form = NA, se.fit = T)
pred_spatial_comp$prod_pred <- exp(temp$fit)
pred_spatial_comp$prod_lwr <- exp(temp$fit - 1.96*temp$se.fit)
pred_spatial_comp$prod_upr <- exp(temp$fit + 1.96*temp$se.fit)

# histogram of sample completeness
(ggplot(data = dt_complete) +
   labs(x = 'Sample completeness', y = 'Count') +
   geom_histogram(aes(x = compness), colour = '#333', fill = 'transparent') + style) +
  
  (f_spatial_dt +
     stat_density_2d(data = vars_single %>% filter(compness < 0.8), aes(x = S, y = prod_kg), colour = '#000', alpha = 0.5) +
     scale_colour_cherulean(palette = 'cheridis', discrete = F, name = 'Completeness') +
     theme(legend.position = 'bottom'))

# comparing model predictions
f_spatial_dt +
  geom_ribbon(data = pred_spatial_comp, aes(x = S, ymin = prod_lwr, ymax = prod_upr), 
              fill = 'blue', alpha = 0.3) +
  geom_line(data = pred_spatial_comp, aes(x = S, y = prod_pred), colour = 'blue') +
  geom_ribbon(data = pred_spatial %>% ungroup %>% group_by(S) %>% 
                summarise(lwr = quantile(prod_S, 0.025),
                          upr = quantile(prod_S, 0.975),
                          prod_S = quantile(prod_S, 0.5)), 
              aes(x = S, ymin = lwr, ymax = upr), 
              fill = 'transparent', colour = '#000', linetype = 'dashed') +
  
  geom_line(data = pred_spatial %>% ungroup %>% group_by(S) %>% 
              summarise(lwr = quantile(prod_S, 0.025),
                        upr = quantile(prod_S, 0.975),
                        prod_S = quantile(prod_S, 0.5)), 
            aes(x = S, y = prod_S), 
            colour = '#000', linewidth = 1)
