

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
  geom_linerange(aes(xmin = act_prod_l, xmax = act_prod_u, y = bin_prod), 
                 linewidth = 0.25, colour = '#666', alpha = 0.1) +
  geom_linerange(aes(ymin = bin_prod_l, ymax = bin_prod_u, x = act_prod), 
                 linewidth = 0.25, colour = '#666', alpha = 0.1) +
  geom_point(aes(x = act_prod, y = bin_prod), 
             shape = 21, size = 2, fill = '#FFF', alpha = 0.6) +
  labs(x = bquote('Productivity'['actual']~' (kg/day)'),
       y = bquote('Productivity'['bin']~' (kg/day)')) +
  style + scale_x_log10() + scale_y_log10()

f_Lsens_delta <- ggplot(data = dt_length) +
  geom_hline(aes(yintercept = 1), linetype = 'dashed', colour = '#AAA') +
  geom_point(aes(x = delta_length, y = delta_prod + 1),
             shape = 21, size = 2, fill = '#FFF', alpha = 0.6) +
  labs(x = '∆Length (cm)',
       y = '∆Productivity + 1 (kg/day)') +
  style + coord_trans(y = 'log10')

f_Lsens_rel <- ggplot(data = dt_length) +
  geom_hline(aes(yintercept = 1), linetype = 'dashed', colour = '#AAA') +
  geom_point(aes(x = delta_length_rel, y = delta_prod + 1),
             shape = 21, size = 2, fill = '#FFF', alpha = 0.6) +
  labs(x = bquote('Length/L'['max']),
       y = '∆Productivity + 1 (kg/day)') +
  style + coord_trans(y = 'log10')

f_Lsens_max <-  ggplot(data = dt_length) +
  geom_hline(aes(yintercept = 1), linetype = 'dashed', colour = '#AAA') +
  geom_point(aes(x = MaxSizeTL, y = delta_prod + 1),
             shape = 21, size = 2, fill = '#FFF', alpha = 0.6) +
  labs(x = bquote('L'['max']),
       y = '∆Productivity + 1 (kg/day)') +
  style + coord_trans(x = 'log10', y = 'log10')

f_Lsens_pair + f_Lsens_delta + f_Lsens_rel + f_Lsens_max
