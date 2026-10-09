
# Supplementary figures

require(ggplot2)
require(patchwork)
require(ggdist)
require(dplyr)

style <- theme_classic(base_size = 13, base_family = 'Lato') + 
  theme(panel.grid = element_blank(),
    axis.text = element_text(colour = '#222222'),
    axis.ticks = element_line(linewidth = 0.3, colour = '#222222'),
    axis.line = element_line(linewidth = 0.3, colour = '#222222'))

# ==================================================================================================
# Spatial model
# ==================================================================================================

# Posterior distributions

f_space_post1 <- ggplot(dr_spatial %>% filter(var_type == 'main')) +
  stat_slab(aes(x = values, y = prettyvar, fill = after_stat(level)), 
            normalize = 'groups', .width = .95, colour = '#222222', linewidth = 0.25) +  
  scale_fill_manual(values = '#ffffff', na.value = '#ffffff') +
  stat_spike(aes(x = values, y = prettyvar, ), at = 'median',
             normalize = 'groups', size = 0, colour = '#222222') +
  geom_vline(aes(xintercept = 0), linetype = 'dashed', colour = '#000000', alpha = 0.4) +
  labs(x = NULL, y = NULL) + style +
  guides(fill = 'none')

f_space_post2 <- ggplot(dr_spatial %>% filter(var_type == 'randommodel')) +
  stat_slab(aes(x = values, y = prettyvar, fill = after_stat(level)), 
            normalize = 'groups', .width = .95, colour = '#222222', linewidth = 0.25) +  
  scale_fill_manual(values = '#ffffff', na.value = '#ffffff') +
  stat_spike(aes(x = values, y = prettyvar), at = 'median',
             normalize = 'groups', size = 0, colour = '#222222') +
  geom_vline(aes(xintercept = 0), linetype = 'dashed', colour = '#000000', alpha = 0.4) +
  labs(x = NULL, y = NULL) + style +
  guides(fill = 'none')

f_space_post <- f_space_post1 / f_space_post2 + plot_layout(heights = c(2, 5))
f_space_post
ggsave(f_space_post, filename = '../figures/supp/f_spatial_post.svg', width = 10, height = 14, units = 'cm')

# ==================================================================================================
# Spatial model, alternate with N
# ==================================================================================================

# Posterior distributions

f_spaceN_post1 <- ggplot(dr_spatialN %>% filter(var_type == 'main')) +
  stat_slab(aes(x = values, y = prettyvar, fill = after_stat(level)), 
            normalize = 'groups', .width = .95, colour = '#222', linewidth = 0.25) +  
  scale_fill_manual(values = '#fff', na.value = '#fff') +
  stat_spike(aes(x = values, y = prettyvar, ), at = 'median',
             normalize = 'groups', size = 0, colour = '#222') +
  geom_vline(aes(xintercept = 0), linetype = 'dashed', colour = '#000', alpha = 0.4) +
  labs(x = NULL, y = NULL) + style +
  guides(fill = 'none')

f_spaceN_post2 <- ggplot(dr_spatialN %>% filter(var_type == 'randommodel')) +
  stat_slab(aes(x = values, y = prettyvar, fill = after_stat(level)), 
            normalize = 'groups', .width = .95, colour = '#222', linewidth = 0.25) +  
  scale_fill_manual(values = '#fff', na.value = '#fff') +
  stat_spike(aes(x = values, y = prettyvar), at = 'median',
             normalize = 'groups', size = 0, colour = '#222') +
  geom_vline(aes(xintercept = 0), linetype = 'dashed', colour = '#000', alpha = 0.4) +
  labs(x = NULL, y = NULL) + style +
  guides(fill = 'none')

f_spaceN_post <- f_spaceN_post1 / f_spaceN_post2 + plot_layout(heights = c(2, 5))
f_spaceN_post
ggsave(f_spaceN_post, filename = '../figures/supp/f_spatialN_post.svg', width = 10, height = 14, units = 'cm')

# ==================================================================================================
# Biogeography group-level effect posteriors
# ==================================================================================================

dr_geo_b0$biogeo <- factor(dr_geo_b0$biogeo, levels = regions)
dr_geo_S$biogeo <- factor(dr_geo_S$biogeo, levels = regions)

f_geo_b0 <- ggplot(data = dr_geo_b0) +
  stat_slab(aes(x = value, y = biogeo, fill = after_stat(level)), 
            normalize = 'groups', .width = .95, colour = '#222', slab_linewidth = 0.3, alpha = 0.5) +  
  scale_fill_manual(values = '#60dcac', na.value = '#fff') +
  stat_spike(aes(x = value, y = biogeo), 
    p_limits = c(.05, .975), at = c('median', 'qi'),
    normalize = 'groups', linewidth = 0.5, size = 0, colour = '#222') +
  geom_vline(aes(xintercept = -3.52), linetype = 'dashed', colour = '#000', alpha = 0.4) +
  labs(x = NULL, y = NULL) + style +
  guides(fill = 'none')

f_geo_S <- ggplot(data = dr_geo_S) +
  stat_slab(aes(x = value, y = biogeo, fill = after_stat(level)), 
            normalize = 'groups', .width = .95, colour = '#222', slab_linewidth = 0.3, alpha = 0.5) +  
  scale_fill_manual(values = '#60dcac', na.value = '#fff') +
  stat_spike(aes(x = value, y = biogeo), 
    p_limits = c(.05, .975), at = c('median', 'qi'),
    normalize = 'groups', linewidth = 0.5, size = 0, colour = '#222') +
  geom_vline(aes(xintercept = 1.46), linetype = 'dashed', colour = '#000', alpha = 0.4) +
  labs(x = NULL, y = NULL) + style +
  guides(fill = 'none')

f_geo_b0 + f_geo_S

# ==================================================================================================
#  Quantile with N
# ==================================================================================================

ggplot() +
  labs(x = 'N', y = bquote("Productivity (kg⋅ha"^-1~"⋅day"^-1~")")) +
  
  stat_bin_2d(data = vars_single, aes(x = N, y = prod_kg, fill = after_stat(count)), bins = 100) +
  # geom_density_2d(data = vars_single, aes(x = S, y = prod_kg), colour = '#fff', bins = 8, linewidth = 0.3) +
  scale_fill_cherulean(name = 'No. of surveys', palette = 'gomphosus', discrete = F, reverse = T, limits = c(1,200)) +
  
  scale_x_log10() + scale_y_log10(breaks = single_y) + 
  style +

  geom_ribbon(data = pred_quantA %>% filter(quantile != 'med'), 
    aes(x = N, ymin = lower, ymax = higher, group = quantile),
              fill = 'transparent', colour = 'black', linetype = 'dashed') +
  geom_line(data = pred_quantA %>% filter(quantile != 'med'), 
    aes(x = N, y = fit, group = quantile), linewidth = .75) 

# ==================================================================================================
# Null abundance compare
# ==================================================================================================

f_A0 <- ggplot(data = dt_null) +
  labs(x = 'S', y = bquote("Productivity"[A0]~"(kg⋅ha"^-1~"⋅day"^-1~")")) +
  stat_bin_2d(aes(x = S, y = prod0_a, fill = after_stat(count)), bins = 100) +
  scale_fill_cherulean(palette = 'cheridis', discrete = F, reverse = T) +
  style + scale_x_log10() + scale_y_log10()


(f_spatial_dt + f_A0) * theme(legend.position = 'bottom') + plot_annotation(tag_levels = 'A')
ggsave(filename = '../figures/f_abund0.svg', width = 20, height = 12, units = 'cm')


# ==============================================================================
# Temporal model
# ==============================================================================

f_time_post1 <- ggplot(dr_timeN %>% filter(var_type == 'fixed')) +
  stat_slab(aes(x = values, y = prettyvar, fill = after_stat(level)), 
            normalize = 'groups', .width = .95, colour = '#000', slab_linewidth = 0.5) +  
  scale_fill_manual(values = '#fff', na.value = '#fff') +
  stat_spike(aes(x = values, y = prettyvar, ), at = 'median',
             normalize = 'groups', size = 0, colour = '#333660') +
  geom_vline(aes(xintercept = 0), linetype = 'dashed', colour = '#000', alpha = 0.4) +
  labs(x = NULL, y = NULL) + style +
  guides(fill = 'none')

f_time_post2 <- ggplot(dr_timeN %>% filter(var_type == 'randommodel')) +
  stat_slab(aes(x = values, y = prettyvar, fill = after_stat(level)), 
            normalize = 'groups', .width = .95, colour = '#000', slab_linewidth = 0.5) +  
  scale_fill_manual(values = '#fff', na.value = '#fff') +
  stat_spike(aes(x = values, y = prettyvar), at = 'median',
             normalize = 'groups', size = 0, colour = '#333660') +
  geom_vline(aes(xintercept = 0), linetype = 'dashed', colour = '#000', alpha = 0.4) +
  labs(x = NULL, y = NULL) + style +
  guides(fill = 'none')

f_time_post <- f_time_post1 / f_time_post2 + plot_layout(heights = c(2, 5))
f_time_post
ggsave(f_time_post, filename = '../figures/supp/f_time_post.svg', width = 10, height = 14, units = 'cm')


# ∆Prod ~ ∆S, log ratio adjacent years =============================================================

f_flux_post1 <- ggplot() +
  # log ratio fluxes
  stat_slab(data = dr_time_flux %>% filter(var_type == 'fixed'),
            aes(x = values, y = prettyvar, fill = after_stat(level)), 
            normalize = 'groups', .width = .95, linewidth = 0.5, colour = '#222') +  
  scale_fill_manual(values = '#fff', na.value = '#fff') +
  stat_spike(data = dr_time_flux %>% filter(var_type == 'fixed'),
             aes(x = values, y = prettyvar, ), at = 'median',
             normalize = 'groups', size = 0, colour = '#222', linewidth = 0.5) +
  geom_vline(aes(xintercept = 0), linetype = 'dashed', colour = '#000', alpha = 0.4) +
  labs(x = NULL, y = NULL) + style +
  guides(fill = 'none')

f_flux_post2 <- ggplot() +
  stat_slab(data = dr_time_flux %>% filter(var_type == 'randommodel'),
            aes(x = values, y = prettyvar, fill = after_stat(level)), 
            normalize = 'groups', .width = .95, linewidth = 0.5, colour = '#222') +  
  scale_fill_manual(values = '#fff', na.value = '#fff') +
  stat_spike(data = dr_time_flux %>% filter(var_type == 'randommodel'),
             aes(x = values, y = prettyvar), at = 'median',
             normalize = 'groups', size = 0, colour = '#222') +
  geom_vline(aes(xintercept = 0), linetype = 'dashed', colour = '#000', alpha = 0.4) +
  labs(x = NULL, y = NULL) + style +
  guides(fill = 'none')

f_flux_post <- f_flux_post1 / f_flux_post2 + plot_layout(heights = c(1,2))
f_flux_post

ggsave(filename = '../figures/supp/f_time_flux_post.svg', width = 10, height = 14, units = 'cm')


## Partial prediction plot =========================================================================
 

f_flux_1 <- ggplot(data = vars_time_long) +
  labs(x = bquote("∆log"[2]~"S"), y = bquote("∆log"[2]~"Productivity (kg⋅ha"^-1~"⋅day"^-1~")")) +
  geom_hline(aes(yintercept = 0), linetype = 'dashed', colour = '#AAA') +
  geom_vline(aes(xintercept = 0), linetype = 'dashed', colour = '#AAA') +
  geom_linerange(aes(xmin = S_slope - 1.96*S_se, 
                     xmax = S_slope + 1.96*S_se, 
                     y = prod_slope, colour = quad), size = 0.2, alpha = 0.3) +
  geom_linerange(aes(ymin = prod_slope - 1.96*prod_se, 
                     ymax = prod_slope + 1.96*prod_se, 
                     x = S_slope, colour = quad), size = 0.2, alpha = 0.3) +
  scale_colour_manual(values = c('#1c5dac','#1c5dac','#60dcac','#60dcac')) +
  # uncertain points
  geom_point(data = vars_time_long %>% filter(!str_detect(quad, 'certain')),
             aes(x = S_slope, y = prod_slope), 
             size = 2, shape = 21, alpha = 0.3, colour = '#444', fill = '#FFF') +
  # certain slopes
  geom_point(data = vars_time_long %>% filter(str_detect(quad, 'certain')),
             aes(x = S_slope, y = prod_slope, fill = quad), 
             size = 2, shape = 21, colour = '#444444AA') +
  scale_fill_manual(values = c('#1c5dac','#60dcac')) +
  geom_ribbon(data = pred_time %>% ungroup %>% group_by(S_slope) %>% 
                summarise(lwr = quantile(deltaProd_S, 0.025),
                          upr = quantile(deltaProd_S, 0.975),
                          deltaProd_S = quantile(deltaProd_S, 0.5)), 
              aes(x = S_slope, ymin = lwr, ymax = upr), 
              colour = '#000', fill = 'transparent', linetype = 'dashed') +
  geom_line(data = pred_time %>% ungroup %>% group_by(S_slope) %>% 
              summarise(lwr = quantile(deltaProd_S, 0.025),
                        upr = quantile(deltaProd_S, 0.975),
                        deltaProd_S = quantile(deltaProd_S, 0.5)), 
            aes(x = S_slope, y = deltaProd_S), 
            colour = '#000', linewidth = .75) +
  geom_abline(aes(intercept = 0, slope = 1), linetype = 'dashed', colour = '#444', alpha = 0.6) +
  guides(colour = 'none', fill = 'none') +
  coord_cartesian(ylim = c(-3.5,3), xlim = c(-2.5, 1.5)) + style

f_flux_x <- ggplot(data = vars_time_long) +
  geom_density(aes(x = S_slope), fill = '#CCC', colour = '#666') +
  geom_vline(aes(xintercept = 0), linetype = 'dashed', colour = '#666') +
  theme_void() + scale_x_continuous(limits = c(-2.5, 1.5))

f_flux_y <- ggplot(data = vars_time_long) +
  geom_density(aes(y = prod_slope), fill = '#CCC', colour = '#666') +
  geom_hline(aes(yintercept = 0), linetype = 'dashed', colour = '#666') +
  theme_void() + scale_y_continuous(limits = c(-3.5,3))

f_flux <- f_flux_x + plot_spacer() + f_flux_1 + f_flux_y + 
  plot_layout(ncol = 2, nrow = 2, widths = c(5,1), heights = c(1,6))
f_flux

# ==================================================================================================
# Sample completeness plots
# ==================================================================================================

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

# ==================================================================================================
# Length binning sensitivity
# ==================================================================================================


# binned ~ actual
f_L_pair <- ggplot(data = dt_length) +
  geom_linerange(aes(xmin = act_prod_l, xmax = act_prod_u, y = bin_prod), 
                 linewidth = 0.25, colour = '#666', alpha = 0.1) +
  geom_linerange(aes(ymin = bin_prod_l, ymax = bin_prod_u, x = act_prod), 
                 linewidth = 0.25, colour = '#666', alpha = 0.1) +
  geom_point(aes(x = act_prod, y = bin_prod), 
             shape = 21, size = 2, fill = '#FFF', alpha = 0.3) +

  # model prediction
  geom_ribbon(data = pred_L1, aes(x = act_prod, ymin = (fit - 1.96*se) , ymax = (fit + 1.96*se)),
    linetype = 'dashed', colour = '#000', linewidth = 0.3, fill = 'transparent') +
  geom_line(data = pred_L1, aes(x = act_prod, y = fit), linewidth = 0.75) +

  labs(x = bquote('Productivity'['actual']~' (kg/day)'),
       y = bquote('Productivity'['bin']~' (kg/day)')) +
  style + scale_x_log10() + scale_y_log10()

f_L_1 <- ggplot() +
  geom_point(data = dt_length, aes(x = delta_length, y = delta_prod), 
             shape = 21, size = 2, fill = '#FFF', alpha = 0.3) +
  
  # model prediction
  geom_ribbon(data = pred_length, aes(x = delta_length, 
              ymin = fit_deltaL - 1.96 * se_deltaL, ymax = fit_deltaL + 1.96 * se_deltaL),
    linetype = 'dashed', colour = '#000', linewidth = 0.3, fill = 'transparent') +
  geom_line(data = pred_length, aes(x = delta_length, y = fit_deltaL), linewidth = 0.75) +

  labs(x = '∆Length (cm)',
       y = '∆Productivity (kg/day)') + style + scale_y_log10() + scale_x_log10()

f_L_2 <- ggplot() +
  geom_point(data = dt_length, aes(x = MaxSizeTL, y = delta_prod), 
             shape = 21, size = 2, fill = '#FFF', alpha = 0.3) +
  
  # model prediction
  geom_ribbon(data = pred_length, aes(x = MaxSizeTL, 
              ymin = fit_Lmax - 1.96 * se_Lmax, ymax = fit_Lmax + 1.96 * se_Lmax),
    linetype = 'dashed', colour = '#000', linewidth = 0.3, fill = 'transparent') +
  geom_line(data = pred_length, aes(x = MaxSizeTL, y = fit_Lmax), linewidth = 0.75) +

  labs(x = 'Maximum length (cm)',
       y = '∆Productivity (kg/day)') + style + scale_y_log10() + scale_x_log10()

f_L_3 <- ggplot() +
  geom_point(data = dt_length, aes(x = Length/MaxSizeTL, y = delta_prod), 
             shape = 21, size = 2, fill = '#FFF', alpha = 0.3) +
  
  # model prediction
  geom_ribbon(data = pred_length, aes(x = rel_Length, 
              ymin = fit_relL - 1.96 * se_relL, ymax = fit_relL + 1.96 * se_relL),
    linetype = 'dashed', colour = '#000', linewidth = 0.3, fill = 'transparent') +
  geom_line(data = pred_length, aes(x = rel_Length, y = fit_relL), linewidth = 0.75) +

  labs(x = bquote('Length/'~L[max]),
       y = '∆Productivity (kg/day)') + style + scale_y_log10()

f_L_pair
f_L_1 + f_L_2 + f_L_3 + plot_annotation(tag_levels = 'A')

# --------------------------------------------------------------------------------------------------
#  Compare without binned length
# --------------------------------------------------------------------------------------------------

nobins <- dt_single %>% filter(is.na(Length_lwr)) %>% distinct(DataID, SamplingEventID)

ggplot() +
  labs(x = 'S', y = bquote("Productivity (kg⋅ha"^-1~"⋅day"^-1~")")) +
  
  stat_bin_2d(data = vars_single %>% 
    filter(DataID %in% nobins$DataID, SamplingEventID %in% nobins$SamplingEventID), 
    aes(x = S, y = prod_kg, fill = after_stat(count)), bins = 100) +

  geom_point(data = vars_single %>% 
    filter(!DataID %in% nobins$DataID, !SamplingEventID %in% nobins$SamplingEventID),
    aes(x = S, y = prod_kg), shape = 21, size = 2, colour = '#222', alpha = 0.25) +

  scale_fill_cherulean(palette = 'cheridis', discrete = F, reverse = T) +
  
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
            colour = '#000', linewidth = 1) +

  scale_x_log10(breaks = single_x) + scale_y_log10(breaks = single_y) + 
  style

f_spatial_dt
