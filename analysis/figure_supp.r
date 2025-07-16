
# Supplementary figures

require(ggplot2)
require(patchwork)
require(ggdist)
require(dplyr)
require(tidytable)


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


# ==============================================================================
# Posterior distributions
# ==============================================================================

f_time_post1 <- ggplot(dr_time_full %>% filter(var_type == 'fixed')) +
  stat_slab(aes(x = values, y = prettyvar, fill = after_stat(level)), 
            normalize = 'groups', .width = .95, colour = '#000', slab_linewidth = 0.5) +  
  scale_fill_manual(values = '#92d4c7', na.value = '#fff') +
  stat_spike(aes(x = values, y = prettyvar, ), at = 'median',
             normalize = 'groups', size = 0, colour = '#333660') +
  geom_vline(aes(xintercept = 0), linetype = 'dashed', colour = '#000', alpha = 0.4) +
  labs(x = NULL, y = NULL) + style +
  guides(fill = 'none')

f_time_post2 <- ggplot(dr_time_full %>% filter(var_type == 'randommodel')) +
  stat_slab(aes(x = values, y = prettyvar, fill = after_stat(level)), 
            normalize = 'groups', .width = .95, colour = '#000', slab_linewidth = 0.5) +  
  scale_fill_manual(values = '#92d4c7', na.value = '#fff') +
  stat_spike(aes(x = values, y = prettyvar), at = 'median',
             normalize = 'groups', size = 0, colour = '#333660') +
  geom_vline(aes(xintercept = 0), linetype = 'dashed', colour = '#000', alpha = 0.4) +
  labs(x = NULL, y = NULL) + style +
  guides(fill = 'none')

f_time_post <- f_time_post1 / f_time_post2 + plot_layout(heights = c(2, 5))
f_time_post


# NULL MODEL  ======================================================================================

# posterior distributions

f_null_post1 <- ggplot(null_draws %>% filter(var_type == 'fixed')) +
  statime_lab(aes(x = values, y = prettyvar, fill = after_stat(level)), 
            normalize = 'groups', .width = .95) +  
  scale_fill_manual(values = '#60dcac', na.value = 'grey90') +
  statime_pike(aes(x = values, y = prettyvar, ), at = 'median',
             normalize = 'groups', size = 0, colour = '#333660') +
  geom_vline(aes(xintercept = 0), linetype = 'dashed', colour = '#000', alpha = 0.4) +
  labs(x = NULL, y = NULL) + style +
  guides(fill = 'none')

f_null_post2 <- ggplot(null_draws %>% filter(var_type == 'randommodel')) +
  statime_lab(aes(x = values, y = prettyvar, fill = after_stat(level)), 
            normalize = 'groups', .width = .95) +  
  scale_fill_manual(values = '#60dcac', na.value = 'grey90') +
  statime_pike(aes(x = values, y = prettyvar), at = 'median',
             normalize = 'groups', size = 0, colour = '#333660') +
  geom_vline(aes(xintercept = 0), linetype = 'dashed', colour = '#000', alpha = 0.4) +
  labs(x = NULL, y = NULL) + style +
  guides(fill = 'none')

f_null_post <- f_null_post1 / f_null_post2 + plot_layout(heights = c(2, 5))
f_null_post

f_null <- ggplot(data = vars_single) +
  labs(x = 'S', y = bquote("Productivity (kg⋅ha"^-1~"⋅day"^-1~")")) +
  geom_point(aes(x = S, y = prod_null), size = 1, shape = 21, 
             colour = '#666', alpha = 0.2) +
  
  # geom_point(aes(x = S, y = prod_kg, colour = biogeo, fill = biogeo), 
  #            shape = 21, size = 2, alpha = 0.5) +
  # scale_colour_cherulean(palette = 'cheridis', discrete = T) +
  # scale_fill_cherulean(palette = 'cheridis', discrete = T) +
  # 
  # geom_ribbon(data = spatial_partialpred %>% ungroup %>% group_by(S) %>% 
  #               summarise(lwr = quantile(prod_S, 0.025),
  #                         upr = quantile(prod_S, 0.975),
  #                         prod_S = quantile(prod_S, 0.5)), 
  #             aes(x = S, ymin = lwr, ymax = upr), 
  #             fill = 'transparent', colour = '#000', linetype = 'dashed') +
  # 
  # geom_line(data = spatial_partialpred %>% ungroup %>% group_by(S) %>% 
  #             summarise(lwr = quantile(prod_S, 0.025),
  #                       upr = quantile(prod_S, 0.975),
  #                       prod_S = quantile(prod_S, 0.5)), 
  #           aes(x = S, y = prod_S), 
  #           colour = '#000', linewidth = 0.75) +
  
  geom_ribbon(data = null_partialpred %>% ungroup %>% group_by(S) %>%
                summarise(lwr = quantile(prod_S, 0.025),
                          upr = quantile(prod_S, 0.975),
                          prod_S = quantile(prod_S, 0.5)),
              aes(x = S, ymin = lwr, ymax = upr),
              fill = 'transparent', colour = '#000', linetype = 'dashed') +

  geom_line(data = null_partialpred %>% ungroup %>% group_by(S) %>%
              summarise(lwr = quantile(prod_S, 0.025),
                        upr = quantile(prod_S, 0.975),
                        prod_S = quantile(prod_S, 0.5)),
            aes(x = S, y = prod_S),
            colour = '#000', linewidth = 0.75) +
  scale_x_log10() + scale_y_log10() + style +
  guides(colour = 'none', fill = 'none')
f_null

# Spatial single time point model ==================================================================

## Posterior distributions

f_space_post1 <- ggplot(spatial_draws %>% filter(var_type == 'fixed')) +
  statime_lab(aes(x = values, y = prettyvar, fill = after_stat(level)), 
            normalize = 'groups', .width = .95) +  
  scale_fill_manual(values = '#60dcac', na.value = 'grey90') +
  statime_pike(aes(x = values, y = prettyvar, ), at = 'median',
             normalize = 'groups', size = 0, colour = '#333660') +
  geom_vline(aes(xintercept = 0), linetype = 'dashed', colour = '#000', alpha = 0.4) +
  labs(x = NULL, y = NULL) + style +
  guides(fill = 'none')

f_space_post2 <- ggplot(spatial_draws %>% filter(var_type == 'randommodel')) +
  statime_lab(aes(x = values, y = prettyvar, fill = after_stat(level)), 
            normalize = 'groups', .width = .95) +  
  scale_fill_manual(values = '#60dcac', na.value = 'grey90') +
  statime_pike(aes(x = values, y = prettyvar), at = 'median',
             normalize = 'groups', size = 0, colour = '#333660') +
  geom_vline(aes(xintercept = 0), linetype = 'dashed', colour = '#000', alpha = 0.4) +
  labs(x = NULL, y = NULL) + style +
  guides(fill = 'none')

f_space_post <- f_space_post1 / f_space_post2 + plot_layout(heights = c(2, 5))
f_space_post

## Partial predictions 




# looking at null models

singlept <- ggplot(data = vars_single) +
  labs(x = 'S', y = bquote("Productivity (kg⋅ha"^-1~"⋅day"^-1~")")) +
  
  geom_point(aes(x = S, y = prod_kg, colour = biogeo),
             shape = 21, size = 1, alpha = 0.1) +
  scale_colour_cherulean(palette = 'cheridis', discrete = T) +
  scale_x_log10() + scale_y_log10() + style

# compare null abundance
(ggplot(data = vars_single) +
    labs(x = 'S', y = bquote("Productivity (kg⋅ha"^-1~"⋅day"^-1~")"), title = 'Original data') +
    
    geom_point(aes(x = S, y = prod_kg, colour = log(totalN)),
               shape = 21, size = 1, alpha = 0.1) + style + 
    scale_x_log10() + scale_y_log10() + 
    scale_colour_continuous(type = 'viridis') +
ggplot() + 
    labs(x = 'S', y = bquote("Productivity (kg⋅ha"^-1~"⋅day"^-1~")"), title = 'Null data') +
    geom_point(data = dt_null, aes(x = S, y = prod0_Amin, colour = log(totalN)), shape = 21, alpha = 0.1) + 
    style + scale_x_log10() + scale_y_log10() +
    scale_colour_continuous(type = 'viridis')) + plot_layout(guides = 'collect')

(ggplot(data = vars_single) +
    labs(x = 'Total abundance', y = bquote("Productivity (kg⋅ha"^-1~"⋅day"^-1~")"), title = 'Original data') +
    
    geom_point(aes(x = totalN, y = prod_kg, colour = log(totalN)),
               shape = 21, size = 1, alpha = 0.1) + style + 
    scale_x_log10() + scale_y_log10() + 
    scale_colour_continuous(type = 'viridis') +
    ggplot() + 
    labs(x = 'Total abundance', y = bquote("Productivity (kg⋅ha"^-1~"⋅day"^-1~")"), title = 'Null data') +
    geom_point(data = dt_null, aes(x = totalN, y = prod0_Amin, colour = log(totalN)), shape = 21, alpha = 0.1) + 
    style + scale_x_log10() + scale_y_log10() +
    scale_colour_continuous(type = 'viridis')) + plot_layout(guides = 'collect')


f_spaceNull <- singlept +
  
  # model estimate
  geom_ribbon(data = spatial_partialpred %>% ungroup %>% group_by(S) %>% 
                summarise(lwr = quantile(prod_S, 0.025),
                          upr = quantile(prod_S, 0.975),
                          prod_S = quantile(prod_S, 0.5)), 
              aes(x = S, ymin = lwr, ymax = upr), 
              fill = 'transparent', colour = '#000', linetype = 'dashed') +
  
  geom_line(data = spatial_partialpred %>% ungroup %>% group_by(S) %>% 
              summarise(lwr = quantile(prod_S, 0.025),
                        upr = quantile(prod_S, 0.975),
                        prod_S = quantile(prod_S, 0.5)), 
            aes(x = S, y = prod_S), 
            colour = '#000', linewidth = 1) +
  
  # null model lines
  # max Kmax
  geom_ribbon(data = nullsp1_partial %>% ungroup %>% group_by(S) %>% 
                summarise(lwr = quantile(prod_S, 0.025),
                          upr = quantile(prod_S, 0.975),
                          prod_S = quantile(prod_S, 0.5)), 
              aes(x = S, ymin = lwr, ymax = upr), 
              fill = '#000', alpha = 0.1) +
  
  geom_line(data = nullsp1_partial %>% ungroup %>% group_by(S) %>% 
              summarise(lwr = quantile(prod_S, 0.025),
                        upr = quantile(prod_S, 0.975),
                        prod_S = quantile(prod_S, 0.5)), 
            aes(x = S, y = prod_S), 
            colour = '#000', alpha = 0.5, linewidth = 0.75) +
  
  # min Kmax
  geom_ribbon(data = nullsp2_partial %>% ungroup %>% group_by(S) %>% 
                summarise(lwr = quantile(prod_S, 0.025),
                          upr = quantile(prod_S, 0.975),
                          prod_S = quantile(prod_S, 0.5)), 
              aes(x = S, ymin = lwr, ymax = upr), 
              fill = '#000', alpha = 0.1) +
  
  geom_line(data = nullsp2_partial %>% ungroup %>% group_by(S) %>% 
              summarise(lwr = quantile(prod_S, 0.025),
                        upr = quantile(prod_S, 0.975),
                        prod_S = quantile(prod_S, 0.5)), 
            aes(x = S, y = prod_S), 
            colour = '#000', alpha = 0.8, linewidth = 0.75) +
  
  # null productivity
  geom_ribbon(data = null_partialpred %>% ungroup %>% group_by(S) %>%
                summarise(lwr = quantile(prod_S, 0.025),
                          upr = quantile(prod_S, 0.975),
                          prod_S = quantile(prod_S, 0.5)),
              aes(x = S, ymin = lwr, ymax = upr),
              fill = '#444444', alpha = 0.2) +
  
  geom_line(data = null_partialpred %>% ungroup %>% group_by(S) %>%
              summarise(lwr = quantile(prod_S, 0.025),
                        upr = quantile(prod_S, 0.975),
                        prod_S = quantile(prod_S, 0.5)),
            aes(x = S, y = prod_S),
            colour = 'grey30', alpha = 0.5, linewidth = 0.75, linetype = 'dashed') +
  
  geom_ribbon(data = nullP2_partialpred %>% ungroup %>% group_by(S) %>%
                summarise(lwr = quantile(prod0_nullP2, 0.025),
                          upr = quantile(prod0_nullP2, 0.975),
                          prod_S = quantile(prod0_nullP2, 0.5)),
              aes(x = S, ymin = lwr, ymax = upr),
              fill = 'blue', alpha = 0.2) +
  
  geom_line(data = nullP2_partialpred %>% ungroup %>% group_by(S) %>%
              summarise(lwr = quantile(prod0_nullP2, 0.025),
                        upr = quantile(prod0_nullP2, 0.975),
                        prod_S = quantile(prod0_nullP2, 0.5)),
            aes(x = S, y = prod_S),
            colour = 'blue', alpha = 0.5, linewidth = 0.75, linetype = 'dashed') +
  
  # null abundance
  geom_ribbon(data = abund0_partial %>% ungroup %>% group_by(S) %>% 
                summarise(lwr = quantile(prod_S, 0.025),
                          upr = quantile(prod_S, 0.975),
                          prod_S = quantile(prod_S, 0.5)), 
              aes(x = S, ymin = lwr, ymax = upr), 
              fill = 'darkred', alpha = 0.1) +
  
  geom_line(data = abund0_partial %>% ungroup %>% group_by(S) %>% 
              summarise(lwr = quantile(prod_S, 0.025),
                        upr = quantile(prod_S, 0.975),
                        prod_S = quantile(prod_S, 0.5)), 
            aes(x = S, y = prod_S), 
            colour = 'darkred', alpha = 0.5, linewidth = 0.75) +
  
  scale_x_continuous(breaks = c(1,5,10,50,100,300,500)) + 
  scale_y_continuous(breaks = c(0,0.001,.01, .1, 1, 10, 50, 100)) +
  coord_trans(y = "log10", x = 'log10') + style +
  guides(colour = 'none', fill = 'none')

# with human gravity
f_space_g <- ggplot(data = vars_single) +
  labs(x = 'Human gravity', y =  bquote("Productivity (kg⋅ha"^-1~"⋅day"^-1~")")) +
  geom_errorbar(aes(x = hgravity + 1, ymin = prod_kg_l, ymax = prod_kg_u),
                linewidth = 0.2, colour = '#666') +
  geom_point(aes(x = hgravity + 1, y = prod_kg, colour = biogeo, fill = biogeo), 
             shape = 21, size = 2, alpha = 0.5) +
  scale_colour_cherulean(palette = 'cheridis', discrete = T) +
  scale_fill_cherulean(palette = 'cheridis', discrete = T) +
  
  geom_ribbon(data = spatial_partialpred %>% ungroup %>% group_by(hgravity) %>% 
                summarise(lwr = quantile(prod_hgrav, 0.025),
                          upr = quantile(prod_hgrav, 0.975),
                          prod_hgrav = quantile(prod_hgrav, 0.5)), 
              aes(x = hgravity + 1, ymin = lwr, ymax = upr), 
              fill = 'transparent', colour = '#000', linetype = 'dashed') +
  geom_line(data = spatial_partialpred %>% ungroup %>% group_by(hgravity) %>% 
              summarise(lwr = quantile(prod_hgrav, 0.025),
                        upr = quantile(prod_hgrav, 0.975),
                        prod_hgrav = quantile(prod_hgrav, 0.5)), 
            aes(x = hgravity + 1, y = prod_hgrav), 
            colour = '#333660', linewidth = 0.75) +
  
  scale_y_continuous(breaks = c(0,0.001,.01, .1, 1, 10, 50, 100)) +
  coord_trans(y = "log10", x = 'log10') + style +
  guides(colour = 'none', fill = 'none')



# Biogeographic single time point model ============================================================

## Posterior distributions 

f_geo_post1 <- ggplot(dr_spatial %>% filter(var_type == 'main')) +
  statime_lab(aes(x = values, y = prettyvar, fill = after_stat(level)), 
            normalize = 'groups', .width = .95, colour = '#666', linewidth = 0.2) +  
  scale_fill_manual(values = '#BBB', na.value = '#FFF') +
  statime_pike(aes(x = values, y = prettyvar, ), at = 'median',
             normalize = 'groups', size = 0, colour = '#333660') +
  geom_vline(aes(xintercept = 0), linetype = 'dashed', colour = '#000', alpha = 0.4) +
  labs(x = NULL, y = NULL) + style +
  guides(fill = 'none')

f_geo_post2 <- ggplot(dr_spatial %>% filter(var_type == 'randommodel')) +
  statime_lab(aes(x = values, y = prettyvar, fill = after_stat(level)), 
            normalize = 'groups', .width = .95, colour = '#666', linewidth = 0.2) +  
  scale_fill_manual(values = '#BBB', na.value = '#FFF') +
  statime_pike(aes(x = values, y = prettyvar), at = 'median',
             normalize = 'groups', size = 0, colour = '#333660') +
  geom_vline(aes(xintercept = 0), linetype = 'dashed', colour = '#000', alpha = 0.4) +
  labs(x = NULL, y = NULL) + style +
  guides(fill = 'none')

f_geo_post <- f_geo_post1 / f_geo_post2 + plot_layout(heights = c(2, 5))
f_geo_post

## Partial predictions

# split by bioregion


ggplot(data = vars_single) +
  labs(x = 'log(S)', y = bquote("log Productivity (kg⋅ha"^-1~"⋅day"^-1~")")) +
  
  geom_line(data = geo_partialpred %>% ungroup %>% group_by(biogeo, S) %>% 
              summarise(lwr = quantile(prod_S, 0.025),
                        upr = quantile(prod_S, 0.975),
                        prod_S = quantile(prod_S, 0.5)), 
            aes(x = S, y = prod_S, colour = biogeo), linewidth = 1.2, alpha = 0.3) +
  
  scale_x_log10() + scale_y_log10() + style



## Residuals =======================================================================================

geo_resid1 <- ggplot(data = spatial_resid) + style +
  labs(x = 'Data', y = 'Model predicted') +
  geom_point(aes(x = prod_kg, y = prod_S), 
             shape = 21, colour = '#000', alpha = 0.2) +
  scale_x_log10() + scale_y_log10()

geo_resid2 <- ggplot(data = spatial_resid) + style +
  labs(x = 'S', y = 'Model residuals') +
  geom_point(aes(x = S, y = resid), 
             shape = 21, colour = '#000', alpha = 0.2) +
  scale_x_log10()

geo_resid3 <- ggplot(data = spatial_resid) + style +
  labs(x = 'S', y = 'Model residuals') +
  geom_point(aes(x = S, y = abs(resid), 
             shape = 21, colour = '#000', alpha = 0.2) +
  scale_x_log10() + scale_y_log10()

geo_resid1 + geo_resid2 + geo_resid3




# ∆Prod ~ ∆S, log ratio adjacent years =============================================================

f_flux_post1 <- ggplot() +
  # log ratio fluxes
  statime_lab(data = dr_time_full %>% filter(var_type == 'fixed'),
            aes(x = values, y = prettyvar, fill = after_stat(level)), 
            normalize = 'groups', .width = .95) +  
  scale_fill_manual(values = 'darkgreen', na.value = 'grey90') +
  statime_pike(data = dr_time_full %>% filter(var_type == 'fixed'),
             aes(x = values, y = prettyvar, ), at = 'median',
             normalize = 'groups', size = 0, colour = 'darkgreen') +
  # full time series
  statime_lab(data = dr_time_full %>% filter(var_type == 'fixed'),
            aes(x = values, y = prettyvar, fill = after_stat(level)), 
            normalize = 'groups', .width = .95) +  
  statime_pike(data = dr_time_full %>% filter(var_type == 'fixed'),
             aes(x = values, y = prettyvar, ), at = 'median',
             normalize = 'groups', size = 0, colour = '#000') +
  geom_vline(aes(xintercept = 0), linetype = 'dashed', colour = '#000', alpha = 0.4) +
  labs(x = NULL, y = NULL) + style +
  guides(fill = 'none')

f_flux_post2 <- ggplot() +
  statime_lab(data = dr_time_full %>% filter(var_type == 'randommodel'),
            aes(x = values, y = prettyvar, fill = after_stat(level)), 
            normalize = 'groups', .width = .95) +  
  scale_fill_manual(values = 'darkgreen', na.value = 'grey90') +
  statime_pike(data = dr_time_full %>% filter(var_type == 'randommodel'),
             aes(x = values, y = prettyvar), at = 'median',
             normalize = 'groups', size = 0, colour = 'darkgreen') +
  
  # full time series
  statime_lab(data = dr_time_full %>% filter(var_type == 'randommodel'),
            aes(x = values, y = prettyvar, fill = after_stat(level)), 
            normalize = 'groups', .width = .95) +  
  statime_pike(data = dr_time_full %>% filter(var_type == 'randommodel'),
             aes(x = values, y = prettyvar, ), at = 'median',
             normalize = 'groups', size = 0, colour = '#000') +
  geom_vline(aes(xintercept = 0), linetype = 'dashed', colour = '#000', alpha = 0.4) +
  labs(x = NULL, y = NULL) + style +
  guides(fill = 'none')

f_flux_post <- f_flux_post1 / f_flux_post2 + plot_layout(heights = c(1,2))
f_flux_post


## Partial prediction plot =========================================================================
 

f_time_1 <- ggplot(data = vars_time_long) +
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

f_time_x <- ggplot(data = vars_time_long) +
  geom_density(aes(x = S_slope), fill = '#CCC', colour = '#666') +
  geom_vline(aes(xintercept = 0), linetype = 'dashed', colour = '#666') +
  theme_void() + scale_x_continuous(limits = c(-2.5, 1.5))

f_time_y <- ggplot(data = vars_time_long) +
  geom_density(aes(y = prod_slope), fill = '#CCC', colour = '#666') +
  geom_hline(aes(yintercept = 0), linetype = 'dashed', colour = '#666') +
  theme_void() + scale_y_continuous(limits = c(-3.5,3))

f_time <- f_time_x + plot_spacer() + f_time_1 + f_time_y + 
  plot_layout(ncol = 2, nrow = 2, widths = c(5,1), heights = c(1,6))
f_time

# NULL SPECIES MODEL ===============================================================================

## Posterior distributions =========================================================================

f_nullsp1_post1 <- ggplot(nullsp1_draws %>% filter(var_type == 'fixed')) +
  statime_lab(aes(x = values, y = prettyvar, fill = after_stat(level)), 
            normalize = 'groups', .width = .95) +  
  scale_fill_manual(values = '#60dcac', na.value = 'grey90') +
  statime_pike(aes(x = values, y = prettyvar, ), at = 'median',
             normalize = 'groups', size = 0, colour = '#333660') +
  geom_vline(aes(xintercept = 0), linetype = 'dashed', colour = '#000', alpha = 0.4) +
  labs(x = NULL, y = NULL) + style +
  guides(fill = 'none')

f_nullsp1_post2 <- ggplot(nullsp1_draws %>% filter(var_type == 'randommodel')) +
  statime_lab(aes(x = values, y = prettyvar, fill = after_stat(level)), 
            normalize = 'groups', .width = .95) +  
  scale_fill_manual(values = '#60dcac', na.value = 'grey90') +
  statime_pike(aes(x = values, y = prettyvar), at = 'median',
             normalize = 'groups', size = 0, colour = '#333660') +
  geom_vline(aes(xintercept = 0), linetype = 'dashed', colour = '#000', alpha = 0.4) +
  labs(x = NULL, y = NULL) + style +
  guides(fill = 'none')

f_nullsp1_post <- f_nullsp1_post1 / f_nullsp1_post2 + plot_layout(heights = c(2, 5))
f_nullsp1_post

## Partial predictions =============================================================================

f_nullsp <- ggplot(data = vars_single) +
  labs(x = 'S', y = bquote("Productivity (kg⋅ha"^-1~"⋅day"^-1~")")) +
  
  geom_point(aes(x = S, y = prod_maxsp), 
             shape = 21, alpha = 0.2, colour = 'grey60') +
  
  geom_point(aes(x = S, y = prod_minsp), 
             shape = 21, alpha = 0.2, colour = 'grey80') +
  
  geom_ribbon(data = nullsp1_partial %>% ungroup %>% group_by(S) %>% 
                summarise(lwr = quantile(prod_S, 0.025),
                          upr = quantile(prod_S, 0.975),
                          prod_S = quantile(prod_S, 0.5)), 
              aes(x = S, ymin = lwr, ymax = upr), 
              fill = 'transparent', colour = 'darkred', linetype = 'dashed') +
  
  geom_line(data = nullsp1_partial %>% ungroup %>% group_by(S) %>% 
              summarise(lwr = quantile(prod_S, 0.025),
                        upr = quantile(prod_S, 0.975),
                        prod_S = quantile(prod_S, 0.5)), 
            aes(x = S, y = prod_S), 
            colour = 'darkred', linewidth = 0.75) +
  
  geom_ribbon(data = nullsp2_partial %>% ungroup %>% group_by(S) %>% 
                summarise(lwr = quantile(prod_S, 0.025),
                          upr = quantile(prod_S, 0.975),
                          prod_S = quantile(prod_S, 0.5)), 
              aes(x = S, ymin = lwr, ymax = upr), 
              fill = 'transparent', colour = 'darkred', linetype = 'dashed') +
  
  geom_line(data = nullsp2_partial %>% ungroup %>% group_by(S) %>% 
              summarise(lwr = quantile(prod_S, 0.025),
                        upr = quantile(prod_S, 0.975),
                        prod_S = quantile(prod_S, 0.5)), 
            aes(x = S, y = prod_S), 
            colour = 'darkred', linewidth = 0.75) +
  
  # geom_ribbon(data = geo_partialpred %>% ungroup %>% group_by(biogeo, S) %>% 
  #               summarise(lwr = quantile(prod_S, 0.025),
  #                         upr = quantile(prod_S, 0.975),
  #                         prod_S = quantile(prod_S, 0.5)), 
  #             aes(x = S, ymin = lwr, ymax = upr, group = biogeo), 
  #             fill = 'transparent', colour = '#000', linetype = 'dashed') +
  
  geom_line(data = spatial_partialpred %>% ungroup %>% group_by(S) %>% 
              summarise(lwr = quantile(prod_S, 0.025),
                        upr = quantile(prod_S, 0.975),
                        prod_S = quantile(prod_S, 0.5)), 
            aes(x = S, y = prod_S), 
            colour = '#000', linewidth = 0.75) +
  
  scale_x_log10() + scale_y_log10() + style +
  guides(colour = 'none', fill = 'none')