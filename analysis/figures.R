
# ProducDiv
# Figure code

require(ggplot2)
require(patchwork)
require(sf)
require(rnaturalearth)
require(ggdist)
require(ggpattern)

# colour palettes
source('https://gist.githubusercontent.com/cherfychow/e9ae890fd16f4c86730748c067feee2b/raw/b2db138ab8164c237129308ea643de78cd54b252/cherulean.R')
style <- theme_classic(base_size = 15, base_family = 'Lato') + 
  theme(panel.grid = element_blank(),
    axis.text = element_text(colour = '#222'),
    axis.ticks = element_line(linewidth = 0.3, colour = '#222'),
    axis.line = element_line(linewidth = 0.3, colour = '#222'))

single_x <- c(1,5,10,50,100,150)
single_y <- c(0,0.001,.01, .1, 1, 10, 100, 500)


# ==================================================================================================
# Figure 1: Dataset map
# ==================================================================================================

proj <- '+proj=eqearth +lon_0=160 +x_0=0 +y_0=0 +R=6371008.7714 +units=m +no_defs +type=crs'

world <- ne_coastline('medium', returnclass = 'sf') %>%
  st_make_valid(.) %>% st_break_antimeridian(., lon_0 = 160)

# make a geo point object with the data
point_rare <- dt_rareass %>% ungroup %>% group_by(DataID, seriesID) %>%
  summarise(Latitude = mean(Latitude), Longitude = mean(Longitude), Duration = diff(range(Year)) + 1) %>% 
  ungroup %>% distinct() %>% 
  left_join(vars_time2, ., by = c('DataID', 'seriesID')) %>% 
  st_as_sf(., coords = c('Longitude', 'Latitude'), crs = 4326)

point_single <- vars_single %>% distinct(DataID, seriesID, SamplingEventID, S) %>% 
  left_join(., distinct(dt_single, DataID, seriesID, SamplingEventID, Latitude, Longitude), 
            by = c('DataID', 'seriesID', 'SamplingEventID')) %>% 
  group_by(DataID, seriesID, Latitude, Longitude) %>% 
  summarise(meanS = mean(S)) %>% 
  filter(is.na(Latitude) == F) %>%
  st_as_sf(., crs = 4326, coords = c("Longitude", "Latitude")) %>%
  st_make_valid %>% st_cast(., to = "MULTIPOINT")

# by dataset ID
f_map_rare <- ggplot() +
  geom_sf(data = world, color='grey75', fill='grey95') +
  geom_sf(data = point_rare, 
          aes(fill = meanS, size = Duration), shape = 21, alpha = 0.5) +
  coord_sf(crs = proj, ylim = c(-4366274.76, 4366274.76)) +
  theme_minimal() + guides(fill = 'none') + theme(legend.position = 'bottom') +
  # scale_fill_cherulean(palette = 'cheridis', discrete = F, name = 'mean(S)', 
  #                      limits = c(1,150), breaks = c(1,10,50,100,150)) +
  scale_fill_gradientn(name = 'mean(S)', colours = cherulean_cols(palette = 'gomphosus')(5), 
                       values = c(1,10,50,100,150)/150, breaks = c(1,10,50,100,150)) +
  scale_size_continuous(breaks = c(2,5,10,15,20,30), range = c(2,6))
f_map_single <- ggplot() +
  geom_sf(data = world, color='grey75', fill='grey95') +
  geom_sf(data = point_single, 
          aes(fill = meanS), shape = 21, size = 2, alpha = 0.5) +
  coord_sf(crs = proj, ylim = c(-4366274.76, 4366274.76)) +
  theme_minimal() + theme(legend.position = 'bottom') +
  scale_fill_gradientn(name = 'mean(S)', colours = cherulean_cols(palette = 'gomphosus')(5), 
                       values = c(1,10,50,100,150)/150, breaks = c(1,10,50,100,150))
  # scale_fill_cherulean(palette = 'cheridis', discrete = F,name = 'mean(S)',
                       # limits = c(1,150), breaks = c(1,10,50,100,150))

dt_biogeo <- dt_assemblage %>% distinct(DataID, SamplingEventID, gridcell) %>% 
  left_join(., grid_sf, by = 'gridcell') %>% group_by(biogeo) %>% 
  summarise(studyN = n_distinct(DataID), eventN = n_distinct(SamplingEventID))
dt_biogeo$biogeo <- factor(dt_biogeo$biogeo, 
                           levels = regions)

f_regionN <- ggplot(data = dt_biogeo) +
  geom_col(aes(y = biogeo, x = eventN, fill = biogeo)) +
  scale_fill_cherulean(palette = 'cheridis', discrete = T) + 
  guides(fill = 'none') + style + scale_x_log10()

(f_map_rare / f_map_single)

f_map <- ggplot() +
  geom_sf(data = world, color='grey75', fill='grey95') +
  geom_sf(data = points, 
          aes(fill = DataID), shape = 21, size = 2.5) +
  coord_sf(crs = proj, ylim = c(-4366274.76, 4366274.76)) +
  theme_minimal(base_size = 13, base_family = 'Helvetica') +
  guides(fill = 'none') +
  scale_fill_cherulean(palette = 'cheridis', discrete = T, alpha = 0.4)



# ==================================================================================================
# Figure 2: Spatial model with quantile + nulls
# ==================================================================================================

# Base single year data

f_spatial_dt <- ggplot() +
  labs(x = 'S', y = bquote("Productivity (kg⋅ha"^-1~"⋅day"^-1~")")) +
  
  stat_bin_2d(data = vars_single, aes(x = S, y = prod_kg, fill = after_stat(count)), bins = 100) +
  # geom_density_2d(data = vars_single, aes(x = S, y = prod_kg), colour = '#fff', bins = 8, linewidth = 0.3) +
  scale_fill_cherulean(palette = 'gomphosus', discrete = F, reverse = T) +
  
  scale_x_log10(breaks = single_x) + scale_y_log10(breaks = single_y) + 
  style

# heat map + spatial partial prediction
f_spatial <- f_spatial_dt + 
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

# quantile regression
f_spatial2 <- f_spatial_dt + 
  geom_ribbon(data = pred_quant %>% filter(quantile != 'med'), 
    aes(x = S, ymin = lower, ymax = higher, group = quantile),
              fill = 'transparent', colour = '#000', linetype = 'dashed') +
  geom_line(data = pred_quant %>% filter(quantile != 'med'), 
    aes(x = S, y = fit, group = quantile), linewidth = .75) 

# null models
f_spatial3 <- f_spatial_dt + 
  # model prediction
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
  
  # null model lines
  # max Kmax
  geom_ribbon(data = pred_sp0_max %>% ungroup %>% group_by(S) %>% 
                summarise(lwr = quantile(prod_S, 0.025),
                          upr = quantile(prod_S, 0.975),
                          prod_S = quantile(prod_S, 0.5)), 
              aes(x = S, ymin = lwr, ymax = upr), 
              fill = '#000', alpha = 0.1) +
  
  geom_line(data = pred_sp0_max %>% ungroup %>% group_by(S) %>% 
              summarise(lwr = quantile(prod_S, 0.025),
                        upr = quantile(prod_S, 0.975),
                        prod_S = quantile(prod_S, 0.5)), 
            aes(x = S, y = prod_S), 
            colour = '#000', alpha = 0.5, linewidth = 0.75) +
  
  # min Kmax
  geom_ribbon(data = pred_sp0_min %>% ungroup %>% group_by(S) %>% 
                summarise(lwr = quantile(prod_S, 0.025),
                          upr = quantile(prod_S, 0.975),
                          prod_S = quantile(prod_S, 0.5)), 
              aes(x = S, ymin = lwr, ymax = upr), 
              fill = '#000', alpha = 0.1) +
  
  geom_line(data = pred_sp0_min %>% ungroup %>% group_by(S) %>% 
              summarise(lwr = quantile(prod_S, 0.025),
                        upr = quantile(prod_S, 0.975),
                        prod_S = quantile(prod_S, 0.5)), 
            aes(x = S, y = prod_S), 
            colour = '#000', alpha = 0.8, linewidth = 0.75) +
  
  # null productivity
  geom_ribbon(data = pred_prod0 %>% ungroup %>% group_by(S) %>%
                summarise(lwr = quantile(prod_S, 0.025),
                          upr = quantile(prod_S, 0.975),
                          prod_S = quantile(prod_S, 0.5)),
              aes(x = S, ymin = lwr, ymax = upr),
              fill = '#000', alpha = 0.1) +
  
  geom_line(data = pred_prod0 %>% ungroup %>% group_by(S) %>%
              summarise(lwr = quantile(prod_S, 0.025),
                        upr = quantile(prod_S, 0.975),
                        prod_S = quantile(prod_S, 0.5)),
            aes(x = S, y = prod_S),
            colour = 'grey30', alpha = 0.5, linewidth = 0.75, linetype = 'dashed') +
  
  # null abundance
  geom_ribbon(data = pred_abund0 %>% ungroup %>% group_by(S) %>% 
                summarise(lwr = quantile(prod_S, 0.025),
                          upr = quantile(prod_S, 0.975),
                          prod_S = quantile(prod_S, 0.5)), 
              aes(x = S, ymin = lwr, ymax = upr), 
              fill = 'darkred', alpha = 0.1) +
  
  geom_line(data = pred_abund0 %>% ungroup %>% group_by(S) %>% 
              summarise(lwr = quantile(prod_S, 0.025),
                        upr = quantile(prod_S, 0.975),
                        prod_S = quantile(prod_S, 0.5)), 
            aes(x = S, y = prod_S), 
            colour = 'darkred', alpha = 0.5, linewidth = 0.75) +

  theme(axis.text.x = element_blank()) + labs(x = NULL)

(f_spatial | (f_spatial2 / f_spatial3)) + plot_layout(ncol = 2, widths = c(2,1))


# ==================================================================================================
# Figure 3: Biogeographical patterns
# ==================================================================================================

unique(vars_single$biogeo)

# to sort in geographical order
regions <- c("w_indian", "e_indian", "central_indopacific", "s_china", "aus_indopacific", 
"central_pacific", "outer_pacific", "outer_e_pacific", "e_pacific",
"caribbean_w_atl", "wc_atlantic", "e_atlantic")

vars_single$biogeo <- factor(vars_single$biogeo, levels = regions)
pred_biogeo$biogeo <- factor(pred_biogeo$biogeo, levels = regions)
ecoR$region <- factor(ecoR$region, levels = regions)
ecoR <- ecoR %>% st_make_valid %>% st_break_antimeridian(., lon_0 = 160)

# biogeographical map

ggplot() +
  geom_sf(data = world, color='grey75', fill='grey95') +
  geom_sf(data = ecoR, 
          aes(fill = region), colour = '#333') +
  coord_sf(crs = proj, ylim = c(-4366274.76, 4366274.76)) +
  theme_minimal(base_family = 'Lato') + guides(fill = 'none') + 
  theme(panel.grid.minor.x = element_blank()) +
  scale_fill_cherulean(palette = 'cheridis', discrete = T)

f_geo <- ggplot(data = vars_single) +
  labs(x = 'S', y = bquote("Productivity (kg⋅ha"^-1~"⋅day"^-1~")")) +
  geom_point(aes(x = S, y = prod_kg, colour = biogeo), 
             shape = 21, size = 1.5, alpha = 0.5, fill = 'transparent') +
  
  geom_ribbon(data = pred_biogeo %>% ungroup %>% group_by(biogeo, S) %>% 
                summarise(lwr = quantile(prod_S, 0.025),
                          upr = quantile(prod_S, 0.975),
                          prod_S = quantile(prod_S, 0.5)), 
              aes(x = S, ymin = lwr, ymax = upr, group = biogeo), 
              fill = 'transparent', colour = '#000', linetype = 'dashed') +
  
  # overall S effect line
  geom_line(data = pred_spatial %>% ungroup %>% group_by(S) %>% 
              summarise(lwr = quantile(prod_S, 0.025),
                        upr = quantile(prod_S, 0.975),
                        prod_S = quantile(prod_S, 0.5)), 
            aes(x = S, y = prod_S), 
            colour = '#000', linewidth = 0.75, alpha = 0.4) +
  
  geom_line(data = pred_biogeo %>% ungroup %>% group_by(biogeo, S) %>% 
              summarise(lwr = quantile(prod_S, 0.025),
                        upr = quantile(prod_S, 0.975),
                        prod_S = quantile(prod_S, 0.5)), 
            aes(x = S, y = prod_S, group = biogeo), 
            colour = '#000', linewidth = 0.75) +

  scale_colour_cherulean(palette = 'cheridis', discrete = T) +
  scale_fill_cherulean(palette = 'cheridis', discrete = T) +
  scale_x_log10() + scale_y_log10() +
  guides(colour = 'none', fill = 'none') +
  style + facet_wrap(~biogeo)

  f_geomap
  f_geo

# ==================================================================================================
# Figure 4: Temporal trends and model
# ==================================================================================================

f_time <- ggplot(data = vars_time_long) +
  labs(x = "∆log S", y = bquote("∆log Productivity (kg⋅ha"^-1~"⋅day"^-1~")")) +
  geom_hline(aes(yintercept = 0), linetype = 'dashed', colour = '#AAA') +
  geom_vline(aes(xintercept = 0), linetype = 'dashed', colour = '#AAA') +
  geom_linerange(aes(xmin = S_slope - 1.96*S_se, 
                     xmax = S_slope + 1.96*S_se, 
                     y = prod_slope, colour = quad), linewidth = 0.2, alpha = 0.3) +
  geom_linerange(aes(ymin = prod_slope - 1.96*prod_se, 
                     ymax = prod_slope + 1.96*prod_se, 
                     x = S_slope, colour = quad), linewidth = 0.2, alpha = 0.3) +
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
            colour = '#000', linewidth = 0.5) +
  geom_abline(aes(intercept = 0, slope = 1), linetype = 'dashed', colour = '#444', alpha = 0.6) +
  guides(colour = 'none', fill = 'none') +
  scale_x_continuous(limits = c(-2.5, 1.5)) +
  scale_y_continuous(limits = c(-3.5, 3)) + style

# f_time_x <- ggplot(data = vars_time_long) +
#   geom_density(aes(x = S_slope), fill = '#CCC', colour = '#666') +
#   geom_vline(aes(xintercept = 0), linetype = 'dashed', colour = '#666') +
#   theme_void() + scale_x_continuous(limits = c(-2.5, 1.5))

# f_time_y <- ggplot(data = vars_time_long) +
#   geom_density(aes(y = prod_slope), fill = '#CCC', colour = '#666') +
#   geom_hline(aes(yintercept = 0), linetype = 'dashed', colour = '#666') +
#   theme_void() + scale_y_continuous(limits = c(-3.5, 3))

# f_time <- f_time_x + plot_spacer() + f_time1 + f_time_y + 
#   plot_layout(ncol = 2, nrow = 2, widths = c(5,1), heights = c(1,6))
# f_time


# Change trajectory bars ===========================================================================

ggplot(data = dt_tchange %>% filter(group == 'full')) +
  geom_col(aes(y = n, fill = quad, x = group), position = 'fill', colour = '#222', linewidth = 0.5) + 
  scale_fill_manual(values = c('#1c5dac77','#1c5dac','#60dcac77','#60dcac'), name = NULL,
                    labels = c('negative', 'negative 95% CI', 'positive', 'positive 95% CI')) +
  theme_minimal() + labs(x = NULL, y = NULL) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
  theme(panel.grid = element_blank(), legend.position = 'bottom')


f_time_slab <- ggplot() +
  stat_slab(data = vars_time_long %>% select(S_slope, prod_slope) %>% 
              pivot_longer(cols = everything(), names_to = 'metric', values_to = 'slopes'), 
            aes(x = slopes, y = metric), normalize = 'groups', 
            fill = '#fff', colour = '#222', linewidth = 0.5) + 
  geom_vline(aes(xintercept = 0), linetype = 'dashed', alpha = 0.5, colour = '#000') + style
  

f_time_bar <- ggplot(data = dt_temp %>% filter(group == 'full')) +
  geom_col_pattern(aes(x = metric, y = n, fill = direction, pattern = delta_certain), 
                   colour = '#222', position = 'fill', width = 0.5) + style +
  geom_hline(aes(yintercept = 0.5), linetype = 'dashed', alpha = 0.5, colour = '#000') +
  labs(x = NULL, y = NULL) +
  scale_pattern_discrete(choices = c('circle', 'stripe'), labels = c('NA', 'certain')) +
  scale_fill_manual(values = c('#3a89e3', '#3a89e355', '#FFF'), labels = c('decreasing', 'increasing', 'no change'))

f_time_slab / f_time_bar


ggplot() +
  labs(x = 'S', y = bquote("Productivity (kg⋅ha"^-1~"⋅day"^-1~")")) +
  # gradient background
  geom_point(data = vars_single, aes(x = S, y = prod_kg), size = 2, colour = '#666') +
  
  # temporal points
  geom_point(data = vars_time2, aes(x = S1, y = P1, fill = quad), 
             size = 2.5, shape = 21, alpha = 0.7, colour = 'grey20') +
  scale_fill_manual(values = c('#1c5dac','#60dcac')) +
  scale_x_log10() + scale_y_log10() + style + theme(legend.position = 'bottom') +
  facet_wrap(~quad)


# ==============================================================================
# Spatiotemporal figure
# ==============================================================================

point_time <- grid_sf %>% select(gridcell, geometry) %>% 
  mutate(geometry = st_centroid(geometry)) %>% # get centroids of each gridcell
  left_join(vars_time_long, ., by = 'gridcell') %>% 
  select(quad, delta_certain, geometry) %>% st_as_sf

point_time$dir <- str_extract(point_time$quad, '^[a-z]')

ggplot() +
  geom_sf(data = world, color='grey75', fill='grey95') +
  geom_sf(data = point_time %>% dplyr::filter(delta_certain == 0), 
          aes(colour = dir), shape = 21, size = 2, alpha = 0.2, fill = 'transparent') +
  geom_sf(data = point_time %>% dplyr::filter(delta_certain == 1), 
          aes(fill = dir), shape = 21, size = 2, alpha = 0.3, colour = '#222') +
  coord_sf(crs = proj, ylim = c(-4366274.76, 4366274.76)) +
  theme_minimal() + guides(colour = 'none') + 
  theme(legend.position = 'bottom', panel.grid = element_blank()) +
  scale_fill_manual(values = c('#1c5dac66','#60dcac66'), name = NULL,
                    labels = c('opposite', 'parallel')) +
  scale_colour_manual(values = c('#1c5dac','#60dcac')) +
  facet_grid(rows = vars(dir))

ggplot() +
  labs(x = 'S', y = bquote("Productivity (kg⋅ha"^-1~"⋅day"^-1~")")) +
  geom_point(data = vars_single, aes(x = S, y = prod_kg), size = 1, colour = '#bbb', alpha = 0.3) +

  geom_point(data = vars_time2, aes(x = meanS, y = meanP, fill = quad),
    size = 2, shape = 21, alpha = 0.5) +
  scale_fill_manual(values = c('#1c5dac66','#60dcac66'), name = NULL,
                    labels = c('opposite', 'parallel')) +
  scale_x_log10(breaks = single_x) + scale_y_log10(breaks = single_y) + style +
  facet_grid(rows = vars(quad))


