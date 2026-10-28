
# ProducDiv
# Figure code

require(ggplot2)
require(patchwork)
require(sf)
require(rnaturalearth)
require(ggdist)
require(ggpattern)

# colour palettes
source('https://gist.githubusercontent.com/cherfychow/e9ae890fd16f4c86730748c067feee2b/raw/218784603602ef3f738dee9006b8e7ae1a70e9e0/cherulean.R')
style <- theme_classic(base_size = 15, base_family = 'Lato') + 
  theme(panel.grid = element_blank(),
    axis.text = element_text(colour = '#222222'),
    axis.ticks = element_line(linewidth = 0.3, colour = '#222222'),
    axis.line = element_line(linewidth = 0.3, colour = '#222222'))

single_x <- c(1,5,10,50,100,150)
single_y <- c(0,0.001,.01, .1, 1, 10, 100, 500)

proj <- '+proj=eqearth +lon_0=160 +x_0=0 +y_0=0 +R=6371008.7714 +units=m +no_defs +type=crs'

world <- ne_coastline('medium', returnclass = 'sf') %>%
  st_make_valid(.) %>% st_break_antimeridian(., lon_0 = 160)
  
# ==================================================================================================
# Figure 1: Spatial pattern + model predictions
# ==================================================================================================

# Base single year data

f_spatial_dt <- ggplot() +
  labs(x = 'S', y = bquote("Productivity (kg⋅ha"^-1~"⋅day"^-1~")")) +
  
  stat_bin_2d(data = vars_single, aes(x = S, y = prod_kg, fill = after_stat(count)), bins = 100) +
  # geom_density_2d(data = vars_single, aes(x = S, y = prod_kg), colour = '#fff', bins = 8, linewidth = 0.3) +
  scale_fill_cherulean(palette = 'gomphosus', discrete = F, reverse = T, limits = c(1,200)) +
  
  scale_x_log10(breaks = single_x) + scale_y_log10(breaks = single_y) + 
  style

# heat map + spatial partial prediction
f_spatial <- f_spatial_dt + 
  geom_ribbon(data = pred_spatial %>% ungroup %>% group_by(S) %>% 
                summarise(lwr = quantile(prod_S, 0.025),
                          upr = quantile(prod_S, 0.975),
                          prod_S = quantile(prod_S, 0.5)), 
              aes(x = S, ymin = lwr, ymax = upr), 
              fill = 'transparent', colour = 'black', linetype = 'dashed') +
  
  geom_line(data = pred_spatial %>% ungroup %>% group_by(S) %>% 
              summarise(lwr = quantile(prod_S, 0.025),
                        upr = quantile(prod_S, 0.975),
                        prod_S = quantile(prod_S, 0.5)), 
            aes(x = S, y = prod_S), 
            colour = 'black', linewidth = 1) +
  guides(fill = "none")

# quantile regression
f_spatial2 <- f_spatial_dt + 
  geom_ribbon(data = pred_quant %>% filter(quantile != 'med'), 
    aes(x = S, ymin = lower, ymax = higher, group = quantile),
              fill = 'transparent', colour = 'black', linetype = 'dashed') +
  geom_line(data = pred_quant %>% filter(quantile != 'med'), 
    aes(x = S, y = fit, group = quantile), linewidth = .75) 

(f_spatial | f_spatial2)

# ==================================================================================================
# Figure 2: Temporal trends and model
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

  # predictions
  geom_ribbon(data = pred_time_full,
              aes(x = S_slope, ymin = prod_S_l, ymax = prod_S_u),
              colour = 'black', fill = 'transparent', linetype = 'dashed') +
  geom_line(data = pred_time_full,
            aes(x = S_slope, y = prod_S),
            colour = 'black', linewidth = 0.5) +
  geom_abline(aes(intercept = 0, slope = 1), linetype = 'dashed', colour = '#444', alpha = 0.6) +
  guides(colour = 'none', fill = 'none') +
  scale_x_continuous(limits = c(-2.5, 1.5)) +
  scale_y_continuous(limits = c(-3.5, 3)) + style

  f_time_N <- ggplot() +
  labs(x = "∆log N", y = bquote("∆log Productivity (kg⋅ha"^-1~"⋅day"^-1~")")) +
  geom_hline(aes(yintercept = 0), linetype = 'dashed', colour = '#AAA') +
  geom_vline(aes(xintercept = 0), linetype = 'dashed', colour = '#AAA') +
  geom_linerange(data = vars_time_long, 
                aes(
                  xmin = N_slope - 1.96*N_se, 
                  xmax = N_slope + 1.96*N_se, 
                  y = prod_slope, colour = quad), 
                linewidth = 0.2, alpha = 0.3) +
  geom_linerange(data = vars_time_long, aes(ymin = prod_slope - 1.96*prod_se, 
                     ymax = prod_slope + 1.96*prod_se, 
                     x = N_slope, colour = quad), linewidth = 0.2, alpha = 0.3) +
  scale_colour_manual(values = c('#1c5dac','#1c5dac','#60dcac','#60dcac')) +
  # uncertain points
  geom_point(data = vars_time_long %>% filter(!str_detect(quad, 'certain')),
             aes(x = N_slope, y = prod_slope), 
             size = 2, shape = 21, alpha = 0.3, colour = '#444', fill = '#FFF') +
  # certain slopes
  geom_point(data = vars_time_long %>% filter(str_detect(quad, 'certain')),
             aes(x = N_slope, y = prod_slope, fill = quad), 
             size = 2, shape = 21, colour = '#444444AA') +
  scale_fill_manual(values = c('#1c5dac','#60dcac')) +

  # partial prediction
  geom_ribbon(data = pred_time_full,
              aes(x = N_slope, ymin = prod_N_l, ymax = prod_N_u),
              colour = 'black', fill = 'transparent', linetype = 'dashed') +
  geom_line(data = pred_time_full,
            aes(x = N_slope, y = prod_N),
            colour = 'black', linewidth = 0.5) +
  geom_abline(aes(intercept = 0, slope = 1), linetype = 'dashed', colour = '#444', alpha = 0.6) +
  guides(colour = 'none', fill = 'none') +
  scale_y_continuous(limits = c(-3.5, 3)) + style

  f_time + f_time_N

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
  geom_vline(aes(xintercept = 0), linetype = 'dashed', alpha = 0.5, colour = 'black') + style
  

f_time_bar <- ggplot(data = dt_temp %>% filter(group == 'full')) +
  geom_col_pattern(aes(x = metric, y = n, fill = direction, pattern = delta_certain), 
                   colour = '#222', position = 'fill', width = 0.5) + style +
  geom_hline(aes(yintercept = 0.5), linetype = 'dashed', alpha = 0.5, colour = 'black') +
  labs(x = NULL, y = NULL) +
  scale_pattern_discrete(choices = c('circle', 'stripe'), labels = c('NA', 'certain')) +
  scale_fill_manual(values = c('#3a89e3', '#3a89e355', '#FFF'), labels = c('decreasing', 'increasing', 'no change'))

f_time_slab / f_time_bar


ggplot() +
  labs(x = 'S', y = bquote("Productivity (kg⋅ha"^-1~"⋅day"^-1~")")) +
  # gradient background
  geom_point(data = vars_single, aes(x = S, y = prod_kg), size = 2, colour = '#666') +
  
  # temporal points
  geom_point(data = vars_time, aes(x = S1, y = P1, fill = quad), 
             size = 2.5, shape = 21, alpha = 0.7, colour = 'grey20') +
  scale_fill_manual(values = c('#1c5dac','#60dcac')) +
  scale_x_log10() + scale_y_log10() + style + theme(legend.position = 'bottom') +
  facet_wrap(~quad)


# ==============================================================================
# Fig. 3: Spatiotemporal figure
# ==============================================================================

require(dggridR)
require(tidytable)
require(stringr)

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

  geom_point(data = vars_time, aes(x = meanS, y = meanP, fill = quad),
    size = 2, shape = 21, alpha = 0.5) +
  scale_fill_manual(values = c('#1c5dac66','#60dcac66'), name = NULL,
                    labels = c('opposite', 'parallel')) +
  scale_x_log10(breaks = single_x) + scale_y_log10(breaks = single_y) + style +
  facet_grid(rows = vars(quad))

# heat map alternative
grid_heat <- dgconstruct(area = 500^2, metric = TRUE, resround = 'nearest')
temp <- grid_sf %>% group_by(gridcell) %>%
  summarise(centroid = st_centroid(geometry))
temp$lon <- st_coordinates(temp$centroid)[,1]
temp$lat <- st_coordinates(temp$centroid)[,2]
temp$centroid <- NULL
summary(temp)

heat_match <- dgGEO_to_SEQNUM(dggs = grid_heat, in_lat_deg = temp$lat, in_lon_deg = temp$lon)
n_distinct(heat_match$seqnum)
# dggrid object to sf + cell subset
grid_heat <- dgcellstogrid(grid_heat, cells = sort(unique(heat_match$seqnum)), return_sf = T)

# match up for aggregating within cell
heat_time <- st_drop_geometry(grid_sf) %>% distinct(gridcell) %>%
  mutate(heatcellID = heat_match$seqnum) %>%
  left_join(., grid_heat, by = c('heatcellID' = 'seqnum'))

heat_time <- vars_time_long %>% select(seriesID, quad, delta_certain, gridcell) %>%
  left_join(., heat_time, by = 'gridcell') %>% 
  mutate(dir = str_extract(quad, '^[a-z]+')) %>% st_as_sf %>%
  # fix cells at the edges for plotting later
  st_wrap_dateline(.) %>% 
  st_break_antimeridian(., lon_0 = 160) %>% st_make_valid(.)

str(heat_time)

# generate standardised counts
cellcounts <- heat_time %>% ungroup %>% group_by(heatcellID) %>% 
  summarise(celltotal = n_distinct(seriesID)) %>%
  st_drop_geometry

heat_dir <- heat_time %>% ungroup %>% group_by(heatcellID, dir) %>%
  summarise(count = n_distinct(seriesID), geometry = first(geometry)) %>% 
  st_as_sf %>% ungroup
heat_dir <- left_join(heat_dir, cellcounts, by = 'heatcellID')
heat_dir$propn <- with(heat_dir, count/celltotal)

heat_dircert <- heat_time %>% ungroup %>% 
  group_by(heatcellID, dir, delta_certain) %>%
  dplyr::summarise(count = n_distinct(seriesID), geometry = first(geometry)) %>% 
  st_as_sf %>% ungroup
heat_dircert <- left_join(heat_dircert, cellcounts, by = 'heatcellID')
heat_dircert$propn <- with(heat_dircert, count/celltotal)

heat_dir <- st_break_antimeridian(heat_dir, lon_0 = 160)
heat_dircert <- st_break_antimeridian(heat_dircert, lon_0 = 160)
grid_heat <- grid_heat %>% filter(seqnum %in% heat_time$heatcellID)

proj <- '+proj=eqearth +lon_0=160 +x_0=0 +y_0=0 +R=6371008.7714 +units=m +no_defs +type=crs'

ggplot() +
  # base
  geom_sf(data = world, color='grey75', fill='grey95') +
  theme_minimal() +

  # all data
  geom_sf(data = grid_heat, fill = 'grey85', alpha = 0.8, colour = '#FFFFFF') +
  geom_sf(data = heat_dir, aes(fill = propn), alpha = 0.75) + 
  scale_fill_cherulean(palette = 'cheridis', discrete = F, name = 'Proportion of surveys') +

  theme(legend.position = 'bottom', panel.grid = element_blank()) +
  coord_sf(crs = proj, ylim = c(-4366274.76, 4366274.76)) +
  facet_grid(rows = vars(dir))
