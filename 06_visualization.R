# Figure 1 ----------------------------------------------------------------

# Load necessary libraries
library(phyloseq) #  For handling and analyzing microbiome data
library(tidyverse) # For data manipulation
library(ggthemes) # For visualization
library(patchwork) # For combining panels into multi-panel figures
library(sf)
library(ggnewscale)
library(ggspatial)
library(httr)
library(jsonlite)
library(rnaturalearth)
library(rnaturalearthdata)


# Panel A — sampling site map --------------------------------------------------

samples_meta_fig1 <- read.csv('data/metadata.csv') %>%
  mutate(site = tolower(site))
meta_sf <- sf::st_as_sf(
  samples_meta_fig1,
  coords = c('long', 'lat'),
  crs = 4326
)

fetch_overpass_coastline <- function(bbox, max_tries = 5) {
  query <- sprintf(
    '[out:json][timeout:25];way["natural"="coastline"](%f,%f,%f,%f); out geom;',
    bbox['ymin'],
    bbox['xmin'],
    bbox['ymax'],
    bbox['xmax']
  )
  for (i in seq_len(max_tries)) {
    resp <- httr::POST(
      'https://overpass-api.de/api/interpreter',
      body = list(data = query),
      encode = 'form'
    )
    txt <- httr::content(resp, 'text', encoding = 'UTF-8')
    j <- tryCatch(
      jsonlite::fromJSON(txt, simplifyVector = FALSE),
      error = function(e) NULL
    )
    if (!is.null(j)) {
      return(j)
    }
    Sys.sleep(5) # Overpass rate-limits repeated requests; back off and retry
  }
  stop(
    'Failed to retrieve coastline data from Overpass API after ',
    max_tries,
    ' attempts.'
  )
}

overpass_json <- fetch_overpass_coastline(c(
  xmin = 13.9600,
  ymin = 40.7280,
  xmax = 13.9680,
  ymax = 40.7350
))

coastline_lines <- lapply(overpass_json$elements, function(el) {
  if (is.null(el$geometry)) {
    return(NULL)
  }
  coords <- do.call(rbind, lapply(el$geometry, function(pt) c(pt$lon, pt$lat)))
  st_linestring(coords)
})
coastline_lines <- coastline_lines[!sapply(coastline_lines, is.null)]
coastline_sf <- st_sf(
  id = seq_along(coastline_lines),
  geometry = st_sfc(coastline_lines, crs = 4326)
)

map_bbox <- st_as_sfc(st_bbox(
  c(xmin = 13.9622, xmax = 13.9660, ymin = 40.7300, ymax = 40.7335),
  crs = 4326
))
coastline_merged <- st_union(coastline_sf) %>% st_line_merge()
coastline_closed <- st_union(coastline_merged, st_boundary(map_bbox))
land_water_polys <- st_polygonize(coastline_closed) %>%
  st_collection_extract('POLYGON')
land_water_sf <- st_sf(
  id = seq_along(land_water_polys),
  geometry = land_water_polys
) %>%
  mutate(area = as.numeric(st_area(geometry))) %>%
  mutate(type = ifelse(area == max(area), 'water', 'land'))

panel_A <- ggplot() +
  geom_sf(data = land_water_sf, aes(fill = type), color = NA) +
  scale_fill_manual(
    values = c(water = '#cfe7f5', land = '#c9a876'),
    guide = 'none'
  ) +
  new_scale_fill() +
  geom_sf(
    data = meta_sf,
    aes(fill = pH),
    shape = 21,
    size = 4,
    color = 'black'
  ) +
  geom_sf_text(
    data = meta_sf,
    aes(label = site),
    nudge_x = 0.00009,
    nudge_y = 0.000038,
    size = 3,
    hjust = 0
  ) +
  scale_fill_distiller(name = 'pH', palette = 'YlGnBu', direction = 1) +
  annotate(
    'text',
    x = 13.9645,
    y = 40.7332,
    label = 'Mediterranean Sea',
    fontface = 'italic',
    size = 4
  ) +
  annotation_north_arrow(
    location = 'tl',
    which_north = 'true',
    pad_x = unit(0.3, 'cm'),
    pad_y = unit(0.3, 'cm'),
    height = unit(1.0, 'cm'),
    width = unit(0.8, 'cm'),
    style = north_arrow_fancy_orienteering
  ) +
  annotation_scale(
    location = 'bl',
    width_hint = 0.32,
    unit_category = 'metric',
    pad_x = unit(0.3, 'cm'),
    pad_y = unit(0.3, 'cm')
  ) +
  scale_x_continuous(breaks = c(13.9625, 13.9635, 13.9645, 13.9655)) +
  coord_sf(
    xlim = c(13.9622, 13.9660),
    ylim = c(40.7300, 40.7335),
    expand = FALSE
  ) +
  theme_classic() +
  xlab('') +
  ylab('') +
  labs(tag = 'A')

italy_sf <- ne_countries(scale = 'large', returnclass = 'sf') %>%
  filter(name_long == 'Italy')
italy_inset <- ggplot(data = italy_sf) +
  geom_sf(fill = 'lightgrey', color = 'black') +
  geom_point(aes(x = 13.96, y = 40.73), shape = 21, size = 4, fill = 'red') +
  theme_void() +
  theme(panel.background = element_rect(fill = 'white', color = 'black'))

panel_A <- panel_A +
  inset_element(
    italy_inset,
    left = 0.65,
    bottom = 0.0,
    right = 1.0,
    top = 0.24,
    ignore_tag = TRUE
  )


#Load raw phyloseq obj (not-transformed data ) from Script 01
raw_data <- readRDS('data/phylo_obj_raw_dat.RDS')

#Load env data
env_data <- read.csv('data/metadata.csv') %>%
  mutate(site = tolower(site)) %>%
  select(site, pH)

#Load alpha div table from Script 01
alpha_tab <- readRDS('data/alpha_div_rarefied.RDS') %>%
  rownames_to_column('site') %>%
  left_join(env_data, by = 'site')


alpha_tab$transect <- ifelse(grepl("^n", alpha_tab$site), "North", "South")


#Derive % ASVs per sample fr each class

rel_asv_table <- psmelt(raw_data) %>%
  select(OTU, Sample, class, Abundance) %>%
  group_by(Sample, OTU) %>%
  mutate(tot_reads_OTU_sample = sum(Abundance)) %>%
  ungroup() %>%
  filter(tot_reads_OTU_sample > 0) %>%
  group_by(Sample, class) %>%
  summarise(OTU_count = n_distinct(OTU)) %>%
  ungroup() %>%
  group_by(Sample) %>%
  mutate(tot_OTU_count = sum(OTU_count)) %>%
  ungroup() %>%
  mutate(rel_OTU_count = OTU_count / tot_OTU_count)

#Generate plot (panel B)

panel_B <- ggplot(
  rel_asv_table,
  aes(
    x = factor(
      Sample,
      levels = c(
        'n1',
        'n2',
        'n3',
        'n4',
        's1',
        's2',
        's3',
        's4',
        's5',
        's6',
        's7',
        's8',
        's9',
        's10'
      )
    ),
    y = rel_OTU_count,
    fill = class
  )
) +
  geom_bar(stat = "identity", position = "stack", col = 'black') +
  scale_fill_tableau(palette = 'Tableau 10') +
  theme_classic() +
  labs(x = "", y = "% ASVs", fill = "Class", tag = 'B') +
  theme(
    axis.title = element_text(face = 'bold'),
    axis.text.x = element_text(size = 13, face = 'bold'),
    axis.text.y = element_text(size = 10),
    legend.text = element_text(size = 11),
    legend.title = element_blank(),
    legend.position = 'none'
  )

#Derive % reads per sample for each class

rel_reads_table <- psmelt(raw_data) %>%
  select(Sample, class, Abundance) %>%
  group_by(Sample, class) %>%
  summarise(total_reads_class = sum(Abundance), .groups = "drop") %>%
  group_by(Sample) %>%
  mutate(total_reads_sample = sum(total_reads_class)) %>%
  mutate(rel_read_abundance = total_reads_class / total_reads_sample) %>%
  ungroup()


#generate plot panel C

panel_C <- ggplot(
  rel_reads_table,
  aes(
    x = factor(
      Sample,
      levels = c(
        'n1',
        'n2',
        'n3',
        'n4',
        's1',
        's2',
        's3',
        's4',
        's5',
        's6',
        's7',
        's8',
        's9',
        's10'
      )
    ),
    y = rel_read_abundance,
    fill = class
  )
) +
  geom_bar(stat = "identity", position = "stack", col = 'black') +
  scale_fill_tableau(palette = 'Tableau 10') +
  guides(fill = guide_legend(nrow = 2)) +
  theme_classic() +
  labs(x = "", y = "% reads", fill = "Class", tag = 'C') +
  theme(
    axis.title = element_text(face = 'bold'),
    axis.text.x = element_text(size = 13, face = 'bold'),
    axis.text.y = element_text(size = 10),
    legend.text = element_text(size = 11),
    legend.title = element_blank(),
    legend.position = 'bottom'
  )


#generate plot panel D

panel_D <- ggplot(
  alpha_tab,
  aes(
    x = factor(
      site,
      levels = c(
        'n1',
        'n2',
        'n3',
        'n4',
        's1',
        's2',
        's3',
        's4',
        's5',
        's6',
        's7',
        's8',
        's9',
        's10'
      )
    ),
    y = diversity_shannon
  )
) +
  geom_line(aes(group = transect), color = "gray30", size = 0.6) +
  geom_point(aes(fill = pH), shape = 21, size = 5) +
  scale_fill_distiller(palette = "YlGnBu", direction = 1) +
  theme_minimal() +
  labs(x = "", y = "Shannon", fill = "Class", tag = 'D') +
  theme(
    axis.title = element_text(face = 'bold'),
    axis.text.x = element_text(size = 13, face = 'bold'),
    axis.text.y = element_text(size = 10),
    legend.text = element_text(size = 11),
    legend.title = element_blank(),
    legend.position = 'none'
  )

#combine plots
side_plot <- panel_B / panel_C / panel_D

figure_1 <- (panel_A | side_plot) +
  plot_layout(widths = c(1.2, 1))

ggsave(
  'figures/figure_1.png',
  figure_1,
  device = 'png',
  dpi = 300,
  width = 12,
  height = 8
)


# Figure 2 ----------------------------------------------------------------

# Okabe-Ito colorblind-safe palette for 3 clusters
cluster_cols <- c('C 1a' = '#0072B2', 'C 1b' = '#E69F00', 'C 2' = '#009E73')


# Load necessary libraries
library(vegan)
library(ggtree)
library(tidyverse)
library(patchwork)
library(ggrepel)
library(ggdist)


samples_meta <- read.csv('data/metadata.csv') %>%
  mutate(site = tolower(site)) %>%
  dplyr::select(1:4)

# Manually define sample groupings from cluster
clust_1a <- c('s2', 's4', 's3', 's5', 's1')
clust_1b <- c('s6', 's8', 'n1', 's9', 's10')
clust_2 <- c('n3', 's7', 'n4', 'n2')


samples_meta_cluster <- samples_meta %>%
  mutate(
    Cluster = ifelse(
      site %in% clust_1a,
      'C 1a',
      ifelse(site %in% clust_1b, 'C 1b', 'C 2')
    )
  )

#panel A

#load cluster res (Script 03)
cluster_res <- readRDS('data/hclust_res.RDS')


p <- ggtree(cluster_res, size = 0.6) %>%
  ggtree::rotate(15) %>%
  ggtree::rotate(16) +
  geom_tiplab(vjust = -1.4, hjust = -0.1, size = 3)

#add cluster metadata to color-code samples
final_cluster_plot <- p %<+%
  samples_meta_cluster +
  layout_dendrogram() +
  geom_tippoint(size = 3, aes(fill = Cluster), shape = 21) +
  scale_fill_manual(values = cluster_cols) +
  theme_dendrogram() +
  theme(
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14, face = 'bold')
  )

final_cluster_plot
#panel B

#load pcoa res (Scrip 03)

pcoa_res <- readRDS('data/PcoA_res.RDS')

#load envfit res (Scrip 03)
env_fit_res <- readRDS('data/envfit_res.RDS')

ordiplot(pcoa_res)
plot(env_fit_res, p.max = 0.05)

#extract pcoa points coordinates

pcoa_coords <- pcoa_res$points %>%
  as.data.frame() %>%
  rownames_to_column('site')

#extract envfit coordinates
arrow_factor <- ordiArrowMul(env_fit_res)
spp.scrs <- as.data.frame(vegan::scores(env_fit_res, display = "vectors")) *
  arrow_factor
spp.scrs <- cbind(
  spp.scrs,
  site = rownames(spp.scrs),
  Pvalues = env_fit_res$vectors$pvals,
  R_squared = env_fit_res$vectors$r
)

spp.scrs <- subset(spp.scrs, Pvalues < 0.05)


# generate convex hulls
pcoa_coords_plot <- pcoa_coords %>%
  mutate(
    Cluster = case_when(
      site %in% clust_1a ~ 'C 1a',
      site %in% clust_1b ~ 'C 1b',
      TRUE ~ 'C 2'
    )
  )


pcoa_coords_chull <- pcoa_coords_plot %>%
  group_by(Cluster) %>%
  slice(chull(V1, V2))

pcoa_plot <- ggplot(pcoa_coords_plot) +
  geom_polygon(
    data = pcoa_coords_chull,
    aes(x = V1, y = V2, fill = Cluster),
    alpha = 0.25,
    inherit.aes = FALSE
  ) +
  geom_point(
    aes(x = V1, y = V2, fill = Cluster),
    shape = 21,
    color = 'black',
    size = 3,
    stroke = 0.4
  ) +
  geom_segment(
    data = spp.scrs,
    aes(x = 0, xend = Dim1, y = 0, yend = Dim2),
    arrow = arrow(length = unit(0.2, "cm"), type = 'closed'),
    linewidth = 0.5,
    color = 'black'
  ) +
  geom_label_repel(
    data = spp.scrs,
    aes(x = Dim1, y = Dim2, label = site),
    size = 2.5,
    label.padding = unit(0.15, 'cm'),
    label.size = 0.5,
    fill = 'white'
  ) +
  scale_fill_manual(values = cluster_cols) +
  coord_fixed() +
  annotate("text", x = -25, y = -60, label = "r² = 0.46", size = 3.5) +
  xlab('Dim 1') +
  ylab('Dim 2') +
  theme_classic() +
  theme(legend.position = 'none')


panel_c <- ggplot(
  samples_meta_cluster,
  aes(fill = Cluster, color = Cluster, x = pH)
) +
  scale_color_manual(values = cluster_cols) +
  scale_fill_manual(values = cluster_cols) +
  stat_slab(alpha = .3) +
  stat_pointinterval(
    position = position_dodgejust(width = .2),
    justification = 0.1
  ) +
  scale_y_continuous(breaks = NULL) +
  geom_vline(xintercept = 7.7, linetype = 3) +
  ylab('') +
  theme_bw() +
  theme(legend.position = 'none')


(final_cluster_plot + pcoa_plot) /
  panel_c +
  plot_annotation(tag_levels = 'A') +
  plot_layout(guides = "collect")


ggsave('figures/figure_2.png', device = 'png', dpi = 300, width = 8, height = 7)

# Figure 3 ----------------------------------------------------------------
# Load necessary libraries
library(phyloseq)
library(tidyverse)
library(patchwork)
library(ggpubr)
library(ggstatsplot)
library(janitor)


# Molecular data

#load row reads counts (Script 01)

raw_data <- readRDS('data/phylo_obj_raw_dat.RDS')


#load sample metadata
samples_meta <- read.csv('data/metadata.csv') %>%
  mutate(site = tolower(site)) %>%
  dplyr::select(1:4)

# Manually define sample groupings from cluster
clust_1a <- c('s2', 's4', 's3', 's5', 's1')
clust_1b <- c('s6', 's8', 'n1', 's9', 's10')
clust_2 <- c('n3', 's7', 'n4', 'n2')

samples_meta_cluster <- samples_meta %>%
  mutate(
    Cluster = ifelse(
      site %in% clust_1a,
      'C 1a',
      ifelse(site %in% clust_1b, 'C 1b', 'C 2')
    )
  ) %>%
  rename('Sample' = 'site')


#calc n calcareous asvs per cluster

n_asvs_cluster_table <- psmelt(raw_data) %>%
  select(OTU, Sample, class, Abundance, shell_type) %>%
  group_by(Sample, OTU) %>%
  mutate(tot_reads_OTU_sample = sum(Abundance)) %>%
  ungroup() %>%
  filter(tot_reads_OTU_sample > 0) %>%
  select(OTU, Sample, , shell_type) %>%
  left_join(samples_meta_cluster, by = 'Sample') %>%
  filter(shell_type == 'calcareous') %>%
  group_by(Sample) %>%
  mutate(OTU_count = n_distinct(OTU)) %>%
  ungroup() %>%
  select(Sample, OTU_count, Cluster) %>%
  distinct()

n_asvs_cluster_table$Cluster <- factor(
  n_asvs_cluster_table$Cluster,
  levels = c('C 1a', 'C 1b', 'C 2')
)

#perform paired test with fdr correction
ggbetweenstats(
  n_asvs_cluster_table,
  Cluster,
  OTU_count,
  plot.type = "box",
  type = "np",
  pairwise.annotation = "p.value",
  pairwise.display = "significant", ## display only significant pairwise comparisons
  p.adjust.method = "fdr", ## adjust p-values for multiple tests using this method
  ylab = 'n Calcareous ASVs',
  results.subtitle = TRUE
)


#select only significant comparison

my_comparisons <- list(c("C 2", "C 1a"))

#visualize res

asv_count <- ggboxplot(
  n_asvs_cluster_table,
  x = "Cluster",
  y = "OTU_count",
  fill = "Cluster",
  palette = "viridis",
  width = 0.5
) +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif") + # Add pairwise comparisons p-value
  # stat_compare_means(label.y = 90)+
  stat_compare_means(label.y = 80, , size = 3.5) +
  scale_fill_manual(values = cluster_cols) +
  # geom_jitter()+
  ylab('n Calcareous ASVs') +
  xlab('') +
  theme_bw() +
  theme(
    axis.title = element_text(face = 'bold'),
    axis.text.x = element_text(size = 13, face = 'bold'),
    axis.text.y = element_text(size = 10),
    legend.position = 'none'
  )


#calc % calcareous asvs per cluster

relab_cluster_table <- psmelt(raw_data) %>%
  select(OTU, Sample, Abundance, shell_type) %>%
  left_join(samples_meta_cluster, by = 'Sample') %>%
  filter(!shell_type == 'not_known') %>%
  group_by(Sample) %>%
  mutate(tot_reads_by_sample = sum(Abundance)) %>%
  ungroup() %>%
  group_by(Sample, shell_type) %>%
  mutate(reads_by_sample_shell_type = sum(Abundance)) %>%
  ungroup() %>%
  filter(shell_type == 'calcareous') %>%
  mutate(
    rel_reads_calcareus = (reads_by_sample_shell_type / tot_reads_by_sample) *
      100
  ) %>%
  select(Sample, rel_reads_calcareus, Cluster) %>%
  distinct()


relab_cluster_table$Cluster <- factor(
  relab_cluster_table$Cluster,
  levels = c('C 1a', 'C 1b', 'C 2')
)

#perform paired test with fdr correction
ggbetweenstats(
  relab_cluster_table,
  Cluster,
  rel_reads_calcareus,
  plot.type = "box",
  type = "np",
  pairwise.display = "significant", ## display only significant pairwise comparisons
  p.adjust.method = "fdr", ## adjust p-values for multiple tests using this method
  ylab = '% Calcareous reads'
)

#select only significant comparison
my_comparisons <- list(c("C 2", "C 1a"))


#visualize res
reads_plot <- ggboxplot(
  relab_cluster_table,
  x = "Cluster",
  y = "rel_reads_calcareus",
  fill = "Cluster",
  palette = "viridis",
  width = 0.5
) +
  scale_fill_manual(values = cluster_cols) +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif") + # Add pairwise comparisons p-value
  # stat_compare_means(label.y = 90)+
  stat_compare_means(label.y = 60, size = 3.5) +
  ylab('% Calcareous reads') +
  xlab('') +
  theme_bw() +
  theme(
    axis.title = element_text(face = 'bold'),
    axis.text.x = element_text(size = 13, face = 'bold'),
    axis.text.y = element_text(size = 10),
    legend.position = 'none'
  )


reads_plot

# Morphological data

#load morph data
morph_data <- read.csv('data/morph_data/morph_counts.csv', header = TRUE) %>%
  clean_names() %>% #clear column names
  mutate(site = tolower(site)) %>% #make sample name consistent
  mutate(calcareous_tests = hy + porcel) #sum  %procellanaceous and hyaline

#add cluster info

morph_data_cluster <- morph_data %>%
  mutate(
    cluster = ifelse(
      site %in% clust_1a,
      'C 1a',
      ifelse(site %in% clust_1b, 'C 1b', 'C 2')
    )
  )


morph_data_cluster$cluster <- factor(
  morph_data_cluster$cluster,
  levels = c('C 1a', 'C 1b', 'C 2')
)

#perform paired test with fdr correction
ggbetweenstats(
  morph_data_cluster,
  cluster,
  calcareous_tests,
  plot.type = "box",
  type = "np",
  pairwise.display = "significant", ## display only significant pairwise comparisons
  p.adjust.method = "fdr", ## adjust p-values for multiple tests using this method
  ylab = 'calcareous_tests'
)


#select only significant comparison
my_comparisons <- list(c("C 2", "C 1a"))

#visualize res
rel_calc_specimens <- ggboxplot(
  morph_data_cluster,
  x = "cluster",
  y = "calcareous_tests",
  fill = "cluster",
  width = 0.5
) +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif") + # Add pairwise comparisons p-value
  # stat_compare_means(label.y = 90)+
  scale_fill_manual(values = cluster_cols) +
  stat_compare_means(label.y = 115, size = 3.5) +
  ylab('% Calcareous specimens') +
  xlab('') +
  theme_bw() +
  theme(
    axis.title = element_text(face = 'bold'),
    axis.text.x = element_text(size = 13, face = 'bold'),
    axis.text.y = element_text(size = 10),
    legend.position = 'none'
  )


#perform paired test with fdr correction  (Fragmentation index)

ggbetweenstats(
  morph_data_cluster,
  cluster,
  fi,
  plot.type = "box",
  type = "np",
  pairwise.display = "significant", ## display only significant pairwise comparisons
  p.adjust.method = "fdr", ## adjust p-values for multiple tests using this method
  ylab = 'FI'
)

#select only significant comparison
my_comparisons <- list(c("C 2", "C 1a"))

#visualize res
fi_plot <- ggboxplot(
  morph_data_cluster,
  x = "cluster",
  y = "fi",
  fill = "cluster",
  width = 0.5
) +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif") + # Add pairwise comparisons p-value
  # stat_compare_means(label.y = 90)+
  stat_compare_means(label.y = 90, size = 3.5) +
  scale_fill_manual(values = cluster_cols) +
  ylab('FI') +
  xlab('') +
  theme_bw() +
  theme(
    axis.title = element_text(face = 'bold'),
    axis.text.x = element_text(size = 13, face = 'bold'),
    axis.text.y = element_text(size = 10),
    legend.position = 'none'
  )


#perform paired test with fdr correction  (% Living Specimens)

ggbetweenstats(
  morph_data_cluster,
  cluster,
  living,
  plot.type = "box",
  type = "np",
  pairwise.display = "significant", ## display only significant pairwise comparisons
  p.adjust.method = "fdr", ## adjust p-values for multiple tests using this method
  ylab = 'Living'
)

#select only significant comparison
my_comparisons <- list(c("C 2", "C 1a"))


#visualize res
living_plot <- ggboxplot(
  morph_data_cluster,
  x = "cluster",
  y = "living",
  fill = "cluster",
  width = 0.5
) +
  # stat_compare_means()+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif") + # Add pairwise comparisons p-value
  # stat_compare_means(label.y = 90)+
  scale_fill_manual(values = cluster_cols) +
  stat_compare_means(label.y = 110, size = 3.5) +
  ylab('% Living specimens') +
  xlab('') +
  theme_bw() +
  theme(
    axis.title = element_text(face = 'bold'),
    axis.text.x = element_text(size = 13, face = 'bold'),
    axis.text.y = element_text(size = 10),
    legend.position = 'none'
  )


#perform paired test with fdr correction  (% Foraminiferal Density)

ggbetweenstats(
  morph_data_cluster,
  cluster,
  fd_d_l,
  plot.type = "box",
  type = "np",
  pairwise.display = "significant", ## display only significant pairwise comparisons
  p.adjust.method = "fdr", ## adjust p-values for multiple tests using this method
  ylab = 'FD(D+L)'
)

#select only significant comparison
my_comparisons <- list(c("C 2", "C 1a"))


#visualize res
fd_dl_plot <- ggboxplot(
  morph_data_cluster,
  x = "cluster",
  y = "fd_d_l",
  fill = "cluster",
  width = 0.5
) +
  # stat_compare_means()+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif") + # Add pairwise comparisons p-value
  # stat_compare_means(label.y = 90)+
  scale_fill_manual(values = cluster_cols) +
  stat_compare_means(label.y = 160, , size = 3.5) +
  ylab('FD(D+L)') +
  xlab('') +
  theme_bw() +
  theme(
    axis.title = element_text(face = 'bold'),
    axis.text.x = element_text(size = 13, face = 'bold'),
    axis.text.y = element_text(size = 10),
    legend.position = 'none'
  )

#combine plots
asv_count +
  reads_plot +
  rel_calc_specimens +
  fi_plot +
  living_plot +
  fd_dl_plot +
  plot_annotation(tag_levels = 'A')


ggsave('figures/figure_3.png', device = 'png', dpi = 300, height = 8, width = 9)


# Figure 4 ----------------------------------------------------------------

# Colorblind-safe colors for pH map categories
ph_map_cols <- c('C1a range' = '#D55E00', 'below 7.7' = '#E69F00')
# Load necessary libraries
library("rnaturalearth")
library("rnaturalearthdata")
library(tidyverse)

#Chronic vs episodic acidification in the Mediterranean (2090–2100)
#
# Rationale: Bio-ORACLE ph_min = decadal minimum pH; ph_mean = decadal mean pH.
# Since mean >= min always, classifying cells by which statistic crosses the
# empirical Ischia thresholds yields ecologically distinct exposure regimes:
#
#   Chronic  → ph_mean <= threshold: mean conditions already at Ischia-like levels
#   Episodic → only ph_min <= threshold: acidification reaches threshold
#              seasonally/episodically but not on average
#   Neither  → both above threshold
#

library(tidyverse)
library(rnaturalearth)
library(rnaturalearthdata)
library(patchwork)


# Shared style ----------------------------------------------------------------

journal_theme <- theme_classic(base_size = 7) +
  theme(
    axis.title = element_text(size = 7, face = 'bold'),
    axis.text = element_text(size = 6, color = 'black'),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 7, face = 'bold'),
    strip.text = element_text(size = 6),
    strip.background = element_rect(fill = 'grey92', color = NA),
    plot.title = element_text(size = 7, face = 'bold'),
    panel.border = element_rect(color = 'black', fill = NA, linewidth = 0.4)
  )

map_theme <- journal_theme +
  theme(
    panel.background = element_rect(fill = '#d6eaf8', color = NA),
    panel.grid.major = element_line(color = 'white', linewidth = 0.2),
    legend.position = 'bottom',
    legend.key.size = unit(0.35, 'cm'),
    legend.direction = 'horizontal'
  )

# Colorblind-safe 3-category palette (Okabe-Ito derived)
exposure_cols <- c(
  'Chronic' = '#D55E00', # vermillion — persistent, most concerning
  'Episodic' = '#E69F00', # amber — seasonal/event-driven
  'Neither' = '#CCCCCC' # light grey — not shown (filtered out below)
)


# Load and prepare data -------------------------------------------------------

# ph_min projections (Script 04 output)
ph_min_data <- readRDS('data/points_biooracle_ph_min.RDS') %>%
  filter(benthic_level == 'depth max') %>%
  filter(decade == '2090-2100') %>%
  select(x, y, scenario_label, ph_vals_min = ph_vals)

# ph_mean projections (Script 05 output)
ph_mean_data <- readRDS('data/points_biooracle_ph_mean.RDS') %>%
  filter(benthic_level == 'depth max') %>%
  filter(decade == '2090-2100') %>%
  select(x, y, scenario_label, ph_vals_mean = ph_vals)


# Full grid for the decade: include all cells from ph_min (which covers more
# area since min pH will always be <= mean pH)
full_grid <- readRDS('data/points_biooracle_ph_min.RDS') %>%
  filter(benthic_level == 'depth max') %>%
  filter(decade == '2090-2100') %>%
  select(x, y, scenario_label, predicted, ph_vals) %>%
  rename(ph_vals_min = ph_vals)

# Join mean classifications onto the full grid (NA where mean does not reach threshold)
combined <- full_grid %>%
  left_join(
    ph_mean_data,
    by = c('x', 'y', 'scenario_label')
  ) %>%
  mutate(
    exposure = case_when(
      !is.na(ph_vals_mean) ~ 'Chronic', # mean pH also within threshold range
      !is.na(ph_vals_min) ~ 'Episodic', # only min pH within threshold range
      TRUE ~ 'Neither'
    )
  ) %>%
  filter(exposure != 'Neither') # plot only cells within threshold range


# Scenario label ordering (mild to severe)
scenario_order <- c(
  'ssp1 2.6',
  'ssp2 4.5',
  'ssp4 6.0',
  'ssp3 7.0',
  'ssp5 8.5'
)
combined$scenario_label <- factor(
  combined$scenario_label,
  levels = scenario_order
)


# Individual scenario maps ----------------------------------------------------

world <- ne_countries(scale = "large", returnclass = "sf")

# Helper: one map panel per scenario; legend only on the last panel.
# Axis tick labels are shown only on the left column (latitude) and the
# bottom-most map of each column (longitude) -- repeating them on every
# panel wastes space that could otherwise go to the map itself.
make_scenario_map <- function(
  scen,
  show_legend = FALSE,
  show_y_axis = FALSE,
  show_x_axis = FALSE
) {
  dat <- combined %>% filter(scenario_label == scen)

  p <- ggplot(data = world) +
    geom_sf(fill = 'grey85', color = 'grey60', linewidth = 0.2) +
    coord_sf(xlim = c(-10, 40), ylim = c(30, 46), expand = FALSE) +
    geom_point(
      data = dat,
      aes(x = x, y = y, color = exposure),
      shape = 16,
      size = 0.9,
      alpha = 0.85
    ) +
    scale_color_manual(
      values = exposure_cols,
      labels = c(
        'Chronic' = 'Chronic (mean pH ≤ 7.7)',
        'Episodic' = 'Episodic (only min pH ≤ 7.7)'
      )
    ) +
    labs(color = 'Exposure regime', title = scen) +
    xlab('') +
    ylab('') +
    map_theme +
    guides(color = guide_legend(override.aes = list(size = 3)))

  if (!show_legend) {
    p <- p + theme(legend.position = 'none')
  }
  if (!show_y_axis) {
    p <- p +
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  }
  if (!show_x_axis) {
    p <- p +
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  }

  p
}

# Layout is "AB/CD/EF": scenarios 1,3,5 sit in the left column (A,C,E),
# 2,4 in the right column (B,D). E (5) and D (4) are the bottom-most map
# in their respective columns (F is the barplot, not a map).
map_panels <- lapply(
  seq_along(scenario_order),
  function(i) {
    make_scenario_map(
      scenario_order[i],
      show_legend = (i == 5),
      show_y_axis = (i %in% c(1, 3, 5)),
      show_x_axis = (i %in% c(4, 5))
    )
  }
)

# Summary table ---------------------------------------------------------------
# Area estimates per scenario and exposure regime.
#
# Each Bio-ORACLE cell is 0.05° × 0.05°. Cell area varies with latitude:
#   area_km2 = (0.05 × 111.32)^2 × cos(lat_rad)
# ≈ 26–31 km² across the Mediterranean (30–46°N).
# Summing cell areas gives a rough but reasonable total area estimate.

combined <- combined %>%
  mutate(cell_area_km2 = (0.05 * 111.32)^2 * cos(y * pi / 180))

exposure_summary <- combined %>%
  group_by(scenario_label, exposure) %>%
  summarise(
    n_cells = n(),
    area_km2 = round(sum(cell_area_km2)),
    .groups = 'drop'
  ) %>%
  pivot_wider(
    names_from = exposure,
    values_from = c(n_cells, area_km2),
    values_fill = 0
  ) %>%
  mutate(
    total_area_km2 = area_km2_Chronic + area_km2_Episodic,
    pct_chronic_area = round(area_km2_Chronic / total_area_km2 * 100, 1)
  ) %>%
  select(
    scenario_label,
    n_cells_Chronic,
    n_cells_Episodic,
    area_km2_Chronic,
    area_km2_Episodic,
    total_area_km2,
    pct_chronic_area
  )

print(exposure_summary)


# Panel G — exposure area barplot ---------------------------------------------

exposure_long <- exposure_summary %>%
  select(scenario_label, area_km2_Chronic, area_km2_Episodic) %>%
  pivot_longer(
    cols = c(area_km2_Chronic, area_km2_Episodic),
    names_to = 'exposure',
    values_to = 'area_km2'
  ) %>%
  mutate(
    exposure = str_remove(exposure, 'area_km2_'),
    exposure = factor(exposure, levels = c('Episodic', 'Chronic')),
    scenario_label = factor(scenario_label, levels = scenario_order),
    area_km2_k = area_km2 / 1000 # express in thousands of km²
  )

panel_G <- ggplot(
  exposure_long,
  aes(x = scenario_label, y = area_km2_k, fill = exposure)
) +
  geom_bar(
    stat = 'identity',
    position = 'stack',
    color = 'black',
    linewidth = 0.2
  ) +
  scale_fill_manual(
    values = exposure_cols,
    labels = c('Chronic' = 'Chronic', 'Episodic' = 'Episodic')
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(
    x = '',
    y = expression('Area (×10'^3 * ' km²)'),
    fill = 'Exposure regime'
  ) +
  journal_theme +
  theme(
    legend.position = 'none',
    axis.text.x = element_text(angle = 35, hjust = 1)
  )


# Combine panels --------------------------------------------------------------

layout <- "
AB
CD
EF
"

# 5 maps (A-E) fill the first 5 slots; barplot (F) fills the natural gap
figure_4 <- wrap_plots(c(map_panels, list(panel_G))) +
  plot_layout(design = layout) +
  plot_annotation(tag_levels = 'A')


# Export figure 4

ggsave(
  filename = 'figures/figure_4.png',
  plot = figure_4,
  width = 183, # mm — Nature/Scientific Reports double-column max width
  height = 170, # mm — Nature/Scientific Reports max figure height
  units = 'mm',
  dpi = 300,
  bg = 'white'
)


# Fig S1  -----------------------------

ggplot(data = world) +
  geom_sf(fill = 'grey85', color = 'grey60', linewidth = 0.2) +
  coord_sf(xlim = c(-10, 40), ylim = c(30, 46), expand = FALSE) +
  geom_point(
    data = combined %>% filter(exposure == 'Chronic'),
    aes(x = x, y = y),
    color = '#D55E00',
    shape = 16,
    size = 0.5,
    alpha = 0.85
  ) +
  facet_wrap(~scenario_label, ncol = 2) +
  labs(title = 'Chronic exposure only (mean pH ≤ threshold) — 2090–2100') +
  xlab('') +
  ylab('') +
  map_theme +
  theme(legend.position = 'none')


ggsave(
  filename = 'figures/figure_S1.png',
  # plot = figure_4,
  width = 180, # mm — full page width for most journals (e.g. Nature, Science)
  height = 140, # mm — adjust if map panels feel too compressed
  units = 'mm',
  dpi = 300,
  bg = 'white'
)
