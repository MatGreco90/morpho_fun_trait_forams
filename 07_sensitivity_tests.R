# Sensitivity tests
#
# Two robustness checks that arose from peer review, now part of the
# permanent analysis pipeline:
#
#   Part 1 — Does the calcareous-ASV cluster pattern (Fig. 3A-B) hold when
#            the two lowest taxonomic-resolution stations (N1, S1; >90%
#            unassigned reads) are excluded?
#
#   Part 2 — Is the projected pH threshold-crossing (Fig. 4; pH <= 7.7 by
#            2090-2100) robust to inter-model (CMIP6 ensemble) disagreement
#            in the Bio-ORACLE v3.0 projections?
#
# Part 1 depends on Script 01's output (data/phylo_obj_raw_dat.RDS) and
# Script 03's cluster assignment. Part 2 is fully self-contained and
# downloads its own data directly from Bio-ORACLE.


library(phyloseq)
library(tidyverse)
library(ggpubr)
library(patchwork)
library(gridExtra)
library(terra)
library(biooracler)

cluster_cols <- c('C 1a' = '#0072B2', 'C 1b' = '#E69F00', 'C 2' = '#009E73')

journal_theme <- theme_classic(base_size = 9) +
  theme(
    axis.title = element_text(size = 9, face = 'bold'),
    axis.text = element_text(size = 8, color = 'black'),
    legend.position = 'none',
    plot.title = element_text(size = 9, face = 'italic', hjust = 0.5),
    plot.tag = element_text(size = 10, face = 'bold')
  )


# ==============================================================================
# Part 1 — Taxonomic-resolution sensitivity: Fig. 3A-B excluding N1 & S1
# ==============================================================================
#
# N1 and S1 have the lowest taxonomic-assignment rates in the dataset (91.5%
# and 92.4% of reads unassigned, respectively). Since shell-type (used to
# identify calcareous ASVs) can only be annotated for taxonomically assigned
# ASVs, these two stations contribute proportionally little information to
# the calcareous-ASV comparisons. This section tests whether excluding them
# changes the reported cluster differences.

raw_data <- readRDS('data/phylo_obj_raw_dat.RDS')

samples_meta <- read.csv('data/metadata.csv') %>%
  mutate(site = tolower(site)) %>%
  dplyr::select(1:4)

# Cluster assignment from Script 03 (as used in 06_visualization.R Fig. 3)
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

melted <- psmelt(raw_data)


## 1a — Per-sample taxonomic / shell-type resolution -------------------------

resolution_tab <- melted %>%
  filter(Abundance > 0) %>%
  group_by(Sample) %>%
  summarise(
    Total_reads = sum(Abundance),
    Total_ASVs = n_distinct(OTU),
    Assigned_reads = sum(Abundance[shell_type != 'not_known']),
    Assigned_ASVs = n_distinct(OTU[shell_type != 'not_known']),
    .groups = 'drop'
  ) %>%
  mutate(
    `% reads assigned` = round(100 * Assigned_reads / Total_reads, 1),
    `% ASVs assigned` = round(100 * Assigned_ASVs / Total_ASVs, 1)
  ) %>%
  left_join(samples_meta_cluster, by = 'Sample') %>%
  dplyr::select(
    Sample, Cluster, Total_reads, Total_ASVs,
    Assigned_reads, Assigned_ASVs, `% reads assigned`, `% ASVs assigned`
  ) %>%
  arrange(`% reads assigned`)

write.csv(resolution_tab, 'data/sensitivity_table_taxonomic_resolution.csv', row.names = FALSE)

table_resolution_plot <- tableGrob(
  resolution_tab,
  rows = NULL,
  theme = ttheme_default(base_size = 9, core = list(fg_params = list(hjust = 0.5, x = 0.5)))
)

ggsave(
  'figures/sensitivity_table_taxonomic_resolution.png',
  table_resolution_plot,
  width = 220, height = 140, units = 'mm', dpi = 300, bg = 'white'
)

cat('\n=== Per-sample taxonomic resolution ===\n')
print(resolution_tab, n = Inf)


## 1b — Fig. 3A/B Kruskal-Wallis comparisons, with and without N1 & S1 -------
# (identical logic to 06_visualization.R Fig. 3A-B, parameterised over an
# excluded-sample set)

build_n_calcareous_asvs <- function(exclude = character(0)) {
  melted %>%
    filter(!Sample %in% exclude) %>%
    select(OTU, Sample, class, Abundance, shell_type) %>%
    group_by(Sample, OTU) %>%
    mutate(tot_reads_OTU_sample = sum(Abundance)) %>%
    ungroup() %>%
    filter(tot_reads_OTU_sample > 0) %>%
    select(OTU, Sample, shell_type) %>%
    left_join(samples_meta_cluster, by = 'Sample') %>%
    filter(shell_type == 'calcareous') %>%
    group_by(Sample) %>%
    mutate(OTU_count = n_distinct(OTU)) %>%
    ungroup() %>%
    select(Sample, OTU_count, Cluster) %>%
    distinct() %>%
    mutate(Cluster = factor(Cluster, levels = c('C 1a', 'C 1b', 'C 2')))
}

build_pct_calcareous_reads <- function(exclude = character(0)) {
  melted %>%
    filter(!Sample %in% exclude) %>%
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
    mutate(rel_reads_calcareus = (reads_by_sample_shell_type / tot_reads_by_sample) * 100) %>%
    select(Sample, rel_reads_calcareus, Cluster) %>%
    distinct() %>%
    mutate(Cluster = factor(Cluster, levels = c('C 1a', 'C 1b', 'C 2')))
}

low_res_stations <- c('n1', 's1')

datasets <- list(
  full_A = build_n_calcareous_asvs(character(0)),
  full_B = build_pct_calcareous_reads(character(0)),
  excl_A = build_n_calcareous_asvs(low_res_stations),
  excl_B = build_pct_calcareous_reads(low_res_stations)
)

kw_summary <- tibble(
  panel = c('A: n Calcareous ASVs', 'B: % Calcareous reads', 'A: n Calcareous ASVs', 'B: % Calcareous reads'),
  dataset = c('All 14 samples', 'All 14 samples', 'Excl. N1 & S1 (n=12)', 'Excl. N1 & S1 (n=12)'),
  data_key = names(datasets)
) %>%
  rowwise() %>%
  mutate(
    value_col = ifelse(grepl('^A', panel), 'OTU_count', 'rel_reads_calcareus'),
    kw = list(kruskal.test(datasets[[data_key]][[value_col]] ~ datasets[[data_key]]$Cluster))
  ) %>%
  mutate(
    chi_sq = round(kw$statistic, 3),
    df = kw$parameter,
    p_value = signif(kw$p.value, 3)
  ) %>%
  ungroup() %>%
  select(panel, dataset, chi_sq, df, p_value)

write.csv(kw_summary, 'data/sensitivity_table_KW_taxonomic_resolution.csv', row.names = FALSE)

cat('\n=== Kruskal-Wallis sensitivity summary (Fig. 3A-B, excl. N1 & S1) ===\n')
print(kw_summary)


## 1c — Figure: full dataset vs. N1/S1-excluded, panels A & B side by side --

make_kw_panel <- function(df, value_col, ylab, title) {
  kw <- kruskal.test(df[[value_col]] ~ df$Cluster)
  kw_label <- paste0('Kruskal-Wallis, p = ', signif(kw$p.value, 2))
  y_max <- max(df[[value_col]], na.rm = TRUE)

  ggboxplot(df, x = 'Cluster', y = value_col, fill = 'Cluster', width = 0.5, outlier.shape = NA) +
    geom_jitter(width = 0.08, size = 1.6, shape = 21, color = 'black', stroke = 0.3) +
    stat_compare_means(
      comparisons = list(c('C 1a', 'C 2')),
      method = 'wilcox.test',
      label = 'p.signif',
      tip.length = 0.02,
      size = 3.5,
      label.y = y_max * 1.12
    ) +
    annotate('text', x = 0.6, y = y_max * 1.28, label = kw_label, hjust = 0, size = 3) +
    scale_fill_manual(values = cluster_cols) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) +
    labs(x = '', y = ylab, title = title) +
    journal_theme
}

p_full_A <- make_kw_panel(datasets$full_A, 'OTU_count', 'n Calcareous ASVs', 'All samples (n = 14)')
p_full_B <- make_kw_panel(datasets$full_B, 'rel_reads_calcareus', '% Calcareous reads', 'All samples (n = 14)')
p_excl_A <- make_kw_panel(datasets$excl_A, 'OTU_count', 'n Calcareous ASVs', 'Excl. N1 & S1 (n = 12)')
p_excl_B <- make_kw_panel(datasets$excl_B, 'rel_reads_calcareus', '% Calcareous reads', 'Excl. N1 & S1 (n = 12)')

sensitivity_figure_taxonomic <- (p_full_A | p_full_B) / (p_excl_A | p_excl_B) +
  plot_annotation(
    title = 'Sensitivity of Fig. 3A-B to exclusion of the lowest taxonomic-resolution stations (N1, S1)',
    tag_levels = 'A',
    theme = theme(plot.title = element_text(size = 11, face = 'bold', hjust = 0.5))
  )

ggsave(
  filename = 'figures/sensitivity_figure_taxonomic_resolution.png',
  plot = sensitivity_figure_taxonomic,
  width = 200, height = 180, units = 'mm', dpi = 300, bg = 'white'
)


# ==============================================================================
# Part 2 — pH threshold-crossing robustness to inter-model (CMIP6) uncertainty
# ==============================================================================
#
# Bio-ORACLE v3.0 (Assis et al., 2024) projections are a 10-model CMIP6
# ensemble mean; the standard deviation across that ensemble (ph_sd) is also
# available per grid cell, quantifying inter-model disagreement. This
# section tests whether the pH threshold-crossing classification underlying
# Fig. 4 (decadal minimum pH within [7.02, 7.69], matching Script 04's
# empirical C1a/C2 cluster range) would survive that disagreement.
#
# Fully self-contained: downloads ph_min and ph_sd directly from Bio-ORACLE
# for the 2090-2100 decade, across SSP scenarios.

threshold <- 7.7
up_threshold_C1a <- 7.61
threshold_minC2 <- 7.69
low_threshold_C1a <- 7.02

BO_ph_layers <- biooracler::list_layers('ph') %>%
  filter(grepl('pH', title)) %>%
  filter(!grepl('surf', dataset_id)) %>%
  filter(!grepl('baseline', dataset_id)) %>%
  filter(grepl('depthmax', dataset_id))  # matches benthic_level used in Fig. 4

bo_ph_layers_id <- BO_ph_layers$dataset_id

extract_ph_min_sd <- function(dataset_id) {
  time <- c('2090-01-01T00:00:00Z', '2090-01-01T00:00:00Z')  # 2090-2100 decade only
  constraints <- list(time, c(30, 45), c(-10, 40))
  names(constraints) <- c('time', 'latitude', 'longitude')

  layers <- download_layers(dataset_id, c('ph_min', 'ph_sd'), constraints)
  cropped_raster <- crop(layers, ext(-10, 40, 30, 46))
  raster2df <- as.data.frame(cropped_raster, xy = TRUE)
  names(raster2df)[3:4] <- c('predicted_min', 'ph_sd')

  raster2df %>%
    mutate(Med_limits = case_when(
      x < 1 & y > 41 ~ 'drop_atl',
      x > 25 & y > 40 ~ 'drop_black',
      .default = 'keep'
    )) %>%
    filter(Med_limits == 'keep') %>%
    filter(x > -5.6) %>%
    mutate(scenario = dataset_id)
}

ph_results <- do.call(rbind, lapply(bo_ph_layers_id, extract_ph_min_sd)) %>%
  mutate(scenario_label = case_when(
    grepl('ssp585', scenario) ~ 'ssp5 8.5',
    grepl('ssp460', scenario) ~ 'ssp4 6.0',
    grepl('ssp370', scenario) ~ 'ssp3 7.0',
    grepl('ssp245', scenario) ~ 'ssp2 4.5',
    grepl('ssp119', scenario) ~ 'ssp1 1.9',
    .default = 'ssp1 2.6'
  ))

scenario_order <- c('ssp1 2.6', 'ssp2 4.5', 'ssp4 6.0', 'ssp3 7.0', 'ssp5 8.5')
ph_results$scenario_label <- factor(ph_results$scenario_label, levels = scenario_order)

# At-risk cells: predicted_min within [7.02, 7.69], matching Script 04's
# chronic/episodic classification -- not a plain <= 7.7 cutoff, to avoid
# extrapolating into pH conditions more extreme than the observed clusters.
at_risk <- ph_results %>%
  filter(!is.na(scenario_label)) %>%   # excludes ssp1-1.9, as in the main analysis
  filter(
    (predicted_min > up_threshold_C1a & predicted_min <= threshold_minC2) |
    (predicted_min <= up_threshold_C1a & predicted_min >= low_threshold_C1a)
  ) %>%
  mutate(
    margin = threshold - predicted_min,
    robust = margin > ph_sd
  )

robustness_summary <- at_risk %>%
  group_by(scenario_label) %>%
  summarise(
    n_cells = n(),
    median_margin = round(median(margin), 3),
    median_ph_sd = round(median(ph_sd, na.rm = TRUE), 3),
    pct_robust = round(100 * mean(robust, na.rm = TRUE), 1),
    .groups = 'drop'
  )

write.csv(robustness_summary, 'data/sensitivity_table_ph_model_uncertainty.csv', row.names = FALSE)

cat('\n=== pH threshold-crossing robustness to inter-model SD, 2090-2100 ===\n')
print(robustness_summary, n = Inf)
cat('\n"pct_robust" = % of at-risk cells per scenario whose margin below 7.7\n')
cat('exceeds the local inter-model standard deviation -- i.e., cells that would\n')
cat('very likely keep their "at risk" status even accounting for the full spread\n')
cat('of the 10-model Bio-ORACLE v3.0 ensemble.\n')

cat('\nSaved:\n',
    ' - data/sensitivity_table_taxonomic_resolution.csv\n',
    ' - figures/sensitivity_table_taxonomic_resolution.png\n',
    ' - data/sensitivity_table_KW_taxonomic_resolution.csv\n',
    ' - figures/sensitivity_figure_taxonomic_resolution.png\n',
    ' - data/sensitivity_table_ph_model_uncertainty.csv\n')
