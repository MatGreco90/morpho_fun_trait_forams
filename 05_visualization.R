
# Figure 1 ----------------------------------------------------------------

# Load necessary libraries
library(phyloseq) #  For handling and analyzing microbiome data
library(tidyverse) # For data manipulation 
library(ggthemes) # For visualization


#Load raw phyloseq obj (not-transformed data ) from Script 01
raw_data <- readRDS('data/phylo_obj_raw_dat.RDS') 

#Load env data
env_data <- read.csv('metadata.csv') %>% 
  mutate(site=tolower(site)) %>% 
  select(site,pH)

#Load alpha div table from Script 01
alpha_tab <- readRDS('data/alpha_div_rarefied.RDS') %>% 
  rownames_to_column('site') %>% 
  left_join(env_data, by='site')


alpha_tab$transect <- ifelse(grepl("^n", alpha_tab$site), "North", "South")


#Derive % ASVs per sample fr each class 

rel_asv_table <- psmelt(raw_data) %>%  
  select(OTU, Sample,class, Abundance) %>% 
  group_by(Sample,OTU) %>% 
  mutate(tot_reads_OTU_sample=sum(Abundance)) %>% 
  ungroup() %>% 
  filter(tot_reads_OTU_sample>0) %>%
  group_by(Sample, class) %>%
  summarise(OTU_count = n_distinct(OTU)) %>% 
  ungroup() %>% 
  group_by(Sample) %>% 
  mutate(tot_OTU_count=sum(OTU_count)) %>% 
  ungroup() %>% 
  mutate(rel_OTU_count=OTU_count/tot_OTU_count)

#Generate plot (panel B)

panel_B <- ggplot(rel_asv_table, aes(x = factor(Sample,
                                                levels=c('n1','n2','n3',
                                                         'n4','s1','s2','s3',
                                                         's4','s5','s6','s7',
                                                         's8','s9','s10')), y = rel_OTU_count, fill = class)) +
  geom_bar(stat = "identity", position = "stack", col='black') +
  scale_fill_tableau(palette = 'Tableau 10')+
  theme_classic() +
  labs(x = "", y = "% ASVs", fill = "Class") +
  theme(axis.title=element_text(face='bold'),
        axis.text.x = element_text(size=13,face = 'bold'),
        axis.text.y=element_text(size=10),
        legend.text = element_text(size=11),
        legend.title =element_blank(),
        legend.position = 'none')

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

panel_C <-ggplot(rel_reads_table, aes(x = factor(Sample,
                                                 levels=c('n1','n2','n3',
                                                          'n4','s1','s2','s3',
                                                          's4','s5','s6','s7',
                                                          's8','s9','s10')), y = rel_read_abundance, fill = class)) +
  geom_bar(stat = "identity", position = "stack", col='black') +
  scale_fill_tableau(palette = 'Tableau 10')+
  theme_classic() +
  labs(x = "", y = "% reads", fill = "Class") +
  theme(axis.title=element_text(face='bold'),
        axis.text.x = element_text(size=13,face = 'bold'),
        axis.text.y=element_text(size=10),
        legend.text = element_text(size=11),
        legend.title =element_blank(),
        legend.position = 'bottom')


#generate plot panel D

panel_D <-ggplot(alpha_tab, aes(x = factor(site,
                                           levels = c('n1','n2','n3',
                                                      'n4','s1','s2','s3',
                                                      's4','s5','s6','s7',
                                                      's8','s9','s10')),
                                y = diversity_shannon)) +
  geom_line(aes(group = transect), color = "gray30", size = 0.6) +
  geom_point(aes(fill = pH), shape = 21, size = 5) +
  scale_fill_viridis_c() +
  theme_minimal()+
  labs(x = "", y = "Shannon", fill = "Class") +
  theme(axis.title=element_text(face='bold'),
        axis.text.x = element_text(size=13,face = 'bold'),
        axis.text.y=element_text(size=10),
        legend.text = element_text(size=11),
        legend.title =element_blank(),
        legend.position = 'bottom')

#combine plots
side_plot <- panel_B/panel_C/panel_D


# Figure 2 ----------------------------------------------------------------

# Load necessary libraries
library(vegan)
library(ggtree)
library(tidyverse)
library(patchwork)
library(ggdist)


samples_meta <- read.csv('metadata.csv') %>% 
  mutate(site=tolower(site)) %>%  
  dplyr::select(1:4)

# Manually define sample groupings from cluster 
clust_1a <- c('s2','s4','s3','s5','s1')
clust_1b <- c('s6','s8','n1','s9','s10')
clust_2 <- c('n3','s7','n4','n2')


samples_meta_cluster <- samples_meta %>% 
  mutate(Cluster=ifelse(site %in% clust_1a,'C 1a',
                        ifelse(site %in% clust_1b,'C 1b','C 2'))) 

#panel A

#load cluster res (Script 03)
cluster_res <- readRDS('data/hclust_res.RDS')


p <- ggtree(cluster_res,size=0.6) %>% 
  rotate(15) %>% 
  rotate(16)+geom_tiplab(vjust = -1.4, hjust =- 0.1,
                         size=3)

#add cluster metadata to color-code samples
final_cluster_plot<-p %<+% samples_meta_cluster + 
  layout_dendrogram() +
  geom_tippoint(size=3, aes(fill=Cluster), shape=21)+
  theme_dendrogram()+
  theme(legend.text = element_text(size=12),
        legend.title = element_text(size=14, face='bold'))


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
spp.scrs <- as.data.frame(vegan::scores(env_fit_res, display = "vectors")) * arrow_factor
spp.scrs <- cbind(spp.scrs, site = rownames(spp.scrs), Pvalues = env_fit_res$vectors$pvals, R_squared = env_fit_res$vectors$r)

spp.scrs <- subset(spp.scrs, Pvalues < 0.05)


# generate convex hulls
pcoa_coords_plot <- pcoa_coords %>% 
  mutate(cluster=ifelse(site %in% clust_1a,'C_1a',
                        ifelse(site %in% clust_1b,'C_1b','C_2'))) 



pcoa_coords_chull <-pcoa_coords_plot %>%  
  group_by(cluster) %>% 
  slice(chull(V1, V2))


#plot pcoa 
pcoa_plot <- ggplot(pcoa_coords_plot) +
  geom_polygon(data = pcoa_coords_chull, alpha=0.5,
               aes(x = V1, y = V2, fill = cluster), inherit.aes = F)+
  geom_point(mapping = aes(x = V1, y = V2, fill=cluster), shape=21,col='black',size=3) +
  coord_fixed() + ## need aspect ratio of 1!
  geom_segment(data = spp.scrs,
               aes(x = 0, xend = Dim1 , y = 0, yend = Dim2 ),
               arrow = arrow(length = unit(0.25, "cm")), colour = "black")+
  geom_text(data = spp.scrs, aes(x = Dim1+5 , y = Dim2+0.05, label = site),
            size = 6)+
  theme_classic()+
  xlab('Dim 1')+
  ylab('Dim 2')+
  guides(fill='none')+
  annotate(geom="text", x=-25, y=-60, label="r2 = 0.46", color="black")+
  theme(legend.position = 'none')


#Panel C

#plot cluster along pH gradient 

panel_c <- ggplot(samples_meta_cluster,aes(fill = Cluster, color = Cluster, x = pH)) +
  stat_slab(alpha = .3) +
  stat_pointinterval(position = position_dodgejust(width = .2), justification = 0.1) +
  scale_y_continuous(breaks = NULL)+
  geom_vline(xintercept = 7.7,linetype = 3)+
  ylab('')+
  theme_bw()+
  theme(legend.position = 'none')


(final_cluster_plot+pcoa_plot)/panel_c+
  plot_annotation(tag_levels = 'A')+
  plot_layout(guides = "collect") 


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
samples_meta <- read.csv('metadata.csv') %>% 
  mutate(site=tolower(site)) %>%  
  dplyr::select(1:4)

# Manually define sample groupings from cluster 
clust_1a <- c('s2','s4','s3','s5','s1')
clust_1b <- c('s6','s8','n1','s9','s10')
clust_2 <- c('n3','s7','n4','n2')

samples_meta_cluster <- samples_meta %>% 
  mutate(Cluster=ifelse(site %in% clust_1a,'C 1a',
                        ifelse(site %in% clust_1b,'C 1b','C 2'))) %>% 
  rename('Sample'='site')


#calc n calcareous asvs per cluster

n_asvs_cluster_table <- psmelt(raw_data) %>% 
  select(OTU, Sample,class, Abundance,shell_type) %>%  
  group_by(Sample,OTU) %>% 
  mutate(tot_reads_OTU_sample=sum(Abundance)) %>% 
  ungroup() %>% 
  filter(tot_reads_OTU_sample>0) %>%
  select(OTU, Sample,,shell_type) %>% 
  left_join(samples_meta_cluster, by='Sample') %>%  
  filter(shell_type=='calcareous') %>% 
  group_by(Sample) %>% 
  mutate(OTU_count = n_distinct(OTU)) %>% 
  ungroup() %>% 
  select(Sample, OTU_count, Cluster) %>% 
  distinct()

n_asvs_cluster_table$Cluster<- factor(n_asvs_cluster_table$Cluster, levels = c('C 1a', 
                                                                               'C 1b',
                                                                               'C 2'))

#perform paired test with fdr correction
ggbetweenstats(n_asvs_cluster_table, Cluster, OTU_count,
               plot.type = "box",
               type = "np",
               pairwise.annotation = "p.value",
               pairwise.display = "significant", ## display only significant pairwise comparisons
               p.adjust.method = "fdr", ## adjust p-values for multiple tests using this method
               ylab = 'n Calcareous ASVs',
               results.subtitle = TRUE)


#select only significant comparison

my_comparisons <-list( c("C 2", "C 1a"))

#visualize res

asv_count <- ggboxplot(n_asvs_cluster_table, x = "Cluster", y = "OTU_count",
                       fill = "Cluster", palette = "viridis",
                       width = 0.5)+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")+ # Add pairwise comparisons p-value
  # stat_compare_means(label.y = 90)+
  stat_compare_means(label.y = 80, , size=3.5)+
  # geom_jitter()+
  ylab('n Calcareous ASVs')+
  xlab('')+
  theme_bw()+
  theme(axis.title=element_text(face='bold'),
        axis.text.x = element_text(size=13,face = 'bold'),
        axis.text.y=element_text(size=10),
        legend.position = 'none')


#calc % calcareous asvs per cluster

relab_cluster_table <- psmelt(raw_data) %>% 
  select(OTU, Sample, Abundance,shell_type) %>%
  left_join(samples_meta_cluster, by='Sample') %>%
  filter(!shell_type=='not_known') %>%
  group_by(Sample) %>%
  mutate(tot_reads_by_sample=sum(Abundance)) %>%  
  ungroup() %>% 
  group_by(Sample,shell_type) %>% 
  mutate(reads_by_sample_shell_type=sum(Abundance)) %>% 
  ungroup() %>% 
  filter(shell_type=='calcareous') %>%
  mutate(rel_reads_calcareus=(reads_by_sample_shell_type/tot_reads_by_sample)*100) %>% 
  select(Sample, rel_reads_calcareus, Cluster) %>% 
  distinct()


relab_cluster_table$Cluster<- factor(relab_cluster_table$Cluster, levels = c('C 1a', 
                                                                             'C 1b',
                                                                             'C 2'))

#perform paired test with fdr correction
ggbetweenstats(relab_cluster_table, Cluster, rel_reads_calcareus,
               plot.type = "box",
               type = "np",
               pairwise.display = "significant", ## display only significant pairwise comparisons
               p.adjust.method = "fdr", ## adjust p-values for multiple tests using this method
               ylab = '% Calcareous reads')

#select only significant comparison
my_comparisons <-list( c("C 2", "C 1a"))


#visualize res
reads_plot <- ggboxplot(relab_cluster_table, x = "Cluster", y = "rel_reads_calcareus",
                        fill = "Cluster", palette = "viridis",
                        width = 0.5)+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")+ # Add pairwise comparisons p-value
  # stat_compare_means(label.y = 90)+
  stat_compare_means(label.y = 60, size=3.5)+
  ylab('% Calcareous reads')+
  xlab('')+
  theme_bw()+
  theme(axis.title=element_text(face='bold'),
        axis.text.x = element_text(size=13,face = 'bold'),
        axis.text.y=element_text(size=10),
        legend.position = 'none')



# Morphological data 


#load morph data 
morph_data <- read.csv('morph_counts.csv', header = TRUE) %>% 
  clean_names() %>% #clear column names 
  mutate(site=tolower(site)) %>% #make sample name consistent
  mutate(calcareous_tests=hy+porcel) #sum  %procellanaceous and hyaline 

#add cluster info

morph_data_cluster <- morph_data %>% 
  mutate(cluster=ifelse(site %in% clust_1a,'C 1a',
                        ifelse(site %in% clust_1b,'C 1b','C 2'))) 


morph_data_cluster$cluster<- factor(morph_data_cluster$cluster, levels = c('C 1a', 
                                                                           'C 1b',
                                                                           'C 2'))

#perform paired test with fdr correction
ggbetweenstats(morph_data_cluster, cluster, calcareous_tests,
               plot.type = "box",
               type = "np",
               pairwise.display = "significant", ## display only significant pairwise comparisons
               p.adjust.method = "fdr", ## adjust p-values for multiple tests using this method
               ylab = 'calcareous_tests')


#select only significant comparison
my_comparisons <-list( c("C 2", "C 1a"))

#visualize res
rel_calc_specimens <- ggboxplot(morph_data_cluster, x = "cluster", y = "calcareous_tests",
                                fill = "cluster",width = 0.5)+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")+ # Add pairwise comparisons p-value
  # stat_compare_means(label.y = 90)+
  stat_compare_means(label.y = 115, size=3.5)+
  ylab('% Calcareous specimens')+
  xlab('')+
  theme_bw()+
  theme(axis.title=element_text(face='bold'),
        axis.text.x = element_text(size=13,face = 'bold'),
        axis.text.y=element_text(size=10),
        legend.position = 'none')


#perform paired test with fdr correction  (Fragmentation index)

ggbetweenstats(morph_data_cluster, cluster, fi,plot.type = "box",
               type = "np",
               pairwise.display = "significant", ## display only significant pairwise comparisons
               p.adjust.method = "fdr", ## adjust p-values for multiple tests using this method
               ylab = 'FI')

#select only significant comparison
my_comparisons <-list( c("C 2", "C 1a"))

#visualize res
fi_plot <- ggboxplot(morph_data_cluster, x = "cluster", y = "fi",
                     fill = "cluster",width = 0.5)+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")+ # Add pairwise comparisons p-value
  # stat_compare_means(label.y = 90)+
  stat_compare_means(label.y = 90, size=3.5)+
  ylab('FI')+
  xlab('')+
  theme_bw()+
  theme(axis.title=element_text(face='bold'),
        axis.text.x = element_text(size=13,face = 'bold'),
        axis.text.y=element_text(size=10),
        legend.position = 'none')


#perform paired test with fdr correction  (% Living Specimens)

ggbetweenstats(morph_data_cluster, cluster, living,
               plot.type = "box",
               type = "np",
               pairwise.display = "significant", ## display only significant pairwise comparisons
               p.adjust.method = "fdr", ## adjust p-values for multiple tests using this method
               ylab = 'Living')

#select only significant comparison
my_comparisons <-list( c("C 2", "C 1a"))


#visualize res
living_plot <- ggboxplot(morph_data_cluster, x = "cluster", y = "living",
                         fill = "cluster",width = 0.5)+
  # stat_compare_means()+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")+ # Add pairwise comparisons p-value
  # stat_compare_means(label.y = 90)+
  stat_compare_means(label.y = 110, size=3.5)+
  ylab('% Living specimens')+
  xlab('')+
  theme_bw()+
  theme(axis.title=element_text(face='bold'),
        axis.text.x = element_text(size=13,face = 'bold'),
        axis.text.y=element_text(size=10),
        legend.position = 'none')


#perform paired test with fdr correction  (% Foraminiferal Density)

ggbetweenstats(morph_data_cluster, cluster, fd_d_l,
               plot.type = "box",
               type = "np",
               pairwise.display = "significant", ## display only significant pairwise comparisons
               p.adjust.method = "fdr", ## adjust p-values for multiple tests using this method
               ylab = 'FD(D+L)')

#select only significant comparison
my_comparisons <-list( c("C 2", "C 1a"))


#visualize res
fd_dl_plot <- ggboxplot(morph_data_cluster, x = "cluster", y = "fd_d_l",
                        fill = "cluster",width = 0.5)+
  # stat_compare_means()+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")+ # Add pairwise comparisons p-value
  # stat_compare_means(label.y = 90)+
  stat_compare_means(label.y = 160,, size=3.5)+
  ylab('FD(D+L)')+
  xlab('')+
  theme_bw()+
  theme(axis.title=element_text(face='bold'),
        axis.text.x = element_text(size=13,face = 'bold'),
        axis.text.y=element_text(size=10),
        legend.position = 'none')

#combine plots
asv_count+reads_plot+ rel_calc_specimens +fi_plot+living_plot+fd_dl_plot+
  plot_annotation(tag_levels = 'A')



# Figure 4 ----------------------------------------------------------------

# Load necessary libraries
library("rnaturalearth")
library("rnaturalearthdata")
library(tidyverse)


#load data from Script 04
points_combined_results_585 <- readRDS('data/points_biooracle_ph_min.RDS') %>% 
  filter(scenario_label=='ssp5 8.5') %>% #select scenario
  filter(benthic_level=='depth max') %>% #select benthic layer
  mutate(ph_vals=case_when(ph_vals=='lower C2'~'below 7.7',
                           .default = 'C1a range'))

world <- ne_countries(scale = "large", returnclass = "sf")


ggplot(data = world) +
  geom_sf(fill= 'gray') +
  coord_sf(xlim = c(-10, 40), ylim = c(30, 46), expand = FALSE)+
  geom_point(data=points_combined_results_585,aes(x,y,col=ph_vals),
             shape=17,size=1, alpha=0.8)+
  facet_wrap(~decade, ncol=2)+
  scale_color_manual(values=c('indianred','darkorange'))+
  theme_classic()+
  labs(col='pH level')+
  xlab('')+
  ylab('')+
  ggtitle('ssp5 8.5')


# Figure 5 ----------------------------------------------------------------


#load data from Script 04
points_combined_results_2090_2100 <- readRDS('data/points_biooracle_ph.RDS') %>%  
  filter(benthic_level=='depth max') %>% #select benthi layer
  filter(decade=='2090-2100')%>% #select decade
  mutate(ph_vals=case_when(ph_vals=='lower C2'~'below 7.7',
                           .default = 'C1a range'))

#plot decade 2090-2100 for each scenario
ggplot(data = world) +
  geom_sf(fill= 'gray') +
  coord_sf(xlim = c(-10, 40), ylim = c(30, 46), expand = FALSE)+
  geom_point(data=points_combined_results_2090_2100,aes(x,y,col=ph_vals),
             shape=17,size=1, alpha=0.8)+
  facet_wrap(~scenario_label, ncol=2)+
  scale_color_manual(values=c('indianred','darkorange'))+
  theme_classic()+
  labs(col='pH level')+
  xlab('')+
  ylab('')+
  ggtitle('Decade: 2090-2100')

