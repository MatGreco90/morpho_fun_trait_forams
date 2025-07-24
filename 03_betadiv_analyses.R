# Load required libraries
library(phyloseq)   # For handling and analyzing microbiome data
library(vegan)      # For ecological analyses, e.g., distance calculations, ordination
library(tidyverse)  # For data manipulation 

# Load CLR (Centered Log-Ratio) transformed phyloseq object (Script 01)

clr_data <- readRDS('data/phylo_obj_clr_transformed.RDS')

# Convert phyloseq object into a long data frame
phylo_2df <- psmelt(clr_data)

# Create ASV (Amplicon Sequence Variant) abundance table: samples as rows, OTUs as columns
clr_asv_table <- phylo_2df %>% 
  select(OTU, Sample, Abundance) %>% 
  pivot_wider(names_from = OTU, values_from = Abundance) %>% 
  column_to_rownames('Sample')

# Compute Euclidean distance between samples based on CLR-transformed abundances

dist_samples <- vegdist(clr_asv_table, method = 'euclidean') 

# Perform hierarchical clustering on the distance matrix using the complete linkage method
clust_samples <-  hclust(dist_samples,method = 'complete')

saveRDS(clust_samples,'data/hclust_res.RDS')

plot(clust_samples)


# Manually define sample groupings from dendrogram for downstream analysis (assumed clusters)
clust_1a <- c('s2','s4','s3','s5','s1')
clust_1b <- c('s6','s8','n1','s9','s10')
clust_2 <- c('n3','s7','n4','n2')


# Extract and deduplicate environmental metadata for each sample
env_data <-phylo_2df %>% 
  select(Sample, pH, O2, depth) %>%
  distinct() 

# Merge abundance and environmental data into a single table

full_table_rda <- clr_asv_table %>% 
  rownames_to_column('Sample') %>% 
  left_join(env_data, by='Sample') %>% 
  column_to_rownames('Sample') 

# Extract environmental matrix
env_mat <- full_table_rda[,1200:1202]

# Perform Principal Coordinates Analysis (PCoA) using the same Euclidean distance matrix
pcoa_result <- cmdscale(dist_samples, eig = TRUE)

saveRDS(pcoa_result,'data/PcoA_res.RDS')

# Visualize the PCoA result (simple ordination plot)
ordiplot(pcoa_result)

# Fit environmental variables to the ordination to assess correlations (999 permutations)
set.seed(23) # Set seed for reproducibility
envfit_result <- envfit(pcoa_result ~ ., data = env_mat,permu=999)
# View significance and direction of environmental vectors in ordination
envfit_result 

saveRDS(envfit_result,'data/envfit_res.RDS')
