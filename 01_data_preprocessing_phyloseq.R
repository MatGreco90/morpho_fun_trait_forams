# Load required libraries
library(tidyverse)     # For data manipulation 
library(microbiome)    # For microbiome-specific transformations (e.g., CLR, alpha diversity)
library(phyloseq)      # For managing microbiome data structures and performing analysis

# Import ASV table
asv_table <- read.csv('asv_table.csv')

# Extract taxonomic information and shell type from the ASV table
tax_info <- asv_table %>% 
  select(asv_id, taxon,shell_type) %>% 
  separate(taxon, into=c('class','order','family','genus','species','variant'), sep=';') %>% 
  column_to_rownames('asv_id')

# Convert taxonomic table to matrix format (required for phyloseq::tax_table)
tax_info <- as.matrix(tax_info)

# Convert abundance table to matrix format 
asv_matrix <- asv_table %>% 
  select(1:15) %>%   
  column_to_rownames('asv_id')

asv_matrix <- as.matrix(asv_matrix)

# Read environmental metadata and standardize site names
env_data <- read.csv('data/metadata.csv') %>% 
  mutate(site=tolower(site)) %>% 
  column_to_rownames('site')

# Create phyloseq components from the data matrices
OTU = otu_table(asv_matrix, taxa_are_rows = TRUE)
TAX = tax_table(tax_info)
samples = sample_data(env_data)

# Combine all components into a single phyloseq object

phylo_obj <- phyloseq(OTU, TAX, samples)

# Save raw (unprocessed) phyloseq object to file

saveRDS(phylo_obj,'data/phylo_obj_raw_dat.RDS')

# Perform rarefaction to even sampling depth (minimum total reads per sample)
phylo_obj_rarefied = rarefy_even_depth(phylo_obj,
                                       rngseed=23, # Set seed for reproducibility
                                       sample.size=min(sample_sums(phylo_obj)), # Use smallest sample size as rarefaction depth
                                       replace=F) # Do not allow replacement of counts during rarefaction

# Save rarefied phyloseq object
saveRDS(phylo_obj_rarefied,'data/phylo_obj_rarefied_dat.RDS')

# Calculate all available alpha diversity metrics on rarefied data
alpha_tab <-microbiome::alpha(phylo_obj_rarefied, index = "all")

# Save alpha diversity table
saveRDS(alpha_tab,'data/alpha_div_rarefied.RDS')

# Apply Centered Log-Ratio (CLR) transformation to ASV counts
phylo_obj_clr <-microbiome::transform(phylo_obj, 'clr')

# Save CLR-transformed phyloseq object
saveRDS(phylo_obj_clr,'data/phylo_obj_clr_transformed.RDS')

