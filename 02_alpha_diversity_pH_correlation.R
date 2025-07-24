# Load required libraries
library(Hmisc) # For rcorr(), which computes correlation matrices with p-values
library(tidyverse) # For data manipulation

# Load environmental data and extract site names and pH values
env_data <- read.csv('metadata.csv') %>% 
  mutate(site=tolower(site)) %>% 
  select(site,pH)

# Load alpha diversity table from Script 01 and merge with pH data
alpha_tab <- readRDS('data/alpha_div_rarefied.RDS') %>% 
  rownames_to_column('site') %>% 
  left_join(env_data, by='site') %>% 
  column_to_rownames('site')

# Define a function to compute Spearman correlation between each diversity index and pH
spear_cor<- function(x){
  my_cor<-rcorr(x,alpha_tab$pH,'spearman')
  my_r<-my_cor$r[2]
  my_p<-my_cor$P[2]
  out<-paste(my_r,my_p, sep = '_')
  return(out)
}

# Apply correlation function to all diversity indices in the table
ph_corr<-alpha_tab %>%
  summarise(across(1:22,spear_cor)) %>% 
  pivot_longer(1:ncol(.), values_to = "val",names_to='diversity_index', names_repair = "unique") %>% 
  separate(val ,into=c('rho','p_val'),sep='_') %>% 
  mutate(across(2:3,as.numeric)) 


# Save the resulting correlation table 
write.csv(ph_corr,'data/Table_S1.csv', row.names = FALSE)




