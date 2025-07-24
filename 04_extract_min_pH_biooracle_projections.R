# Load necessary libraries
library(tidyverse)   # Data manipulation and reshaping
library(terra)       # Raster data handling and spatial operations
library(sf)          # Simple features for spatial data (not explicitly used here but loaded)
library(biooracler)  # Interface to Bio-ORACLE marine environmental data


# List available pH-related Bio-ORACLE layers, filtering to exclude surface layers and baseline datasets
BO_ph_layers <- biooracler::list_layers("ph") %>% 
  filter(grepl('pH',title)) %>% 
  filter(!grepl('surf',dataset_id)) %>% 
  filter(!grepl('baseline',dataset_id))


# Extract the dataset IDs for the filtered layers
bo_ph_layers_id <- BO_ph_layers$dataset_id

# Download pH min ---------------------------------------------------------

# Define a function to download and process minimum pH projections for a given dataset
extract_ph_projections_points <- function(x) {
  
  dataset_id <- x
  
  # Define constraints for downloading data: time range and spatial bounding box
  time = c('2020-01-01T00:00:00Z', '2090-01-01T00:00:00Z')
  
  latitude = c(30, 45)
  longitude = c(-10, 40)
  
  constraints = list(time, latitude, longitude)
  names(constraints) = c("time", "latitude", "longitude")
  
  variables = "ph_min" # Variable to download: minimum pH projections
  
  # Download raster layers for the dataset with given constraints
  
  layers <- download_layers(dataset_id, variables, constraints)
  
  # Define bounding box (extent) for cropping to Mediterranean Sea
  
  bbox_mediterranean <- ext(-10, 40, 30, 46)  # Rough bbox for the Mediterranean
  
  # Crop raster layers to Mediterranean bounding box
  cropped_raster <- crop(layers, bbox_mediterranean)
  
  # Convert raster data to data frame with coordinates 
  raster2df <- as.data.frame(cropped_raster, xy = TRUE)
  
  # Define pH thresholds for categorizing projected pH values
  
  threshold_minC2 <- 7.69
  up_threshold_C1a <- 7.61
  low_threshold_C1a <- 7.02
  
  
  # Reshape and categorize data based on pH thresholds and spatial filters
  
  df_threshold <- raster2df %>% 
    pivot_longer(3:10, names_to = 'var', values_to = 'predicted') %>% 
    mutate(ph_vals = case_when(
      predicted > up_threshold_C1a & predicted <= threshold_minC2 ~ 'lower C2',  # new condition for 7.61 < x < 7.69
      predicted <= up_threshold_C1a & predicted >= low_threshold_C1a ~ 'upper C1a quartile',  # upper quartile
      .default = 'drop')) %>% 
    filter(!ph_vals=='drop') %>% 
    mutate(year=case_when(grepl('_1',var)~2020,
                          grepl('_2',var)~2030,
                          grepl('_3',var)~2040,
                          grepl('_4',var)~2050,
                          grepl('_5',var)~2060,
                          grepl('_6',var)~2070,
                          grepl('_7',var)~2080,
                          .default=2090)) %>% 
    # Filter out points outside Mediterranean proper based on coordinate rules
    mutate(Med_limits=case_when(x<1 & y>41~'drop_atl', # Drop Atlantic region north of 41 latitude west of 1 longitude
                                x>25 & y>40~'drop_black', # Drop Black Sea region east of 25 longitude north of 40 latitude
                                .default = 'keep')) %>%  
    filter(Med_limits=='keep') %>% # Keep only Mediterranean points
    filter(x>(-5.6)) # Additional longitude filter for Mediterranean
  
  # Add a scenario identifier column for downstream analysis
  res_df <- df_threshold %>%
    mutate(scenario=dataset_id)
  
  return(res_df)
}

# Apply the extraction function to all pH-related datasets to get a list of data frames
points_results_list <- lapply(bo_ph_layers_id, extract_ph_projections_points)


# Combine all the data frames in the list into a single data frame
points_combined_results <- do.call(rbind, points_results_list) %>% 
  mutate(scenario_label=case_when(grepl('ssp585',scenario)~'ssp5 8.5',
                                  grepl('ssp460',scenario)~'ssp4 6.0',
                                  grepl('ssp370',scenario)~'ssp3 7.0',
                                  grepl('ssp245',scenario)~'ssp2 4.5',
                                  grepl('ssp119', scenario)~'ssp1 1.9',
                                  .default = 'ssp1 2.6')) %>% 
  mutate(benthic_level=case_when(grepl('depthmax',scenario)~'depth max',
                                 grepl('depthmin',scenario)~'depth min',
                                 grepl('depthmean',scenario)~'depth mean')) %>% 
  mutate(decade=case_when(year==2020~'2020-2030',
                          year==2030~'2030-2040',
                          year==2040~'2040-2050',
                          year==2050~'2050-2060',
                          year==2060~'2060-2070',
                          year==2070~'2070-2080',
                          year==2080~'2080-2090',
                          .default = '2090-2100'))

# Save the final combined dataset as an RDS file 
saveRDS(points_combined_results,'data/points_biooracle_ph_min.RDS')
