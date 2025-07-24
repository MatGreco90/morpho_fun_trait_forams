# morpho_fun_trait_forams
This repo contains the R scripts used to run the analysis and figures presented in the manuscript "Shaping shells: predicting the shifts in Mediterranean benthic calcifiers' coastal assemblages under future ocean acidification scenarios"

## Scripts
**`01_data_preprocessing_phyloseq.R`**  
_Prepares and formats ASV community data using the `phyloseq` package. Steps include importing raw OTU/ASV tables, metadata integration to generate `phyloseq` object, and performing rarefaction. The script also calculates alpha diversity metrics and applies CLR (centered log-ratio) transformation._

**`02_alpha_diversity_pH_correlation.R`**  
_Assesses the Spearman correlation between alpha diversity indices and pH levels._

**`03_betadiv_analyses.R`**  
_Performs beta diversity analyses, including PCoA, hierarchical clustering, and BIOENV._

**`04_extract_min_pH_biooracle_projections.R`**  
_Extracts minimum projected pH values from Bio-ORACLE environmental datasets under future climate scenarios._

**`05_visualization.R`**  
_Generates final figures for the manuscript._
