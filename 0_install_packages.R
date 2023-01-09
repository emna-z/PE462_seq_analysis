
# tidyverse ---------------------------------------------------------------
install.packages("tidyverse")

# phyloseq ----------------------------------------------------------------
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("phyloseq")

# miaverse ---------------------------------------------------------------------
BiocManager::install("mia")
BiocManager::install("microbiome")

# Hmisc -------------------------------------------------------------------

install.packages("Hmisc")



