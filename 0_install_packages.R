
# tidyverse ---------------------------------------------------------------
install.packages("tidyverse")

# phyloseq ----------------------------------------------------------------
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("phyloseq")

# mia ---------------------------------------------------------------------
BiocManager::install("mia")

# Hmisc -------------------------------------------------------------------

install.packages("Hmisc")

