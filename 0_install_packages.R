
# tidyverse ---------------------------------------------------------------
install.packages("tidyverse")

# Devtools ----------------------------------------------------------------

install.packages("devtools")

# phyloseq ----------------------------------------------------------------
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("phyloseq")

# miaverse ----------------------------------------------------------------
BiocManager::install("mia")
BiocManager::install("microbiome")

# Hmisc -------------------------------------------------------------------

install.packages("Hmisc")


# GT suite ----------------------------------------------------------------

install.packages("gt")
install.packages("gtExtras")
install.packages("gtsummary")

# scales ------------------------------------------------------------------

install.packages("scales")

# plots and palettes ------------------------------------------------------

install.packages("RColorBrewer")
install.packages("ggthemes")
install.packages("Polychrome")
install.packages("ggpubr")
install.packages("ggpubr")
install.packages("ggiraph")
install.packages("plotly")
devtools::install_github("david-barnett/microViz")
