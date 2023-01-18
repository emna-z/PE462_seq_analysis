##################################################
## Project: PE462
## Script purpose: install the packages necessary
## for the analysis workflow
## Date: 17-01-2023
## Author: E. Zeghal
##################################################

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

# Tables ------------------------------------------------------------------

install.packages("kableExtra")
install.packages("huxtable")
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
install.packages("ggtext")
install.packages("ggnewscale")
install.packages("ggh4x")
