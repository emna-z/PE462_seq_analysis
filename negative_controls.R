##################################################
## Project: PE462
## Script purpose: Negative controls exploration
## Date: 14-01-2023
## Author: E. Zeghal
##################################################



# init libraries ----------------------------------------------------------


library("tidyverse")
library("Hmisc")
library("phyloseq")
library("metagMisc")
library("mia")
library("miaViz")
library("microbiome")

# import data -------------------------------------------------------------


tax_pseq <- read.csv("../tax_tab_physeq.csv", row.names = 1) %>%
  as.matrix() %>%
  tax_table()

asv_pseq <- otu_table(as.matrix(read.csv("../asv_tab_physeq.csv", row.names = 1)), taxa_are_rows = T)

map <- read.csv("../map_from_tidy.csv",row.names = 1, na.strings = c("NA", ""), stringsAsFactors = T) |> 
  mutate(across(.cols = c("polymer_station","pol_photo_station"), 
                ~str_replace_all(.,c("C05"= "open_water", "C13" = "coast"))))

map_pseq <- sample_data(map)

pseq <- merge_phyloseq(asv_pseq, tax_pseq, map_pseq) 

pseq_neg <- pseq %>% 
  subset_samples(polymer %in% c("Neg_PCR"))



phyla <- plot_bar(pseq_neg, fill="Phylum") + 
  geom_bar(aes(color = Phylum, fill = Phylum), stat="identity", position="stack") +
  labs(x = "", y = "Relative Abundance\n") +
  theme(panel.background = element_blank())

plotly::ggplotly(phyla)













tse <- mia::makeTreeSEFromPhyloseq(pseq_neg)

tse_r <- transformCounts(tse,
                         assay_name = "counts",
                         method = "relabundance")

rowData(tse)$Phylum %>% table()%>% sort() %>% kable("pipe")

mapTaxonomy(tse_r)

head(getPrevalence(tse, detection = 1/100, sort = TRUE, as_relative = TRUE))
head(getPrevalence(tse, detection = 1, sort = TRUE, assay_name = "counts", as_relative = FALSE))


tse_phylum <- agglomerateByRank(tse_r, rank ="Phylum", onRankOnly=TRUE)
top_taxa <- getTopTaxa(tse_phylum,top = 12, assay_name = "relabundance")

# Renaming the "Phylum" rank to keep only top taxa and the rest to "Other"
phylum_renamed <- lapply(rowData(tse)$Phylum,
                         function(x){if (x %in% top_taxa) {x} else {"Other"}})
rowData(tse)$Phylum <- as.character(phylum_renamed)

plotAbundance(tse, rank = "Phylum",
              order_rank_by="abund")
