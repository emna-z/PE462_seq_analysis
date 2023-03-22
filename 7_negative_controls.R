
# phyloseq route ----------------------------------------------------------


tax_pseq <- read.csv("../tax_tab_physeq.csv", row.names = 1) %>%
  as.matrix() %>%
  tax_table()

asv_pseq <- otu_table(as.matrix(read.csv("../asv_tab_physeq.csv", row.names = 1)), taxa_are_rows = T)


map <- read.csv("../map_from_tidy.csv",row.names = 1, na.strings = c("NA", ""), stringsAsFactors = T) %>% 
  mutate(across(.cols = c("polymer_station","pol_photo_station"), 
                ~str_replace_all(.,c("C05"= "open NS water", "C13" = "coastal NS water", "_"= " "))))

map_pseq <- sample_data(map)

pseq <- merge_phyloseq(asv_pseq, tax_pseq, map_pseq) %>% 
  subset_samples(station %nin% c("mock_DNA"))

library(MicEco)
ps_venn(pseq, "station", weight = T)

# 
# library(devtools)
# install_github("Russel88/MicEco")

pseq_neg <- pseq %>% 
  subset_samples(polymer %in% c("Neg_PCR")) %>% 
  transform_sample_counts(function(x) x / sum(x) ) %>% 
  
  
  plot_bar(pseq_neg, fill = "Phylum")




# tidy route --------------------------------------------------------------


t <- read_csv("../tidyPE462_rel_abund_calc.csv")
