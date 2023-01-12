##################################################
## Project: PE462
## Script purpose: Importing raw data & cleaning
## Date: 04-01-2023
## Author: E. Zeghal
##################################################


# load libraries ----------------------------------------------------------

library("magrittr")
library("tidyverse")
library("phyloseq")
library("mia")
library("microbiome")
library("Hmisc")

# raw data import ---------------------------------------------------------

tax1 <- read.delim("../both_lanes/asv/taxonomy_dada2/representative_seq_set_tax_assignments.txt",
                             row.names = 1, na.strings = "NA")
tax2 <- read.delim("../unpaired_PR2/asv/taxonomy_dada2/representative_seq_set_tax_assignments.txt",
                   row.names = 1, na.strings = "NA")

rownames(tax2) <- paste0("asv.",(length(rownames(tax1))+1):(length(rownames(tax1))+length(rownames(tax2))))

# for ease of handling we'll change the rank Division to Phylum for the Eukaryotes table
# and add a Supergroup column to the prokaryotes column by duplicating the kingdom column content
ranks <- c ("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
colnames(tax1) <-  ranks
tax1 <- tax1 %>% mutate(Supergroup = Kingdom, .after = Kingdom)
colnames(tax2)[3] <- "Phylum"
tax2 <- tax2 %>% 
  mutate( across(.cols = everything(), ~str_replace_all(.,c("XXX"="","XX"="","_X"= "", "_" = " "))))
tax <- rbind(tax1,tax2) %>%
  mutate(across(.cols = everything(), ~str_squish(.))) %>%
  as.matrix() %>% 
  tax_table()

otu1 <- as.matrix(read.delim("../both_lanes/asv/asv_table.txt", row.names = 1))
otu2 <- read.delim("../unpaired_PR2/asv/asv_table.txt", row.names = 1)
rownames(otu2) <- rownames(tax2)
samples_missing <- setdiff(colnames(otu1),colnames(otu2))
for (i in 1:length(samples_missing)) {otu2[samples_missing[i]] <- 0 }
otu2 <- as.matrix (otu2)
otu <-rbind(otu1,otu2) %>%  otu_table(taxa_are_rows = TRUE)
map <- sample_data(read.csv("../mapfile_PE462_corrected_final.csv", row.names = 1, na.strings = c("NA", "")))
physeq_object <-  merge_phyloseq(otu, tax, map)            

# data cleaning 1 ----------------------------------------------------------
# summarize_phyloseq(physeq_object)

# Delete poorly assigned ASVs
get_taxa_unique(physeq_object, "Kingdom") 
physeq_object <- subset_taxa(physeq_object, !is.na(Kingdom) & Kingdom %in% c("Bacteria", "Archaea", "Eukaryota"))


# Unify unassigned --------------------------------------------------------

only_unassigned <- function(x) {
  if_else(x %in% c("uncultured","Uncultured","metagenome","Metagenome","unknown","Unknown","NA"),"unassigned", x)
}

taxo <- as.data.frame(physeq_object@tax_table)%>% 
  replace(is.na(.),"unassigned") %>% 
  mutate( across(.cols = everything(), ~only_unassigned(.)))

taxo <- tax_table(as.matrix(taxo))

physeq_object <- merge_phyloseq(physeq_object@otu_table, taxo, map)

# Delete chloroplast and Mitochodria
any((get_taxa_unique(physeq_object, "Order") == "Chloroplast"))
any((get_taxa_unique(physeq_object, "Family") == "Mitochondria"))
physeq_object <- physeq_object %>% 
  subset_taxa(!Order %in% c("Chloroplast")) %>%
  subset_taxa(!Family %in% c("Mitochondria"))

# alpha diversity ---------------------------------------------------------
#microbiome package
# alpha_tab <-microbiome::alpha(physeq_object, index = "all")
# alpha_tab <- alpha_tab %>% cbind(data.frame(physeq_object@sam_data))
# # write.csv(alpha_tab, file = "../alpha_div_indexes_microbiome_package_with_unassigned_phyla.csv")
# #phyloseq package
# a_div <- phyloseq::estimate_richness(physeq_object)
# write.csv(a_div, file = "../alpha_div_phyloseq.csv")


# data cleaning 2 ---------------------------------------------------------

# keep asv with abundance higher than 2 in at least 3 samples
physeq <- physeq_object %>% 
  filter_taxa( function(x) sum(x > 2) >= 3, TRUE) %>% 
  subset_samples(timepoint_days %nin% c("15")) %>% # Filter out timepoint 15 (incomplete) 
  subset_taxa(!Phylum == "unassigned")

physeq

# a <-microbiome::alpha(physeq, index = "all")
# a <- a %>% cbind(data.frame(physeq@sam_data))
# write.csv(a, file = "../alpha_div_after_cleaning2_no_unassigned_phyla.csv")

# Melt data all -----------------------------------------------------------
source("./tidy_psmelt.R")
tidy_physeq_asv <- tidy_psmelt(physeq)

tidy_physeq_asv$Species <- if_else(
    (!tidy_physeq_asv$Genus=="unassigned"&!tidy_physeq_asv$Species=="unassigned"),
    str_c(tidy_physeq_asv$Genus," ",tidy_physeq_asv$Species),
    str_c(tidy_physeq_asv$Genus," sp."))

  tidy_physeq_asv$detail <-if_else(
  tidy_physeq_asv$material=="wood", 
  str_c(tidy_physeq_asv$material,"_",tidy_physeq_asv$station,"_",tidy_physeq_asv$timepoint_days), 
  tidy_physeq_asv$detail)

tidy_physeq_asv$polymer <- if_else(tidy_physeq_asv$material=="wood", 
                                   str_c(tidy_physeq_asv$material), 
                                   tidy_physeq_asv$polymer)
tidy_physeq_asv$replicate <- NULL

# write_csv(tidy_physeq_asv, "../tidyPE462.csv")

# Merge replicates and calculate relative abundances ---------------------

t3 <- tidy_physeq_asv  %>% group_by(Sample) %>% mutate(Sample_rel_abund = Abundance / sum(Abundance)) %>% #relative abundance of each otu per sample
  ungroup() %>%
  group_by(detail) %>% mutate( rep_rel_abund = Sample_rel_abund / sum(Sample_rel_abund)) %>% #relative abundance of each otu per number of samples in replicates
  ungroup() %>% 
  #Kingdom_section
  group_by(Sample, Kingdom) %>% 
  mutate(Kingdom_rel_abund_Sample = sum(Sample_rel_abund)) %>%  #Kingdom relative abundance per sample 
  ungroup() %>% 
  group_by(detail, Kingdom) %>% 
  mutate(Kingdom_st_dev_abund_samples = sd(Kingdom_rel_abund_Sample)) %>% # standard dev of Kingdom relative abundances between replicates of detail (ployner_timepoint_days_treatment)
  mutate(Kingdom_rep_rel_abund = sum(rep_rel_abund)) %>% #Kingdom relative abundance per samples of desc 
  ungroup() %>%
  #Phylum_section
  group_by(Sample, Phylum) %>% 
  mutate(Phylum_rel_abund_Sample = sum(Sample_rel_abund)) %>%
  ungroup() %>% 
  group_by(detail, Phylum) %>% 
  mutate(st_dev_Phylum_abund = sd(Phylum_rel_abund_Sample)) %>%
  mutate(Phyla_rep_rel_abund = sum(rep_rel_abund)) %>%
  ungroup() %>% 
  #Class_section
  group_by(Sample, Class) %>% 
  mutate(Class_rel_abund_Sample = sum(Sample_rel_abund)) %>%
  ungroup() %>% 
  group_by(detail, Class) %>% 
  mutate(st_dev_Class_abund = sd(Class_rel_abund_Sample)) %>%
  mutate(Class_rep_rel_abund = sum(rep_rel_abund)) %>%
  ungroup() %>% 
  #Order_section
  group_by(Sample, Order) %>% 
  mutate(Order_rel_abund_Sample = sum(Sample_rel_abund)) %>%
  ungroup() %>% 
  group_by(detail, Order) %>% 
  mutate(st_dev_Order_abund = sd(Order_rel_abund_Sample)) %>%
  mutate(Order_rep_rel_abund = sum(rep_rel_abund)) %>%
  ungroup() %>% 
  #Family_section
  group_by(Sample, Family) %>% 
  mutate(Family_rel_abund_Sample = sum(Sample_rel_abund)) %>%
  ungroup() %>% 
  group_by(detail, Family) %>% 
  mutate(st_dev_Family_abund = sd(Family_rel_abund_Sample)) %>%
  mutate(Family_rep_rel_abund = sum(rep_rel_abund)) %>%
  ungroup() %>% 
  #Genus_section
  group_by(Sample, Genus) %>% 
  mutate(Genus_rel_abund_Sample = sum(Sample_rel_abund)) %>%
  ungroup() %>% 
  group_by(detail, Genus) %>% 
  mutate(st_dev_Genus_abund = sd(Genus_rel_abund_Sample)) %>%
  mutate(Genus_rep_rel_abund = sum(rep_rel_abund)) %>%
  ungroup() %>% 
  #Species_section
  group_by(Sample, Species) %>% 
  mutate(Species_rel_abund_Sample = sum(Sample_rel_abund)) %>%
  ungroup() %>% 
  group_by(detail, Species) %>% 
  mutate(st_dev_Species_abund = sd(Species_rel_abund_Sample)) %>%
  mutate(Species_rep_rel_abund = sum(rep_rel_abund)) %>%
  ungroup()  

# write_csv(t3, "../tidyPE462_rel_abund_calc.csv")