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
library("Hmisc")

# raw data import ---------------------------------------------------------

tax1 <- read.delim("../both_lanes/asv/taxonomy_dada2/representative_seq_set_tax_assignments.txt",
                             row.names = 1, na.strings = "NA")
tax2 <- read.delim("../unpaired_PR2/asv/taxonomy_dada2/representative_seq_set_tax_assignments.txt",
                   row.names = 1, na.strings = "NA")
# tax2 <- tax2 %>% filter(!is.na(Kingdom) & Kingdom == "Eukaryota") #gets rid of "Eukaryota:plas", "Bacteria" and NA
rownames(tax2) <- paste0("asv.",(length(rownames(tax1))+1):(length(rownames(tax1))+length(rownames(tax2))))

# for ease of handling we'll change the rank Division to Phylum for the Eukaryotes table
# and add a Supergroup column to the prokaryotes column by duplicating the kingdom column content
ranks <- c ("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
colnames(tax1) <-  ranks
tax1 <- tax1 %>% mutate(Supergroup = Kingdom)
colnames(tax2)[3] <- "Phylum"
tax2 <- as.matrix(tax2)
tax1 <- as.matrix(tax1)
tax <- rbind(tax1,tax2) %>% tax_table()

otu1 <- as.matrix(read.delim("../both_lanes/asv/asv_table.txt", row.names = 1))
otu2 <- read.delim("../unpaired_PR2/asv/asv_table.txt", row.names = 1)
rownames(otu2) <- rownames(tax2)
samples_missing <- setdiff(colnames(otu1),colnames(otu2))
for (i in 1:length(samples_missing)) {otu2[samples_missing[i]] <- 0 }
otu2 <- as.matrix (otu2)
otu <-rbind(otu1,otu2) %>%  otu_table(taxa_are_rows = TRUE)
map <- sample_data(read.csv("../mapfile_PE462_corrected_final.csv", row.names = 1, na.strings = c("NA", "")))
physeq_object <-  merge_phyloseq(otu, tax, map)            
