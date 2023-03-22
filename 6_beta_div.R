library("mia")
# library("MicrobiomeStat")
library("tidyverse")
library("Hmisc")
library("vegan")
library("scater")
library("phyloseq")
library("FactoMineR")
library("factoextra")
# library("metagMisc")
# library("finalfit")
# library("kableExtra")
# library("knitr")

set.seed(42)

# source("1_data_import.R")
# tax_physeq <- as.data.frame(physeq@tax_table)
# colnames(tax_physeq) <- colnames(physeq@tax_table@.Data)
# tax_physeq$Species <- if_else(!tax_physeq$Species=="unassigned",
#                               str_c(tax_physeq$Genus," ",tax_physeq$Species), tax_physeq$Species)
# write.csv(tax_physeq, "../tax_tab_physeq.csv")

tax_pseq <- read.csv("../tax_tab_physeq.csv", row.names = 1) %>%
  as.matrix() %>%
  tax_table()

# asv_physeq <- as.data.frame(physeq@otu_table)
# colnames(asv_physeq) <- colnames(physeq@otu_table@.Data)
# write.csv(asv_physeq, "../asv_tab_physeq.csv")

asv_pseq <- otu_table(as.matrix(read.csv("../asv_tab_physeq.csv", row.names = 1)), taxa_are_rows = T)

# map <- read.csv("../mapfile_PE462_corrected_final.csv", na.strings = c("NA", ""))
# map$detail <- if_else(map$material=="wood", str_c(map$material,"_",map$station,"_",map$timepoint_days), map$detail)
# map$polymer <- if_else(map$material=="wood", str_c(map$material), map$polymer)
# polymer_station <-  str_c(map$polymer, "_", map$station)
# polymer_photo <- str_c(map$polymer, "_", map$treatment)
# pol_photo_station <- str_c(map$polymer, "_", map$treatment, "_", map$station)
# map <- map %>% add_column(polymer_station) %>% add_column(polymer_photo) %>% add_column(pol_photo_station)
# map$polymer_station <- if_else(map$material=="wood", str_c(map$material,"_",map$station), map$polymer_station)
# map$polymer_photo <- if_else(map$material=="wood", str_c(map$material), map$polymer_photo)
# map$pol_photo_station <- if_else(map$material=="wood", str_c(map$material,"_",map$station), map$pol_photo_station)
# write_csv(map,"../map_from_tidy.csv")

map <- read.csv("../map_from_tidy.csv",row.names = 1, na.strings = c("NA", ""), stringsAsFactors = T) %>% 
  mutate(across(.cols = c("polymer_station","pol_photo_station"), 
                      ~str_replace_all(.,c("C05"= "open_water", "C13" = "coast"))))

map_pseq <- sample_data(map)

pseq <- merge_phyloseq(asv_pseq, tax_pseq, map_pseq) 

pseq_neg <- pseq %>% 
  subset_samples(polymer %in% c("Neg_PCR"))

pseq <- pseq %>% 
  subset_samples(station %in% c("C05","C13")) %>%
  subset_samples(timepoint_days %nin% c("15","0"))
pseq@sam_data$timepoint_days <- fct_relevel(pseq@sam_data$timepoint_days ,c("5", "10", "30", "45"))
  
# pseq.rel <- microbiome::transform(pseq, "compositional")

pseq <- transform_sample_counts(pseq, function(x) x/sum(x))
ps_clr <- microbiome::transform(pseq, "clr")
iDist <- distance(pseq,"bray")

ordin <- ordinate(ps_clr, "PCoA", "jaccard")
pal <- c("#f7d694","#c9aaef","#8e400c","#2f5fc6","#6d7745","#44c0ff","#da9055","#2f8d7e","#daed89","#3a1493")
Tol_muted <- c('#DDDDDD','#88CCEE', '#44AA99', '#117733', '#332288', '#DDCC77', '#999933','#CC6677', '#882255', '#AA4499')
p <-  plot_ordination (pseq, ordin, color = "polymer_station", shape = "timepoint_days", label = "station")


p+geom_point(size = 2.5)+
  scale_shape(name="days")+
  scale_color_manual(values = pal)+
  theme_minimal()+ggtitle(label = "PCoA on overall dataset - ASV level", subtitle = "distance: Bray-Curtis") +
  theme(strip.text.x = element_text(size=12, face="bold"),strip.text.y = element_text(size=12, face="bold"),
        legend.title = element_text(face = "bold"), 
        legend.title.align = 0.4,
        legend.text = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold"),
        axis.text.y = element_text(face = "bold"),
        axis.text.x = element_text(face = "bold"))


ordin2 <- ordinate(ps_clr,  "RDA","euclidean")

p <-  plot_ordination (pseq, ordin2, color = "polymer_station", shape = "timepoint_days")
Tol_muted <- c('#88CCEE', '#44AA99', '#117733', '#332288', '#DDCC77', '#999933','#CC6677', '#882255', '#AA4499', '#DDDDDD')

p+geom_point(size = 3)+
  scale_shape(name="days")+
  scale_color_manual(values = pal)+
  theme_minimal()+ggtitle(label = "PCA overall", subtitle = "Clr transformed - euclidean distance") +
  theme(strip.text.x = element_text(size=12, face="bold"),strip.text.y = element_text(size=12, face="bold"),
        legend.title = element_text(face = "bold"), 
        legend.title.align = 0.4,
        legend.text = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold"),
        axis.text.y = element_text(face = "bold"),
        axis.text.x = element_text(face = "bold"))+
  guides(color = guide_legend(
    override.aes=list(shape = 18)))

ps_clr <- microbiome::transform(pseq, "clr")
ps_hel <- microbiome::transform(pseq, "hellinger")
ord_clr <- phyloseq::ordinate(ps_clr, "RDA")
phyloseq::plot_scree(ord_clr) + 
  geom_bar(stat="identity", fill = "blue") +
  labs(x = "\nAxis", y = "Proportion of Variance\n")

clr1 <- ord_clr$CA$eig[1] / sum(ord_clr$CA$eig)
clr2 <- ord_clr$CA$eig[2] / sum(ord_clr$CA$eig)
phyloseq::plot_ordination(pseq, ord_clr, type="station", color="timepoint_days") + 
  geom_point(size = 2) +
  coord_fixed(clr2 / clr1) +
  stat_ellipse(aes(group =timepoint_days ), linetype = 2) + facet_grid(~station)





df  <-  as(sample_data(pseq), "data.frame")
d  <-  distance(ps_clr, method = "aitchison")
adon <- adonis2(d ~  timepoint_days+station+ polymer+treatment, df, by = "margin"  )
adon 
# write.csv(adon, "adonis2_margin_clr_euclidean.csv")
adon_seq <-  adonis2(d ~  timepoint_days+station+ polymer+treatment, df, by = "terms" )
adon_seq

#############time effect##################
pseq_coast <- subset_samples(ps_clr, station %in% c("C13"))
d_coast <- vegdist(pseq_coast, "euclidean")
adon <- adonis2(d_coast ~  timepoint_days+polymer+treatment,(as(sample_data(pseq_coast), "data.frame")), by = "margin")
adon 
pseq_os <- subset_samples(ps_clr, station %in% c("C05"))
d_os <- distance(pseq_os, "euclidean")
adon <- adonis2(d_os ~  timepoint_days+polymer+treatment,(as(sample_data(pseq_os), "data.frame")), by = "margin")
adon 
base::source("./permanova_pairwise.R")
l1 <- permanova_pairwise(d_os, pseq_os@sam_data$timepoint_days )
kable(l1, "pipe")
l2 <- permanova_pairwise(d_coast, pseq_coast@sam_data$timepoint_days )
kable(l2, "pipe")

l1$pval <- cell_spec(l1$pval, background_as_tile = TRUE, background = if_else((!is.na(l1$pval) & l1$pval<0.05),"#91EFDE",if_else((!is.na(l1$pval) & l1$pval>0.05 & l1$pval<0.1),"#EFE049" ,"#BBBBBB" )))
l1$p.adj <- cell_spec(l1$p.adj, background_as_tile = TRUE, background = if_else((!is.na(l1$p.adj) & l1$p.adj<0.05),"#91EFDE",if_else((!is.na(l1$p.adj) & l1$p.adj>0.05 & l1$p.adj<0.1),"#EFE049" ,"#BBBBBB" )))
kbl(l1, escape = F, align = "c") %>%
  kable_classic("basic", full_width = F, html_font = "Arial")

################
pseq_coast <- subset_samples(pseq_coast, polymer %nin% c("wood"))
pseq_os <- subset_samples(pseq_os, polymer %nin% c("wood"))



lc <- metagMisc::phyloseq_sep_variable(pseq_coast, "timepoint_days", drop_zeroes = F)
lc2 <- lapply(seq_along(lc), function(i)
  metagMisc::phyloseq_sep_variable(lc[[i]], "polymer", drop_zeroes = F))

lc
samdata_c <- lapply(seq_along(lc), function(i)
  as(sample_data(lc[[i]]), "data.frame"))

dist_c <- lapply(seq_along(lc), function(i)
  distance(lc[[i]], "bray"))
perm_pol_c <- lapply(seq_along(dist_c), function(i)
  permanova_pairwise(dist_c[[i]], lc[[i]]@sam_data$treatment ))# treatment ain't significant

########################################
#Generate data.frame with OTUs and metadata
ps_wilcox <- data.frame(t(data.frame(phyloseq::otu_table(ps_clr))))
ps_wilcox$Status <- phyloseq::sample_data(ps_clr)$detail
#Define functions to pass to map
wilcox_model <- function(df){
  wilcox.test(abund ~ Status, data = df)
}
wilcox_pval <- function(df){
  wilcox.test(abund ~ Status, data = df)$p.value
}
#Create nested data frames by OTU and loop over each using map 
wilcox_results <- ps_wilcox %>%
  gather(key = OTU, value = abund, -Status) %>%
  group_by(OTU) %>%
  nest() %>%
  mutate(wilcox_test = map(data, wilcox_model),
         p_value = map(data, wilcox_pval))                       
#Show results
head(wilcox_results)