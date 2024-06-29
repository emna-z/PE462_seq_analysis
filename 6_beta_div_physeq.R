##################################################
## Project: PE462
## Script purpose: Beta diversity, 
## ordination and community exploration
## Date: 25-01-2023
## Author: E. Zeghal
##################################################


# libraries ---------------------------------------------------------------

library("tidyverse")
library("Hmisc")
library("vegan")
library("mia")
library("scater")
library("phyloseq")
library("metagMisc")
library("RColorBrewer")


set.seed(42)



# Load data ---------------------------------------------------------------

# load the final phyloseq object from which the tidy table was created
source("1_data_import.R")

#the tax_table, otu_table and sample_data from that object will be corrected in the same way the tidy table was

## correct and load tax table ---------------------------------------------

# tax_physeq <- as.data.frame(physeq@tax_table)
# colnames(tax_physeq) <- colnames(physeq@tax_table@.Data)
# tax_physeq$Species <- if_else(
#   (!tax_physeq$Genus=="unassigned"&!tax_physeq$Species=="unassigned"),
#   str_c(tax_physeq$Genus," ",tax_physeq$Species),
#   str_c(tax_physeq$Genus," sp."))
# write.csv(tax_physeq, "../tax_tab_physeq.csv")

tax_pseq <- read.csv("../tax_tab_physeq.csv", row.names = 1) %>%
  as.matrix() %>%
  tax_table()

## correct and load asv table ---------------------------------------------

# asv_physeq <- as.data.frame(physeq@otu_table)
# colnames(asv_physeq) <- colnames(physeq@otu_table@.Data)
# write.csv(asv_physeq, "../asv_tab_physeq.csv")

asv_pseq <- otu_table(as.matrix(read.csv("../asv_tab_physeq.csv", row.names = 1)), taxa_are_rows = T)

## correct and load sample data table -------------------------------------

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
                ~str_replace_all(.,c("C05"= "open NS water", "C13" = "coastal NS water", "_"= " "))))

map_pseq <- sample_data(map)

# phyloseq object ---------------------------------------------------------

pseq <- merge_phyloseq(asv_pseq, tax_pseq, map_pseq) 

pseq_neg <- pseq %>% 
  subset_samples(polymer %in% c("Neg_PCR"))

pseq <- pseq %>% 
  subset_samples(timepoint_days %nin% c("15")) %>% 
  subset_samples(station %in% c("C05","C13")) %>%
  subset_samples(timepoint_days %nin% c("15","0")) %>% 
  subset_taxa(!Phylum == "unassigned") 

pseq@sam_data$timepoint_days <- fct_relevel(pseq@sam_data$timepoint_days ,c("5", "10", "30", "45"))

pseq_no_wood <-  pseq %>% 
  subset_samples(material %nin% c("wood"))

# CLR transformation ------------------------------------------------------

ps_clr <- microbiome::transform(pseq, "clr", pseudocount=1)


ordin <- ordinate(ps_clr,  "RDA")

p <-  plot_ordination (ps_clr, ordin, color = "polymer_station", shape = "timepoint_days")

# brewer.pal(10,"Paired") %>% rev()
pal <- c("#6A3D9A", "#CAB2D6", "#FF7F00", "#FDBF6F", "#E31A1C",
                  "#FB9A99", "#33A02C", "#B2DF8A", "#1F78B4", "#A6CEE3")

p+geom_point(size = 3)+
  scale_shape(name="days")+
  scale_color_manual(name="polymer/station", values = pal)+
  xlim(-5, 5)+
  ylim(-5, 5)+
  theme_minimal()+
  theme(strip.text.x = element_text(size=12, face="bold"),
        strip.text.y = element_text(size=12, face="bold"),
        legend.title = element_text(face = "bold", size=12),
        legend.text = element_text(face = "bold", size=11),
        axis.title.y = element_text(face = "bold", size=12),
        axis.title.x = element_text(face = "bold", size=12),
        axis.text.y = element_text(face = "bold"),
        axis.text.x = element_text(face = "bold"))+
  guides(color = guide_legend(
    override.aes=list(shape = 18)))


# RCLR transform ----------------------------------------------------------

# ps_rclr <- microbiome::transform(pseq, "rclr")
# 
# ordin <- ordinate(ps_rclr,  "RDA","robust.aitchison")
# 
# p <-  plot_ordination (ps_rclr, ordin, color = "polymer_station", shape = "timepoint_days")
# pal <- c("#f7d694","#c9aaef","#8e400c","#2f5fc6","#6d7745","#44c0ff","#da9055","#4b0082","#daed89","#3a1493")
# 
# p+geom_point(size = 3)+
#   scale_shape(name="days")+
#   scale_color_manual(name="polymer/station", values = pal)+
#   theme_minimal()+ggtitle(label = "PCA overall - rClr transformed - raitchison distance", subtitle = "") +
#   theme(strip.text.x = element_text(size=12, face="bold"),strip.text.y = element_text(size=12, face="bold"),
#         legend.title = element_text(face = "bold"), 
#         legend.title.align = 0.4,
#         legend.text = element_text(face = "bold"),
#         axis.title.y = element_text(face = "bold"),
#         axis.title.x = element_text(face = "bold"),
#         axis.text.y = element_text(face = "bold"),
#         axis.text.x = element_text(face = "bold"))+
#   guides(color = guide_legend(
#     override.aes=list(shape = 18)))

# No Wood -----------------------------------------------------------------

ps_clr_no_wood <- microbiome::transform(pseq_no_wood, "clr")

ordin <- ordinate(ps_clr_no_wood,  "RDA","robust.aitchison")

p <-  plot_ordination (ps_clr_no_wood, ordin, color = "polymer_station", shape = "timepoint_days")

brewer.pal(10,"Paired") %>% rev()
pal <- c("#6A3D9A", "#CAB2D6", "#FF7F00", "#FDBF6F", "#E31A1C",
                  "#FB9A99", "#33A02C", "#B2DF8A", "#1F78B4", "#A6CEE3")

p+geom_point(size = 3)+
  # stat_ellipse(aes (group=interaction(pseq_no_wood@sam_data$timepoint_days,pseq_no_wood@sam_data$station)),  linetype = 2) +
  scale_shape(name="days")+
  scale_color_manual(name="polymer/station", values = pal)+
  theme_minimal()+ggtitle(label = "PCA only plastics - robust aitchison distance", subtitle = "") +
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


# dispersion pseq ---------------------------------------------------------
# aitch <- vegdist(t(pseq@otu_table@.Data), method = "aitchison", pseudocount = 1 )

d1 <- distance(ps_clr_no_wood, method = "euclidean")
disp_aitch <- betadisper(d1, ps_clr_no_wood@sam_data$timepoint_days)
plot(disp_aitch)
boxplot(disp_aitch, main = "", xlab = "")
anova(disp_aitch)

perm <- permutest(disp_aitch, pairwise = TRUE)
perm
HSD <- TukeyHSD(disp_aitch)
HSD
plot(HSD)
pstat <- permustats(perm)
densityplot(pstat, scales = list(x = list(relation = "free")))
qqmath(pstat, scales = list(relation = "free"))


# ANOSIM ------------------------------------------------------------------
anosim_time <- with(ps_clr_no_wood@sam_data, anosim(d1, timepoint_days))
summary(anosim_time)

anosim_station <- with(ps_clr_no_wood@sam_data, anosim(d1, station))
summary(anosim_station)

anosim_polymer <- with(ps_clr_no_wood@sam_data, anosim(d1, polymer))
summary(anosim_polymer)

anosim_treatment <- with(ps_clr_no_wood@sam_data, anosim(d1, treatment))
summary(anosim_treatment)
# disp pseq no wood -------------------------------------------------------

# The p-value of permutest < 0.05 ---> the difference in variance is significant


# aitch <- vegdist(t(pseq_no_wood@otu_table@.Data), method = "aitchison", pseudocount = 1 )
# disp_aitch <- betadisper(aitch, pseq_no_wood@sam_data$polymer)
# plot(disp_aitch)
# anova(disp_aitch)
# perm <- permutest(disp_aitch, pairwise = TRUE)
# perm
# HSD <- TukeyHSD(disp_aitch)
# plot(HSD)
# pstat <- permustats(perm)
# densityplot(pstat, scales = list(x = list(relation = "free")))
# qqmath(pstat, scales = list(relation = "free"))

#robust.aitchison
# raitch <- vegdist(t(pseq_no_wood@otu_table@.Data), method = "robust.aitchison")
# disp_raitch <- betadisper(raitch, pseq_no_wood@sam_data$polymer)
# plot(disp_raitch)
# 
# anova(disp_raitch, p.adjust.method="BH")
# 
# perm <- permutest(disp_raitch, pairwise = TRUE, p.adjust.method="BH")
# perm
# HSD <- TukeyHSD(disp_raitch, p.adjust.method="BH")
# plot(HSD)
# pstat <- permustats(perm)
# summary(pstat)
# densityplot(pstat, scales = list(x = list(relation = "free")))
# qqmath(pstat, scales = list(relation = "free"))

# PERMANOVA ---------------------------------------------------------------

df <- as(sample_data(ps_clr), "data.frame")

# aitchison
d1 <- distance(ps_clr, method = "euclidean")

permanova1 <- adonis2(d1 ~ timepoint_days+station+polymer+treatment, df, by = "margin"  )
# permanova1 

# adonis2(raitch ~ timepoint_days+station+polymer+treatment, df, by = "margin"  )
# adonis2(raitch ~ timepoint_days+station+polymer_photo, df, by = "margin"  )
# adonis2(raitch ~ timepoint_days, df, by = "margin"  )
# adonis2(raitch ~ station, df, by = "margin"  )
# adonis2(raitch ~ polymer, df, by = "margin"  )
# adonis2(raitch ~ treatment, df, by = "margin"  )


df <- as(sample_data(ps_clr_no_wood), "data.frame")

# aitchison
d1 <- distance(ps_clr_no_wood, method = "euclidean")

permanova1 <- adonis2(d1 ~ timepoint_days+station+polymer+treatment, df, by = "margin"  )
permanova1

a1 <- glm(d1 ~ timepoint_days + station + polymer + treatment, family= "poisson", df)
summary(a1)

# adonis2(raitch ~ timepoint_days+station+polymer+treatment, df, by = "margin"  )
# adonis2(raitch ~ timepoint_days+station+polymer_photo, df, by = "margin"  )
# adonis2(raitch ~ timepoint_days, df, by = "margin"  )
# adonis2(raitch ~ station, df, by = "margin"  )
# adonis2(raitch ~ polymer, df, by = "margin"  )
# adonis2(raitch ~ treatment, df, by = "margin"  )

# robust aitchison
# d2 <- vegdist(t(pseq_no_wood@otu_table@.Data), method = "robust.aitchison")
# permanova2 <- adonis2(d2 ~ timepoint_days+station+polymer+treatment, df, by = "margin"  )
# permanova2 
# 
# 
# 
# # time --------------------------------------------------------------------
# 
# pseq_no_wood_early <-  pseq_no_wood %>%
#   subset_samples(timepoint_days %in% c("45"))
# 
# raitch <- vegdist(t(pseq_no_wood_early@otu_table@.Data), method = "robust.aitchison")
# 
# df <- as(sample_data(pseq_no_wood_early), "data.frame")
# 
# 
# adonis2(raitch ~ timepoint_days+station+polymer+treatment, df, by = "margin"  )
# adonis2(raitch ~ timepoint_days+station+polymer_photo, df, by = "margin"  )
# adonis2(raitch ~ timepoint_days, df, by = "margin"  )
# adonis2(raitch ~ station, df, by = "margin"  )
# adonis2(raitch ~ polymer, df, by = "margin"  )
# adonis2(raitch ~ treatment, df, by = "margin"  )

