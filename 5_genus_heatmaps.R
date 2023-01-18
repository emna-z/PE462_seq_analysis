##################################################
## Project: PE462
## Script purpose: Negative controls exploration
## Date: 14-01-2023
## Author: E. Zeghal
##################################################


library("tidyverse")
library("Hmisc")
library("RColorBrewer")
library("mapsf")
library("ggtext")
library("glue")
library("ggpubr")
library("ggh4x")

# k <- read_csv("../Kingdom_PE462.csv")
# p <- read_csv("../phyla_PE462.csv")
# o <-read_csv("../order_PE462.csv")


# load genera table -------------------------------------------------------

g_all <- read_csv("../genus_PE462.csv")
g_all <- g_all %>% filter(timepoint_days%nin% c("mock_DNA","Neg_PCR")) %>% 
  filter(station%nin% c("mock_DNA","Neg_PCR"))%>% 
  mutate(station = str_replace_all(station,c("C13" = "coastal station"))) %>% 
  mutate(station = str_replace_all(station,c("C05"= "open water station"))) 


# extract top 5 genera per condition --------------------------------------

top5 <- g_all%>% 
  filter(Genus %nin% c("unassigned")) %>% 
  mutate(across(c(timepoint_days, detail),factor))%>% distinct() %>% 
  group_by(detail) %>% slice_max(order_by = Genus_rep_rel_abund, n = 4)

top5 <- top5 %>% 
  filter(Genus %nin% c("unassigned"))


# Gradient palettes -------------------------------------------------------

blu_grad <- c("#0F2A5B","#253E6A","#3C5178","#526587","#687996",
                       "#7F8CA4","#95A0B3","#ABB4C2","#C2C7D0","grey95") %>%
                         rev()

# Recode treatment --------------------------------------------------------

breaks <- levels(as.factor(top5$Genus))

hcb <- read_lines("../Hydrocarbon_degraders_sorted_22_08.txt")
pdb <- read_lines("../PlasticDB_HCB_genera.txt")

condition <- if_else(breaks %in% intersect(pdb,hcb) ,"#6712bf" ,
                     (if_else(breaks %in% pdb , "#12bf67", (if_else(breaks %in% hcb , "#ff5722", "black")))))

# Coastal station ---------------------------------------------------------

coast <- g_all %>% filter(station %in% c("coastal station")) %>% 
  filter (Genus %in% (unique(top5$Genus)))

# just to add t0 to all polymers, glass_fiber needs to be replaced by each plastic polymer factor and duplicated
pe0 <- coast %>% filter(polymer %in% c("Glass_fiber")) %>% 
  mutate(polymer = str_replace_all(polymer,c("Glass_fiber"="PE"))) %>% 
  mutate(treatment = str_replace_all(treatment,c( "no" ="no_UV"))) 
peuv0 <-  coast %>% filter(polymer %in% c("Glass_fiber")) %>% 
  mutate(polymer = str_replace_all(polymer,c("Glass_fiber"="PE"))) %>% 
  mutate(treatment = str_replace_all(treatment,c( "no" ="UV")))
pet0 <- coast %>% filter(polymer %in% c("Glass_fiber")) %>% 
  mutate(polymer = str_replace_all(polymer,c("Glass_fiber"="PET"))) %>% 
  mutate(treatment = str_replace_all(treatment,c( "no" ="no_UV"))) 
petuv0 <-  coast %>% filter(polymer %in% c("Glass_fiber")) %>% 
  mutate(polymer = str_replace_all(polymer,c("Glass_fiber"="PET"))) %>% 
  mutate(treatment = str_replace_all(treatment,c( "no" ="UV")))
ps0 <- coast %>% filter(polymer %in% c("Glass_fiber")) %>% 
  mutate(polymer = str_replace_all(polymer,c("Glass_fiber"="PS"))) %>% 
  mutate(treatment = str_replace_all(treatment,c( "no" ="no_UV"))) 
psuv0 <-  coast %>% filter(polymer %in% c("Glass_fiber")) %>% 
  mutate(polymer = str_replace_all(polymer,c("Glass_fiber"="PS"))) %>% 
  mutate(treatment = str_replace_all(treatment,c( "no" ="UV")))
nylon0 <- coast %>% filter(polymer %in% c("Glass_fiber")) %>% 
  mutate(polymer = str_replace_all(polymer,c("Glass_fiber"="Nylon"))) %>% 
  mutate(treatment = str_replace_all(treatment,c( "no" ="no_UV"))) 
nylonuv0 <-  coast %>% filter(polymer %in% c("Glass_fiber")) %>% 
  mutate(polymer = str_replace_all(polymer,c("Glass_fiber"="Nylon"))) %>% 
  mutate(treatment = str_replace_all(treatment,c( "no" ="UV")))
wood0 <-  coast %>% filter(polymer %in% c("Glass_fiber")) %>% 
  mutate(polymer = str_replace_all(polymer,c("Glass_fiber"="wood")))

# join all duplicated t0 to original table
coast2 <- coast %>% filter(polymer %nin% c("Glass_fiber")) %>%
  bind_rows (pe0, peuv0, pet0, petuv0, ps0, psuv0, nylon0, nylonuv0, wood0) %>% 
  mutate(treatment = if_else(grepl( "no", treatment), "none", treatment)) %>% 
  mutate_if(is.character,as.factor) 



# heatmap genera coastal station ------------------------------------------

coast_h <- ggplot(data = coast2, mapping = aes( x = fct_relevel(timepoint_days,  c("0","5", "10","30","45")),
                                          y = fct_relevel(Genus,rev), fill = Genus_rep_rel_abund)) +

    geom_tile() +
    scale_fill_gradientn(name = "RA", colours = blu_grad, limits=c(0,0.75) )+ 
    facet_nested( ~ polymer+treatment, nest_line = element_line(linetype = 1, linewidth = 1)) +
    xlab(label = "incubation time (days)") +
    scale_y_discrete("Genus", breaks = breaks) +
    theme_minimal() +
    theme(strip.placement = "outside",
          plot.title = element_text(hjust = 0.5),
          axis.title.y = element_blank(),
          axis.text.y = element_text(face = "bold.italic", colour = condition),
          axis.text.x = element_text(face = "bold"),
          strip.background = element_blank(),
          ggh4x.facet.nestline = element_line(colour = c("#CC6677")),
          strip.text.x = element_text(face = "bold"),
          legend.title = element_text(face = "bold"),
          legend.position="right")

# coast_h


# Open water station ------------------------------------------------------

os <- g_all %>% filter(station %in% c("open water station")) %>% 
  filter (Genus %in% (unique(top5$Genus)))

pe0 <- os %>% filter(polymer %in% c("Glass_fiber")) %>% 
  mutate(polymer = str_replace_all(polymer,c("Glass_fiber"="PE"))) %>% 
  mutate(treatment = str_replace_all(treatment,c( "no" ="no_UV"))) 
peuv0 <-  os %>% filter(polymer %in% c("Glass_fiber")) %>% 
  mutate(polymer = str_replace_all(polymer,c("Glass_fiber"="PE"))) %>% 
  mutate(treatment = str_replace_all(treatment,c( "no" ="UV")))
pet0 <- os %>% filter(polymer %in% c("Glass_fiber")) %>% 
  mutate(polymer = str_replace_all(polymer,c("Glass_fiber"="PET"))) %>% 
  mutate(treatment = str_replace_all(treatment,c( "no" ="no_UV"))) 
petuv0 <-  os %>% filter(polymer %in% c("Glass_fiber")) %>% 
  mutate(polymer = str_replace_all(polymer,c("Glass_fiber"="PET"))) %>% 
  mutate(treatment = str_replace_all(treatment,c( "no" ="UV")))
ps0 <- os %>% filter(polymer %in% c("Glass_fiber")) %>% 
  mutate(polymer = str_replace_all(polymer,c("Glass_fiber"="PS"))) %>% 
  mutate(treatment = str_replace_all(treatment,c( "no" ="no_UV"))) 
psuv0 <-  os %>% filter(polymer %in% c("Glass_fiber")) %>% 
  mutate(polymer = str_replace_all(polymer,c("Glass_fiber"="PS"))) %>% 
  mutate(treatment = str_replace_all(treatment,c( "no" ="UV")))
nylon0 <- os %>% filter(polymer %in% c("Glass_fiber")) %>% 
  mutate(polymer = str_replace_all(polymer,c("Glass_fiber"="Nylon"))) %>% 
  mutate(treatment = str_replace_all(treatment,c( "no" ="no_UV"))) 
nylonuv0 <-  os %>% filter(polymer %in% c("Glass_fiber")) %>% 
  mutate(polymer = str_replace_all(polymer,c("Glass_fiber"="Nylon"))) %>% 
  mutate(treatment = str_replace_all(treatment,c( "no" ="UV")))
wood0 <-  os %>% filter(polymer %in% c("Glass_fiber")) %>% 
  mutate(polymer = str_replace_all(polymer,c("Glass_fiber"="wood")))

os2 <- os %>% filter(polymer %nin% c("Glass_fiber")) %>%
  bind_rows (pe0, peuv0, pet0, petuv0, ps0, psuv0, nylon0, nylonuv0, wood0) %>% 
  mutate(treatment = if_else(grepl( "no", treatment), "none", treatment)) %>% 
  mutate_if(is.character,as.factor)


# heatmap genera open water -----------------------------------------------


openw_h <- ggplot(data = coast2, mapping = aes( x = fct_relevel(timepoint_days,  c("0","5", "10","30","45")),
                                                y = fct_relevel(Genus,rev), fill = Genus_rep_rel_abund)) +
  
  geom_tile() +
  scale_fill_gradientn(name = "RA", colours = blu_grad, limits=c(0,0.75) )+ 
  facet_nested( ~ polymer+treatment, nest_line = element_line(linetype = 1, linewidth = 1)) +
  xlab(label = "incubation time (days)") +
  scale_y_discrete("Genus", breaks = breaks) +
  theme_minimal() +
  theme(strip.placement = "outside",
        plot.title = element_text(hjust = 0.5),
        axis.title.y = element_blank(),
        axis.text.y = element_text(face = "bold.italic", colour = condition),
        axis.text.x = element_text(face = "bold"),
        strip.background = element_blank(),
        ggh4x.facet.nestline = element_line(colour = c("#CC6677")),
        strip.text.x = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.position="right") 
# openw_h



# all conditions top 4 tab ------------------------------------------------
# 
# top4_all <- bind_rows(os2,coast2) 
# top4_h <- ggplot(data = top4_all, mapping = aes( x = fct_relevel(timepoint_days,  c("0","5", "10","30","45")),
#                                                 y = fct_relevel(Genus,rev), fill = Genus_rep_rel_abund)) +
#   
#   geom_tile() +
#   scale_fill_gradientn(name = "RA", colours = blu_grad, limits=c(0,0.75) )+ 
#   facet_nested( ~ station+polymer+treatment, nest_line = element_line(linetype = 1, linewidth = 1)) +
#   xlab(label = "incubation time (days)") +
#   scale_y_discrete("Genus", breaks = breaks) +
#   theme_minimal() +
#   theme(strip.placement = "outside",
#         plot.title = element_text(hjust = 0.5),
#         axis.title.y = element_blank(),
#         axis.text.y = element_text(face = "bold.italic", colour = condition),
#         axis.text.x = element_text(face = "bold"),
#         strip.background = element_blank(),
#         ggh4x.facet.nestline = element_line(colour = c("#CC6677")),
#         strip.text.x = element_text(face = "bold"),
#         legend.title = element_text(face = "bold"),
#         legend.position="top") 
# 
# top4_h
