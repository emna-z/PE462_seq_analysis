##################################################
## Project: PE462
## Script purpose: heatmaps for genera RA
## Date: 14-01-2023
## Author: E. Zeghal
##################################################


# load libraries ----------------------------------------------------------

library("tidyverse")
library("Hmisc")
library("scales")
library("RColorBrewer")
library("mapsf")
library("ggtext")
library("glue")
library("ggpubr")
library("ggh4x")


# load genera table -------------------------------------------------------

g_all <- read_csv("../genus_PE462.csv")
g_all <- g_all %>% filter(timepoint_days%nin% c("mock_DNA","Neg_PCR")) %>% 
  filter(station%nin% c("mock_DNA","Neg_PCR"))%>% 
  mutate(polymer = fct_relevel(polymer, c("PE", "PS", "Nylon", "PET", "wood" ))) %>% 
  mutate(station = str_replace_all(station,c("C13" = "coastal NS water"))) %>% 
  mutate(station = str_replace_all(station,c("C05"= "open NS water"))) 


# extract top 3 genera per condition --------------------------------------

top3 <- g_all%>% 
  filter(Genus %nin% c("unassigned")) %>% 
  mutate(across(c(timepoint_days, detail),factor))%>% distinct() %>% 
  group_by(detail) %>% slice_max(order_by = Genus_rep_rel_abund, n = 3)
top3 <- top3 %>% 
  filter(Genus %nin% c("unassigned"))


# Gradient palettes -------------------------------------------------------
# display.brewer.all()
# 
# pal <- colorRampPalette(c("#E69F00", "#0072B2", "#5913bc"))(30)
# 
# cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
#                    "#F0E442", "#56B4E9", "#D55E00", "#CC79A7")
# 
# blu_grad <- c("#0F2A5B","#253E6A","#3C5178","#526587","#687996",
#                        "#7F8CA4","#95A0B3","#ABB4C2","#C2C7D0","grey95") %>%
#                          rev()

# Recode treatment --------------------------------------------------------

breaks <- levels(as.factor(top3$Genus))

hcb <- read_lines("../Hydrocarbon_degraders_sorted_22_08.txt")
pdb <- read_lines("../../genera_names_plasticDB.txt")
pdb <- unique(pdb)

condition <- if_else(breaks %in% intersect(pdb,hcb) ,"#b34772" ,
                     (if_else(breaks %in% pdb , "#009E73", (if_else(breaks %in% hcb , "#0072B2", "black")))))

# Coastal station ---------------------------------------------------------

coast <- g_all %>% filter(station %in% c("coastal NS water")) %>% 
  filter (Genus %in% (unique(top3$Genus)))

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
  mutate(treatment = if_else(grepl( "no", treatment), "no UV", treatment)) %>% 
  mutate_if(is.character,as.factor) 


# Open water station ------------------------------------------------------

os <- g_all %>% filter(station %in% c("open NS water")) %>% 
  filter (Genus %in% (unique(top3$Genus)))

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
  mutate(treatment = if_else(grepl( "no", treatment), "no UV", treatment)) %>% 
  mutate_if(is.character,as.factor)

# full and subsets data tibble --------------------------------------------

full <- bind_rows(coast2, os2)
plastics <- full %>% filter(polymer %in% c("PE", "PS", "Nylon", "PET"))
#c_back <- full %>% filter(polymer %in% c("PE", "PS", "wood"))
#ha <- full %>% filter(polymer %in% c("Nylon", "PET", "wood"))

# heatmap platics only ----------------------------------------------------

plastics_h <- ggplot(data = plastics, mapping = aes( x = fct_relevel(timepoint_days,  c("0","5", "10","30","45")),
                                              y = fct_relevel(Genus,rev), fill = Genus_rep_rel_abund)) +
  geom_tile() +
  scale_fill_gradientn(name = "relative abundance",colours =c("grey90","orange", "#380282") ,
                       values = rescale(x = c(0, 0.05, 0.8), from = c(0, 0.8)),limits=c(0,0.8))+ 
  facet_nested( ~ fct_relevel(polymer, c("PE", "PS", "Nylon", "PET"))+station+treatment,
                nest_line = element_line(linetype = 1, linewidth = 1)) +
  xlab(label = "incubation time (days)") +
  scale_y_discrete("Genus", breaks = breaks) +
  theme_minimal() +
  theme(strip.placement = "outside",
        plot.title = element_text(hjust = 0.5),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 12, face = "bold.italic", colour = condition),
        axis.text.x = element_text(size = 11, face = "bold", angle = 90, hjust = 1, vjust= 0.5),
        axis.title.x = element_text(face = "bold", vjust= -1,hjust = 0.5 ),
        strip.background = element_blank(),
        ggh4x.facet.nestline = element_line(colour = c("#890100")),
        strip.text.x = element_text(size= 12, face = "bold"),
        legend.title = element_text(face = "bold", vjust = 0.85),
        legend.position="bottom") 

# plastics_h

# ggexport(plast_h, "heatmap.pdf")
# heatmap genera open water -----------------------------------------------


# openw_h <- ggplot(data = os2, mapping = aes( x = fct_relevel(timepoint_days,  c("0","5", "10","30","45")),
#                                                 y = fct_relevel(Genus,rev), fill = Genus_rep_rel_abund)) +
#   
#   geom_tile() +
#   scale_fill_gradientn(name = "RA", colours = blu_grad, limits=c(0,0.75) )+ 
#   facet_nested( ~ polymer+treatment, nest_line = element_line(linetype = 1, linewidth = 1)) +
#   xlab(label = "incubation time (days)") +
#   scale_y_discrete("Genus", breaks = breaks) +
#   theme_minimal() +
#   theme(strip.placement = "outside",
#         plot.title = element_text(hjust = 0.5),
#         axis.title.y = element_blank(),
#         axis.text.y = element_text(face = "bold.italic", colour = condition),
#         axis.text.x = element_text(face = "bold"),
#         strip.background = element_blank(),
#         ggh4x.facet.nestline = element_line(colour = c("darkgoldenrod3")),
#         strip.text.x = element_text(face = "bold"),
#         legend.title = element_text(face = "bold"),
#         legend.position="right") 

# openw_h

# heatmap genera coastal station ------------------------------------------

# coast_h <- ggplot(data = coast2, mapping = aes( x = fct_relevel(timepoint_days,  c("0","5", "10","30","45")),
#                                           y = fct_relevel(Genus,rev), fill = Genus_rep_rel_abund)) +
# 
#     geom_tile() +
#     scale_fill_gradientn(name = "RA", colours = blu_grad, limits=c(0,0.75) )+ 
#     facet_nested( ~ polymer+treatment, nest_line = element_line(linetype = 1, linewidth = 1)) +
#     xlab(label = "incubation time (days)") +
#     scale_y_discrete("Genus", breaks = breaks) +
#     theme_minimal() +
#     theme(strip.placement = "outside",
#           plot.title = element_text(hjust = 0.5),
#           axis.title.y = element_blank(),
#           axis.text.y = element_text(face = "bold.italic", colour = condition),
#           axis.text.x = element_text(face = "bold"),
#           strip.background = element_blank(),
#           ggh4x.facet.nestline = element_line(colour = c("darkgoldenrod3")),
#           strip.text.x = element_text(face = "bold"),
#           legend.title = element_text(face = "bold"),
#           legend.position="right")

# coast_h
