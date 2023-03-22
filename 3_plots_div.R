##################################################
## Project: PE462
## Script purpose: taxonomic ranks relative 
## abundace barplots
## Date: 14-01-2023
## Author: E. Zeghal
##################################################

library("scales")
library("RColorBrewer")
library("ggthemes")
library("Polychrome")
library("tidyverse")
library("Hmisc")
library("ggpubr")


# palettes ----------------------------------------------------------------

rcartocolor::carto_pal(12,"Safe")
safe <- c ("#CC6677", "#DDCC77","#88CCEE" , "#117733", 
                    "#332288", "#44AA99", "#661100","#882255", 
                    "#6699CC", "#999933", "#888888","#AA4499" )
# show_col(safe)

light_poly <- c("#FD3216", "#00FE35", "#6A76FC", "#FED4C4", "#FE00CE", 
                         "#0DF9FF", "#F6F926", "#FF9616", "#479B55", "#EEA6FB", "#DC587D", 
                         "#D626FF", "#6E899C", "#00B5F7", "#B68E00", "#C9FBE5", "#FF0092", 
                         "#22FFA7", "#E3EE9E", "#86CE00", "#BC7196", "#7E7DCD", "#FC6955", "#E48F72")
                         
colors_M2 <- c( "#3d405b","#ce796b", 
                         "#25a18e", "#BC412B", "#95D9DA", "#B10F2E", "#0E273C",
                         "#E3FDD8", "#353535", "#e7ad99", "#0F8B8D", "#7ae582",
                         "#F2AF29",   "#94d2bd", "#606c38","#772e25",
                         "#0047E0", "#344e41","#6c584c", "#5f0f40", "#D7F171", "#c89f9c" )
# swatch(colors_M2)

kelly1 <-  c("#be0032", "#f3c300",   "#875692",   "#f38400" ,  "#a1caf1" ,   "#c2b280" ,  
                      "#848482" ,  "#008856"  , "#e68fac"  ,"#0067a5" ,  "#f99379"  ,
                      "#604e97" ,  "#f6a600"  , "#b3446c"  , "#dcd300"  ,   "#25a18e" , 
                      "#882d17" , "#e25822"  , "#2b3d26", "lightgrey"  , "#654522", "#8db600")
                      

# Kingdom -----------------------------------------------------------------

k <- read_csv("../Kingdom_PE462.csv")
k <- k %>% 
  filter(station %in% c("C05","C13")) %>%
  filter(timepoint_days%nin% c("0")) %>% 
  mutate(station = str_replace_all(station,c("C05"= "open NS water", "C13" = "coastal NS water"))) %>% 
  mutate(polymer_photo = str_replace_all(polymer_photo,c("_"= " "))) %>% 
  mutate(across(c(polymer,polymer_photo,timepoint_days,station),factor)) %>%
  distinct() 
k %>% head()

# k %>% 
#   group_by(Kingdom) %>% 
#   summarise(mean = round(mean(Kingdom_rep_rel_abund)*100, digits = 1), 
#             stdev= round(sd(Kingdom_rep_rel_abund)*100, digits = 1))


# kingdoms RA barplot
king_plot <- ggplot(k, aes(x=fct_relevel(polymer_photo, c("PE no UV", "PE UV", "PS no UV", "PS UV","Nylon no UV", "Nylon UV", "PET no UV", "PET UV", "wood")),
                           y= Kingdom_rep_rel_abund, fill=Kingdom))+
  geom_bar(stat="identity", position="stack") +
  scale_fill_manual(values=safe)+
  scale_y_continuous(labels=percent)+
  theme_minimal()+
  xlab("")+
  ylab("Relative Abundance")+
  facet_grid (fct_relevel(timepoint_days, c('5',"10", "30", "45"))~station)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, face = "bold"),
        axis.title.y = element_text(face = "bold") ,
        strip.text.x = element_text(size=12, face="bold"),
        strip.text.y = element_text(size=10, face="bold", angle = 0),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(face = "bold"),
        panel.grid.major = element_line(colour = NA),
        panel.grid.minor = element_line(colour = NA))

# ggexport(king_plot,filename = "../kingdom_rel_abund_3.pdf")


# Phyla -------------------------------------------------------------------

 
phyla <- read_csv("../phyla_PE462.csv")

## transform for writing ---------------------------------------------------

# p1 <- phyla %>% 
#   filter(station %in% c("C05","C13")) %>%
#   filter(timepoint_days%nin% c("0")) %>%
#   select(-Kingdom) %>% 
#   select(-st_dev_Phylum_abund) %>% 
#   mutate(station = str_replace_all(station,c("C05"= "open NS water", "C13" = "coastal NS water"))) %>% 
#   mutate_if(is.character,as.factor) %>% 
#   group_by(station, timepoint_days, Phylum) %>% 
#   summarise(mean = round(mean(Phyla_rep_rel_abund)*100, digits = 1), 
#             stdev= round(sd(Phyla_rep_rel_abund)*100, digits = 1))
# 
# bact_phyla <- unique((phyla %>% filter(Kingdom == "Bacteria"))$Phylum)

p2 <- phyla %>%
  filter(station %in% c("C05","C13")) %>%
  filter(timepoint_days%nin% c("0")) %>%
  select(-Kingdom) %>%
  select(-st_dev_Phylum_abund) %>%
  mutate(station = str_replace_all(station,c("C05"= "open NS water", "C13" = "coastal NS water"))) %>%
  mutate_if(is.character,as.factor) %>%
  group_by(station, timepoint_days,polymer, Phylum) %>%
  summarise(mean = round(mean(Phyla_rep_rel_abund)*100, digits = 1),
            stdev= round(sd(Phyla_rep_rel_abund)*100, digits = 1))
p3 <- phyla %>%
  filter(station %in% c("C05","C13")) %>%
  filter(timepoint_days%nin% c("0")) %>%
  select(-Kingdom) %>%
  select(-st_dev_Phylum_abund) %>%
  mutate(station = str_replace_all(station,c("C05"= "open NS water", "C13" = "coastal NS water"))) %>%
  mutate_if(is.character,as.factor) %>%
  group_by( timepoint_days,polymer, Phylum) %>%
  summarise(mean = round(mean(Phyla_rep_rel_abund)*100, digits = 1),
            stdev= round(sd(Phyla_rep_rel_abund)*100, digits = 1))

## transform for barplot ---------------------------------------------------


p <- phyla %>% 
  filter(station %in% c("C05","C13")) %>%
  filter(timepoint_days%nin% c("0")) %>%
  mutate(Phylum, Phylum = if_else(Kingdom %nin% c("Bacteria"), str_c("others <5%"), Phylum)) %>% 
  select(-Kingdom) %>% 
  select(- st_dev_Phylum_abund) %>% 
  mutate(station = str_replace_all(station,c("C05"= "open NS water", "C13" = "coastal NS water"))) %>% 
  mutate(polymer_photo = str_replace_all(polymer_photo,c("_"= " ")))%>% 
  mutate(across(c(polymer,polymer_photo,timepoint_days,station),factor)) %>% 
  mutate(Phylum, Phylum = if_else(Phyla_rep_rel_abund < 0.05, str_c("others <5%"), Phylum)) %>% 
  mutate(Phylum = fct_relevel(Phylum,"others <5%", after = Inf))

# the multiple copies of phyla that are now labelled "others <5%" produce a distorted effect when plotting
# to avoid that we'll sum the RA of all "others <5%" for each condition and replace the multiple 
# entries by a single one

psum <- p %>% 
  filter(Phylum %in% "others <5%") %>%
  group_by(detail) %>%
  mutate(Phyla_rep_rel_abund = sum(Phyla_rep_rel_abund)) %>%
  distinct() #psum is the sum of RA of others <5% per condition

p <- p %>% 
  filter(Phylum %nin% "others <5%") %>%
  bind_rows(psum) %>% # the multiple others <5% are now replaced by a single one per condition
  mutate(Phylum = fct_relevel(Phylum,"others <5%", after = Inf)) 


# phylum RA barplot
phyla_plot <- ggplot(p, aes(x=fct_relevel(polymer_photo, c("PE no UV", "PE UV", "PS no UV", "PS UV","Nylon no UV", "Nylon UV", "PET no UV", "PET UV", "wood")),
                                          y= Phyla_rep_rel_abund, fill=Phylum))+
  geom_bar(stat="identity", position="stack")+
  scale_fill_manual(values =  safe)+
  theme_minimal()+ 
  scale_y_continuous(labels=percent)+
  theme_minimal()+
  xlab("")+
  ylab("Relative Abundance")+
  facet_grid (fct_relevel(timepoint_days, c('5',"10", "30", "45"))~station)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, face = "bold"),
        axis.title.y = element_text(face = "bold") ,
        strip.text.x = element_text(size=12, face="bold"),
        strip.text.y = element_text(size=10, face="bold", angle = 0),
        legend.title = element_text(face = "bold"),
        legend.title.align = 0.4,
        legend.text = element_text(face = "bold"),
        panel.grid.major = element_line(colour = NA),
        panel.grid.minor = element_line(colour = NA))

# ggexport(ph,filename = "../phylum_rel_abund.pdf")


# Orders ------------------------------------------------------------------

or <- read_csv("../order_PE462.csv")


## transform for writing ---------------------------------------------------

o1 <- or %>% 
  filter(station %in% c("C05","C13")) %>%
  filter(timepoint_days%nin% c("0")) %>%
  select(-Kingdom) %>% 
  select(-st_dev_Order_abund) %>% 
  mutate(station = str_replace_all(station,c("C05"= "open NS water", "C13" = "coastal NS water"))) %>% 
  mutate_if(is.character,as.factor) %>% 
  group_by(station, timepoint_days, Order) %>% 
  summarise(mean = round(mean(Order_rep_rel_abund)*100, digits = 1), 
            stdev= round(sd(Order_rep_rel_abund)*100, digits = 1))

o2 <- or %>% 
  filter(station %in% c("C05","C13")) %>%
  filter(timepoint_days%nin% c("0")) %>%
  select(-Kingdom) %>% 
  select(-st_dev_Order_abund) %>% 
  mutate(station = str_replace_all(station,c("C05"= "open NS water", "C13" = "coastal NS water"))) %>% 
  mutate_if(is.character,as.factor) %>% 
  group_by( timepoint_days, polymer,Order) %>% 
  summarise(mean = round(mean(Order_rep_rel_abund)*100, digits = 1), 
            stdev= round(sd(Order_rep_rel_abund)*100, digits = 1))

## transform for barplots -------------------------------------------------


o <- or %>% 
  filter(station %in% c("C05","C13")) %>%
  filter(timepoint_days%nin% c("0")) %>%
  select(-Kingdom) %>% 
  select(-st_dev_Order_abund) %>% 
  mutate(station = str_replace_all(station,c("C05"= "open NS water", "C13" = "coastal NS water"))) %>% 
  mutate(polymer_photo = str_replace_all(polymer_photo,c("_"= " ")))%>% 
  mutate(across(c(polymer,polymer_photo,timepoint_days,station),factor))%>%  
  mutate(Order, Order = if_else(Order_rep_rel_abund < 0.05, str_c("unassigned or <5%"), Order)) %>% 
  mutate(Order, Order = if_else(Order == "unassigned", str_c("unassigned or <5%"), Order)) %>% 
  distinct()

# Like phyla, the multiple copies of orders that are now labelled "unassigned or <5%" produce a 
# distorted effect when plotting to avoid that we'll sum the RA of all "others <5%" for each 
# condition and replace the multiple entries by a single one
  
o5 <- o %>% 
  filter(Order %in% "unassigned or <5%") %>%
  group_by(detail) %>% 
  mutate(Order_rep_rel_abund = sum(Order_rep_rel_abund)) %>%
  distinct()
  
of <- o %>% filter(Order %nin% "unassigned or <5%") %>% 
  bind_rows(o5) %>% 
  mutate(Order = fct_relevel(Order,"unassigned or <5%", after = Inf)) %>% 
  distinct()

# orders RA barplot
p_o <- ggplot(of, aes(x=fct_relevel(polymer_photo, c("PE no UV", "PE UV", "PS no UV", "PS UV","Nylon no UV", "Nylon UV", "PET no UV", "PET UV", "wood")), 
                      y= Order_rep_rel_abund, fill=Order))+
  geom_bar(stat="identity", position="stack")+ scale_fill_manual(values = kelly1)+
  scale_y_continuous(labels=percent)+
  theme_minimal()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab("")+ylab("Relative Abundance")+
  facet_grid (fct_relevel(timepoint_days, c('5',"10", "30", "45"))~station)+ 
  theme(strip.text.x = element_text(size=12, face="bold"),strip.text.y = element_text(size=10, face="bold", angle = 0),
        legend.title = element_text(face = "bold"), 
        legend.title.align = 0.4,
        legend.text = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        axis.text.x = element_text(face = "bold"),
        panel.grid.major = element_line(colour = NA),
        panel.grid.minor = element_line(colour = NA))
order_plot <- p_o + guides(fill=guide_legend(ncol =1))
order_plot
# ggpubr::ggexport(order_plot, "../orders.pdf")

