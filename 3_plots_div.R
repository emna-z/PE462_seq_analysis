library("scales")
library("RColorBrewer")
library("ggthemes")
library("Polychrome")
library("tidyverse")
library("Hmisc")
library("ggpubr")

# Kingdom -----------------------------------------------------------------

k <- read_csv("../Kingdom_PE462.csv")
k <- k %>% filter(station %in% c("C05","C13")) %>%
  filter(timepoint_days%nin% c("0")) %>% 
  mutate(station = str_replace_all(station,c("C05"= "open water station", "C13" = "coastal station"))) %>% 
  mutate(polymer_photo = str_replace_all(polymer_photo,c("_"= " "))) %>% 
  mutate(across(c(polymer,polymer_photo,timepoint_days,station),factor)) %>%
  distinct() 
k %>% head()

king_plot <- ggplot(k, aes(x=polymer_photo, y= Kingdom_rep_rel_abund, fill=Kingdom))+
  geom_bar(stat="identity", position="stack") +
  scale_fill_manual(values=c("#DDCC77","#88CCEE","#CC6677"))+
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
        legend.text = element_text(face = "bold"))

# ggexport(king,filename = "../kingdom_rel_abund_2.pdf")


# Phyla -------------------------------------------------------------------

 
phyla <- read_csv("../phyla_PE462.csv")
p <- phyla %>% filter(station %in% c("C05","C13")) %>%
  filter(timepoint_days%nin% c("0")) %>% 
  select(-Kingdom) %>% 
  select(- st_dev_Phylum_abund) %>% 
  mutate(station = str_replace_all(station,c("C05"= "open water station", "C13" = "Coastal station"))) %>% 
  mutate(polymer_photo = str_replace_all(polymer_photo,c("_"= " ")))%>% 
  mutate(across(c(polymer,polymer_photo,timepoint_days,station),factor)) %>% 
  mutate(Phylum, Phylum = if_else(Phyla_rep_rel_abund < 0.05, str_c("others <5%"), Phylum)) %>% 
  mutate(Phylum = fct_relevel(Phylum,"others <5%", after = Inf))

psum <- p %>% filter(Phylum %in% "others <5%") %>%
  group_by(detail) %>%
  mutate(Phyla_rep_rel_abund = sum(Phyla_rep_rel_abund)) %>%
  distinct()

p <- p %>% filter(Phylum %nin% "others <5%") %>%
bind_rows(psum) %>%
  mutate(Phylum = fct_relevel(Phylum,"others <5%", after = Inf))

rcartocolor::carto_pal(12,"Safe")

safe <- c ("#CC6677", "#DDCC77","#88CCEE" , "#117733", 
           "#332288", "#AA4499", "#44AA99", "#999933", 
           "#882255","#6699CC", "#661100", "#888888")
show_col(safe)

phyla_plot <- ggplot(p, aes(x=polymer_photo, y= Phyla_rep_rel_abund, fill=Phylum))+
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
        legend.text = element_text(face = "bold"))

# ggexport(ph,filename = "../phylum_rel_abund.pdf")

# Orders ------------------------------------------------------------------


or <- read_csv("../order_PE462.csv")
o <- or %>% filter(station %in% c("C05","C13")) %>%
  filter(timepoint_days%nin% c("0")) %>%
  select(-Kingdom) %>% 
  select(-st_dev_Order_abund) %>% 
  mutate(station = str_replace_all(station,c("C05"= "open water station", "C13" = "Coastal station"))) %>% 
  mutate(polymer_photo = str_replace_all(polymer_photo,c("_"= " ")))%>% 
  mutate(across(c(polymer,polymer_photo,timepoint_days,station),factor))%>%  
  mutate(Order, Order = if_else(Order_rep_rel_abund < 0.05, str_c("unassigned or <5%"), Order)) %>% 
  mutate(Order, Order = if_else(Order == "unassigned", str_c("unassigned or <5%"), Order)) %>% 
  distinct()
  
o5 <- o %>% filter(Order %in% "unassigned or <5%") %>%
  group_by(detail) %>% 
  mutate(Order_rep_rel_abund = sum(Order_rep_rel_abund)) %>%
  distinct()
  
of <- o %>% filter(Order %nin% "unassigned or <5%") %>% 
  bind_rows(o5) %>% 
  mutate(Order = fct_relevel(Order,"unassigned or <5%", after = Inf)) %>% 
  
  distinct()

t5 <- o %>% filter(Order %in% c("Enterobacterales", "Campylobacterales")) %>% 
  filter (station == "Coastal station") %>% 
  group_by(timepoint_days, Order) %>% 
  mutate(Order_time_rel_abund = mean(Order_rep_rel_abund)) %>% 
  mutate(Order_time_stdv = sd(Order_rep_rel_abund))

unique(o$Order)

sky <- Polychrome::sky.colors()
names(sky) <- NULL
swatch(sky)

polychrome_light <- light.colors()

light <- c("#FD3216", "#00FE35", "#6A76FC", "#FED4C4", "#FE00CE", 
"#0DF9FF", "#F6F926", "#FF9616", "#479B55", "#EEA6FB", "#DC587D", 
"#D626FF", "#6E899C", "#00B5F7", "#B68E00", "#C9FBE5", "#FF0092", 
"#22FFA7", "#E3EE9E", "#86CE00", "#BC7196", "#7E7DCD", "#FC6955", "#E48F72")

light_l <- light[c()]
swatch(light)

colors_M2 <- c( "#3d405b","#ce796b", 
                "#25a18e", "#BC412B", "#95D9DA", "#B10F2E", "#0E273C",
                "#E3FDD8", "#353535", "#e7ad99", "#0F8B8D", "#7ae582",
                "#F2AF29",   "#94d2bd", "#606c38","#772e25",
                "#0047E0", "#344e41","#6c584c", "#5f0f40", "#D7F171", "#c89f9c" )
swatch(colors_M2)

p_o <- ggplot(of, aes(x=polymer_photo, y= Order_rep_rel_abund, fill=Order))+
  geom_bar(stat="identity", position="stack")+ scale_fill_manual(values = light)+
  scale_y_continuous(labels=percent)+
  theme_minimal()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab("")+ylab("Relative Abundance")+
  facet_grid (fct_relevel(timepoint_days, c('5',"10", "30", "45"))~station)+ 
  theme(strip.text.x = element_text(size=12, face="bold"),strip.text.y = element_text(size=10, face="bold", angle = 0),
        legend.title = element_text(face = "bold"), 
        legend.title.align = 0.4,
        legend.text = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        axis.text.x = element_text(face = "bold"))
p_o <- p_o + guides(fill=guide_legend(ncol =1))

# ggpubr::ggexport(p_o, "../orders.pdf")

