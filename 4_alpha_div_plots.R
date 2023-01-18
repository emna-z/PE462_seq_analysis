library(tidyverse)
library(Hmisc)
library(ggpubr)
library(rcartocolor)


# Load data ---------------------------------------------------------------


alpha_tab <- read_csv("../alpha_div_indexes_microbiome_package_with_unassigned_phyla.csv")

alpha_tab <- alpha_tab %>%
  filter(station %in% c("C05", "C13")) %>% filter(material %nin% c("Glass_fiber")) %>% 
  filter(timepoint_days %nin% c("15")) 

alpha_tab$polymer <- if_else(alpha_tab$material=="wood", str_c(alpha_tab$material), alpha_tab$polymer)



# Simpson eveness ---------------------------------------------------------


div_simp <- alpha_tab %>%
  mutate_if(is.character, as.factor) %>% 
  select(detail, treatment, material,polymer,
         timepoint_days ,station, diversity_gini_simpson )%>% 
  group_by(treatment,polymer,station,timepoint_days) %>% 
  summarise(mean=mean(diversity_gini_simpson),
            n = n(),
            sd=sd(diversity_gini_simpson)) %>% 
  mutate(station = recode(station,
                   "C05" = "open water station",
                   "C13" = "coastal station")) %>% 
  mutate(timepoint_days=fct_relevel(timepoint_days, "5","10", "30", "45"))
            


# Shannon diversity -------------------------------------------------------


div_shan <- alpha_tab %>%
  mutate_if(is.character, as.factor) %>% 
  select(detail, treatment, material,polymer,
         timepoint_days ,station, diversity_shannon )%>% 
  group_by(treatment,polymer,station,timepoint_days) %>% 
  summarise(mean=mean(diversity_shannon),
            n = n(),
            sd=sd(diversity_shannon)) %>% 
  mutate(station = recode(station,
                          "C05" = "open water station",
                          "C13" = "coastal station")) %>% 
  mutate(timepoint_days=fct_relevel(timepoint_days, "5","10", "30", "45"))


# Graphs ------------------------------------------------------------------
# display_carto_all(colorblind_friendly = TRUE)
# rcartocolor::display_carto_pal(name = "Safe", n = 12)

## eveness simpson --------------------------------------------------------

simpson <-ggplot(div_simp,
                 aes(x = reorder(timepoint_days, as.numeric(timepoint_days)), 
                     y = mean, color=treatment)) +
  geom_point(size = 3, position=position_dodge(width=0.5))+
  geom_errorbar(aes(ymin=mean-sd,
                    ymax=mean+sd),
                width = .2,
                position=position_dodge(width=0.5)) +
  scale_color_manual(values=c("#661100","#88CCEE" , "#f3c300" ),
                     labels = c("no" = "none",
                                "no_UV" = " not UV treated",
                                "UV" = "UV treated")) +
  theme_minimal()+
  theme(strip.text.x = element_text(colour = "black", face = "bold", size = 12),
        strip.text.y = element_text(colour = "black", face = "bold", size = 10),
        strip.background = element_rect(color = NA),
        axis.text.x = element_text(colour = "black", face = "bold"),
        axis.text.y = element_text(colour = "black", face = "bold"),
        axis.title.x = element_text(colour = "black", face = "bold"),
        axis.title.y = element_text(colour = "black", face = "bold"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y =  element_line(colour = "grey80", linewidth = 0.08, linetype = "dotted"),
        panel.grid.minor.y = element_blank(),
        legend.position = "top",
        legend.title = element_text(face = "bold"),
        legend.text = element_text(face = "bold"))+
  facet_grid(polymer ~ station)+
  labs(x="incubation time (days)", 
       y = "Gini-Simpson diversity index")
# simpson


# ggexport(simpson,filename = "./plots/simpson_index.pdf")

## Shannon diversity -----------------------------------------------------


shannon <-ggplot(div_shan,
                 aes(x = reorder(timepoint_days, as.numeric(timepoint_days)), 
                     y = mean, color=treatment)) +
  geom_point(size = 3, position=position_dodge(width=0.5))+
  geom_errorbar(aes(ymin=mean-sd,
                    ymax=mean+sd),
                width = .2,
                position=position_dodge(width=0.5)) +
  scale_color_manual(values=c("#661100","#88CCEE" , "#f3c300"),
                     labels = c("no" = "none",
                                "no_UV" = " not UV treated",
                                "UV" = "UV treated")) +
  theme_minimal()+
  theme(strip.text.x = element_text(colour = "black", face = "bold", size = 12),
        strip.text.y = element_text(colour = "black", face = "bold", size = 10),
        strip.background = element_rect(color = NA),
        axis.text.x = element_text(colour = "black", face = "bold"),
        axis.text.y = element_text(colour = "black", face = "bold"),
        axis.title.x = element_text(colour = "black", face = "bold"),
        axis.title.y = element_text(colour = "black", face = "bold"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y =  element_line(colour = "grey80", linewidth = 0.08, linetype = "dotted"),
        panel.grid.minor.y = element_blank(),
        legend.position = "top",
        legend.title = element_text(face = "bold"),
        legend.text = element_text(face = "bold"))+
  facet_grid(polymer ~ station)+
  labs(x="incubation time (days)", 
       y = "Shannon diversity index")
shannon

# ggexport(shannon,filename = "./plots/shannon_index.pdf")

