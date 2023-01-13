######libraries#####
library(tidyverse)
library(Hmisc)

# import of final file created in script 1_data_import.R ------------------
# t3 <- read_csv("../tidyPE462_rel_abund_calc.csv")


# Kingdom -----------------------------------------------------------------
Kingdom <- t3  %>%  
  select(station,timepoint_days,treatment,polymer,pol_photo_station, polymer_station,polymer_photo,
          material, detail,Kingdom, Kingdom_rep_rel_abund,Kingdom_st_dev_abund_samples)%>% 
  distinct() 
# write_csv(Kingdom, "../Kingdom_PE462.csv")

# Phylum ------------------------------------------------------------------
Phyla <- t3 %>%
  select(station,timepoint_days,treatment,polymer,pol_photo_station, polymer_station,
                       polymer_photo, material, detail, Kingdom,Phylum, Phyla_rep_rel_abund ,st_dev_Phylum_abund)%>% 
                       distinct() 
# write_csv(Phyla, "../phyla_PE462.csv")


# Class -------------------------------------------------------------------
Class <-t3 %>% 
  select(station,timepoint_days,treatment,polymer,pol_photo_station, polymer_station,
                      polymer_photo, material, detail,Kingdom,
                      Class, Class_rep_rel_abund, st_dev_Class_abund )%>% 
distinct()  
# write_csv(Class, "../class_PE462.csv")


# Order -------------------------------------------------------------------
Order <- t3 %>% 
  select(station,timepoint_days,treatment,polymer,pol_photo_station, polymer_station,
                       polymer_photo, material, detail,Kingdom,
                       Order, Order_rep_rel_abund, st_dev_Order_abund )%>% 
  distinct() 
# write_csv(Order, "../order_PE462.csv")



# Family ------------------------------------------------------------------
Family <- t3 %>% 
  select(station,timepoint_days,treatment,polymer, polymer_photo, material, detail,
                        Kingdom, Family, Family_rep_rel_abund, st_dev_Family_abund )%>% 
  distinct()
# write_csv(Family, "../family_PE462.csv")

# Genus -------------------------------------------------------------------
Genus <- t3%>% 
  select(station,timepoint_days,treatment,polymer,pol_photo_station, polymer_station,
                      polymer_photo, material, detail, Kingdom,Genus, Genus_rep_rel_abund, st_dev_Genus_abund )%>%
  distinct()
# write_csv(Genus, "../genus_PE462.csv")


# Species -----------------------------------------------------------------
Species <- t3 %>% 
  select(station,timepoint_days,treatment,polymer,pol_photo_station, polymer_station,
                         polymer_photo, material, detail,Kingdom,
                         Species, Species_rep_rel_abund, st_dev_Species_abund)%>%
  distinct()
# write_csv(Species, "../species_PE462.csv")
