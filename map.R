# install.packages(c("cowplot", "googleway", "ggplot2", "ggrepel", 
# "ggspatial", "libwgeom", "sf", "rnaturalearth", "rnaturalearthdata"))


# World -------------------------------------------------------------------

library("ggplot2")

library("rnaturalearth")
library("rnaturalearthdata")

world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)


library("sf")
world_points<- st_centroid(world)
world_points <- cbind(world, st_coordinates(st_centroid(world$geometry)))



if(!file.exists("maps/ne_10m_rivers_europe.shp")){
  ne_download(scale = 10, type = "rivers_lake_centerlines", category = "physical", 
              destdir = "maps/", load = FALSE) # major rivers
  ne_download(scale = 10, type = "lakes", category = "physical", 
              destdir = "maps/", load = FALSE) # major lakes
}


rivers <- ne_load(scale = 10, type = "rivers_lake_centerlines", destdir = "maps", returnclass = "sf")
lakes <- ne_load(scale = 10, type = "lakes", destdir = "maps", returnclass = "sf")


# North sea map -----------------------------------------------------------

library("ggspatial")

ggplot(data = world) +
  geom_sf(fill = "#FAE5D3") +
  geom_sf(data = rivers, colour = scales::alpha("#3498DB", 0.3),linewidth = 0.2) + 
  geom_sf(data = lakes, fill = "#D6EAF8") + 
  xlab("Longitude") + ylab("Latitude") +
  theme(axis.title.x = element_text(face = "bold", size= 15),
        axis.title.y = element_text(face = "bold", size = 15),
        axis.text.x = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 14, face = "bold"),
        panel.grid.major = element_line(colour = NA),
        panel.grid.minor = element_line(colour = NA),
        panel.background = element_rect(fill = "#D6EAF8"))+
  annotation_scale(location = "tr", width_hint = 0.5 ) +
  annotation_north_arrow(location = "tr", which_north = "true", 
                         pad_x = unit(0.75, "in"), pad_y = unit(0.5, "in"),
                         style = north_arrow_fancy_orienteering) +
  coord_sf(xlim = c(-3, 9), ylim = c(50, 57), expand = FALSE) +
  geom_text(data= world_points,aes(x=X, y=Y, label=name), size = 5,
            color = "#943126", fontface = "bold", check_overlap = FALSE) +
  annotate(geom = "text", x = 2, y = 56, label = "North Sea (NS)", 
           fontface = "bold.italic", color = "#000080", size =10) +
  geom_point( aes(x = 4, y = 52.085), color = "#D32F2F", size = 4, shape= 18) +
  geom_text(label="coastal NS station", x = 4, y = 52.5, color = "#D32F2F", check_overlap = T, size = 7)+
  geom_point( aes(x = 3.7, y = 54.7), color = "#D32F2F", size = 4, shape= 18) +
  geom_text(label="open NS water station", x = 3.7, y = 54.9, color = "#D32F2F", check_overlap = T, size = 7)
  

# saving ------------------------------------------------------------------


# ggoceanMaps -------------------------------------------------------------
# 
# remotes::install_github("MikkoVihtakari/ggOceanMapsData")
# install.packages("ggOceanMaps")
# 
# library(ggOceanMaps)
# #limits are given longitude min/max, latitude min/max
# basemap(limits = c(-3, 9, 50, 57),
#         bathymetry = TRUE,
#         glaciers = FALSE)
