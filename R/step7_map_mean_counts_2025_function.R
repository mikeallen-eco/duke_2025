# function to plot average raw counts over time +/- 2 SE

library(terra)
library(sf)
library(maptiles)
library(ggplot2)
library(ggspatial)

map_mean_counts <- function(data = d, 
                             alpha = "GRSP", 
                             pts = "data/pt_coords_2010_2025.csv",
                             fields_file = "data/Kaufman_Skeet.geojson",
                             year_to_map = 2025){
  
  if(alpha == "GRSP") spname <- "Grasshopper Sparrow"
  if(alpha == "EAME") spname <- "Eastern Meadowlark"
  if(alpha == "RWBL") spname <- "Red-winged Blackbird"
  if(alpha == "BOBO") spname <- "Bobolink"
  
  # load point coordinates
  points_df <- read.csv(pts)
  
  # get complete list of points to fill in zero samples
  allpts <- data %>% 
    select(field, year, pt) %>% 
    distinct()
  
  # format raw count data
  site_data <- data %>% 
    filter(species %in% alpha) %>%
    filter(dist < 101 | is.na(dist)) %>% 
    group_by(field, year, pt) %>% 
    tally() %>%
    right_join(allpts %>% select(year, field, pt) %>% distinct(), 
               by = join_by(pt, year, field)) %>%
    replace_na(list(n = 0)) %>%
    arrange(pt) %>%
    mutate(dens = n/(4*3.1415),
           pt = substr(pt, 4, 10)) %>%
    filter(year %in% paste0("Y", substr(year_to_map,3,4))) %>%
    left_join(points_df,
              by = join_by(pt)) %>%
    st_as_sf(coords = c("longitude", "latitude"), crs = 4326)
  
    # load fields polygons
    fields <- read_sf(fields_file)

    # Bounding box with buffer
    site_bbox <- st_bbox(site_data)
    buffer <- 0.004
    
    minx <- round(min(st_coordinates(site_data)[,1]),3) - buffer
    maxx <- round(max(st_coordinates(site_data)[,1]),3) + buffer
    xbreaks <- seq(minx, maxx, by = 0.005)
  
map <- ggplot() +
    geom_sf(data = fields) +
    geom_sf(data = site_data, aes(color = dens), shape = 16, size = 5, stroke = 1.2) +
    coord_sf(
        xlim = c(site_bbox["xmin"] - buffer, site_bbox["xmax"] + buffer),
        ylim = c(site_bbox["ymin"] - buffer, site_bbox["ymax"] + buffer),
        expand = FALSE
      ) +
      annotation_scale(location = "br", width_hint = 0.3,
                       pad_y = unit(0.6, "cm"),
                       pad_x = unit(0.6, "cm"),
                       text_col = "black"
      ) +
      annotation_north_arrow(
        location = "tl", #which_north = "true",
        style = north_arrow_fancy_orienteering,
        height = unit(0.6, "cm"), width = unit(0.6, "cm"),
        pad_y = unit(0.6, "cm"),
        pad_x = unit(0.6, "cm")
        # text_col = "white"
      ) +
      scale_color_viridis_c(option = "viridis") + # , end = 0.85
      scale_x_continuous(
        labels = scales::label_number(accuracy = 0.001),
        breaks = xbreaks
      ) +
      scale_y_continuous(
        labels = scales::label_number(accuracy = 0.001)
      ) +
  coord_sf(crs = 26918) +
      labs(x = "", y = "", color = "Mean\ncount\nper ha") +
      theme_bw() +
      ggtitle(spname) +
      theme(
        # legend.position = "none",
        plot.title = element_text(size = 14, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        axis.text.y = element_text(size = 9)
      )

if(alpha == "GRSP"){
map <- map + geom_text(
  data = site_data,
  aes(label = pt, geometry = geometry),
  stat = "sf_coordinates",
  color = "black",
  size = 1.5,
  nudge_y = 0.001
)}
  
  return(map)
  
}
  
