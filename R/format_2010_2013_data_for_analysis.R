# function to read in raw 2010-2013 Duke Farms grassland bird data and format for distance sampling

library(dplyr)
library(lubridate)

format_2010_2013_data_for_analysis <- function(path = "data/Duke_grassland_birds_2010_2013.csv"){

# read in and combine 2010, 2012, and 2013 data
# read in 2010-2013 data, exclude points outside of Kaufman and Skeet
data <- read.csv(path) %>%
  filter(!PointId %in% c("DUKE_27", "DUKE_28", "DUKE_29",
                         "DUKE_30", "DUKE_31", "DUKE_16",
                         "DUKE_18", "DUKE_22", "DUKE_23",
                         "DUKE_14", "DUKE_15", "DUKE_17")) %>%
  # format period field
  mutate(period = case_when(month(VisitDate) %in% 5 ~ 1,
                            month(VisitDate) %in% 6 & day(VisitDate) < 16 ~ 2,
                            month(VisitDate) %in% 6 & day(VisitDate) > 15 ~ 3,
                            month(VisitDate) %in% 7 & day(VisitDate) < 2 ~ 3,
                            month(VisitDate) %in% 7 & day(VisitDate) > 1 ~ 4,
                            TRUE ~ NA),
         # format observer field
         observer = case_when(UserId %in% "MikeA" ~ "Mike Allen",
                              UserId %in% "kristin" ~ "Kristin M",
                              UserId %in% "laurastern" ~ "Laura Stern",
                              TRUE ~ NA),
         # format field field
         field = case_when(PointId %in% c("DUKE_01", "DUKE_02", "DUKE_03",
                                          "DUKE_04", "DUKE_05", "DUKE_06",
                                          "DUKE_24", "DUKE_25", "DUKE_26") ~ "S",
                           TRUE ~ "K"),
         # format sp field
         sp = case_when(sp %in% "Bobolink" ~ "BOBO",
                        sp %in% "Grasshopper Sparrow" ~ "GRSP",
                        sp %in% "Eastern Meadowlark" ~ "EAME",
                        sp %in% "Red-winged Blackbird" ~ "RWBL",
                        sp %in% "American Kestrel" ~ "AMKE",
                        TRUE ~ "Other")) %>%
  # convert to decimal time h.hh
  mutate(time = as.numeric(substr(VisitTime, 1,1)) + 
           as.numeric(substr(VisitTime, 2,3))/60) %>%
  mutate(pt = paste0(year, PointId)) %>%
  mutate(obsperpt = paste0(observer, period, pt),
         perpt = paste0(period, pt)) %>%
  select(observer, year, date = VisitDate, time, pt, field, period, 
         obsperpt, perpt, species = sp, num = n, distcat = dist)

return(data)

}
