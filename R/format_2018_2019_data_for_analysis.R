# subset and format 2018-2019 data for comparison with 2024-2025 points
  # selects the subset of points that correspond spatially to the 2024-2025 points
  # and formats them for distance sampling

format_2018_2019_data_for_analysis <- function(path = "data/Grassland_Birds_2018-2019_data.csv"){
  
# create a lookup table between standard 19 Kaufman & Skeet point count locations and nearest 2018-2019 survey points
survey.pt.lookup <- data.frame(pt = c("DUKE_01", "DUKE_01", 
                                      "DUKE_02", "DUKE_02",
                                      "DUKE_03", "DUKE_03",
                                      "DUKE_04", "DUKE_04",
                                      "DUKE_05", "DUKE_05",
                                      "DUKE_06", "DUKE_06",
                                      "DUKE_07", "DUKE_07",
                                      "DUKE_08", "DUKE_08",
                                      "DUKE_09", "DUKE_09",
                                      "DUKE_10", "DUKE_10",
                                      "DUKE_11", "DUKE_11",
                                      "DUKE_12", "DUKE_12",
                                      "DUKE_13", "DUKE_13",
                                      "DUKE_19", "DUKE_19",
                                      "DUKE_20", "DUKE_20",
                                      "DUKE_21", "DUKE_21",
                                      "DUKE_24", "DUKE_24",
                                      "DUKE_25", "DUKE_25",
                                      "DUKE_26", "DUKE_26"
),
pt18 = c("SR06", "SY04", 
         "SR24", "SY30",
         "SR19", "SY25",
         "SY46", "SR43",
         "SY41", "SR36",
         "SR57", "SY56",
         "KY36", "KR40",
         "KY31", "KR35",
         "KR50", "KY49",
         "KY23", "KR26",
         "KY17", "KR15",
         "KY05", "KR03",
         "KY08", "KR10",
         "KY42", "KR44",
         "KY25", "KR28",
         "KY47", "KR33",
         "SR08", "SY11",
         "SR22", "SY27",
         "SY44", "SR46")
)  

# read in the 2018-2019
u <- read.csv(path) %>%
  rename(pt18 = pt) 

# format data to match 2024-2025 (including timing and point locations)
u.fin <- u %>%
  # associate the 2024-2025 point name (pt) with the corresponding 2018-2019 point name (pt18)
  left_join(survey.pt.lookup, by = join_by(pt18)) %>%
  # remove surveys at points and survey periods not corresponding with the 2024-2025 surveys
  filter(!is.na(pt),
         per %in% 3:8) %>%
  # format time and date fields
  mutate(time = hour(hm(start)) + minute(hm(start))/60,
         date = mdy(date), 
         year = case_when(year(date) %in% 2018 ~ "Y18",
                          year(date) %in% 2019 ~ "Y19"),
         pt = paste0(year, pt)) %>%
  dplyr::rename(species = sp, observer = Collector) %>%
  # recode the period to match those used in 2024-2025
  mutate(period = case_when(month(date) %in% 5 ~ 1,
                            month(date) %in% 6 & day(date) < 16 ~ 2,
                            month(date) %in% 6 & day(date) > 15 ~ 3,
                            month(date) %in% 7 & day(date) < 2 ~ 3,
                            month(date) %in% 7 & day(date) > 1 ~ 4,
                            TRUE ~ NA),
         field = substr(pt18, 1,1),
         obsperpt = paste0(observer, period, pt),
         perpt = paste0(period, pt)) %>%
  select(observer, year, date, time, pt, field, 
         period, obsperpt, perpt, species, num, dist)

return(u.fin)

}