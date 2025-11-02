# script to run through line-by-line to help proof grassland bird data

library(dplyr)
library(lubridate)

# read data
# b <- read.csv("/Users/mikea/Documents/research/duke_2025/data/Grassland_Birds_2025_data_proofing_20250919b.csv")
b <- read.csv("~/Documents/research/duke_2025/data/20250520_grassland_birds_2025_Allen/Grassland_Birds_2025_data.csv")
pt_data <- read.csv("/Users/mikea/Documents/research/duke_2025/data/pt_coords_2024.csv")

# Project_Name	latitude	longitude	Time_Date_Stamp	Collector	pt	period	recorder	
# year	date	temp.s	temp.e	wind.s	wind.e	sky.s	sky.e	field	ID	point	start	
# species	num	minute	distance	az	detect	sex	voc	breed	notes

  ###########################
  ### quality checks 2025 ###
  ###########################

# number of survey per point
b %>%
  select(point, date) %>%
  distinct() %>%
  group_by(point) %>%
  tally()

# examine all species names for spelling
sort(unique(b$species))
b <- b %>%
  mutate(species = case_when(species %in% "UNKN" ~ "UNBI",
                             TRUE ~ species))
sort(unique(b$species))

# check 4 surveys per point
b %>%
  select(point, date) %>%
  distinct() %>%
  group_by(point) %>%
  tally()

# check reasonableness of other values
sort(unique(b$Collector))
sort(unique(b$period))
sort(unique(b$recorder))
sort(unique(b$year))
sort(unique(b$date))
sort(unique(b$temp.s))
sort(unique(b$temp.e))
sort(unique(b$wind.s))
sort(unique(b$wind.e))
sort(unique(b$sky.s))
sort(unique(b$sky.e))
sort(unique(b$field))
sort(unique(b$ID))
sort(unique(b$point)); length(unique(b$point))
sort(unique(b$start))
sort(unique(b$num))
sort(unique(b$minute))
sort(unique(b$distance))
sort(unique(b$az))
sort(unique(b$detect))
sort(unique(b$sex))
table(b$voc) # fixed 3 voc=O
table(b$breed) # fixed 2 "chuck" entries

# check consistency of voc and sex fields
  # e.g., all voc=S or SC are also sex=M (i.e., singing = male)
sex_check <- b %>%
  filter(voc %in% c("S", "SC"),
         !sex %in% "M")

# check consistency of detect and voc fields
detect_voc_check <- b %>%
  filter(detect %in% c("S"),
         !voc %in% "N")

detect_voc_check2 <- b %>%
  filter(detect %in% c("SH","H"),
         voc %in% "N")

# read in proofed 2025 data and format columns
gl2025 <- b %>%
  mutate(pt = case_when(nchar(point) %in% 1 ~ paste0("DUKE_0", point),
                   TRUE ~ paste0("DUKE_", point))) %>%
  select(-latitude, -longitude) %>%
  left_join(pt_data, by = join_by(pt)) %>%
  mutate(
    date_parsed = mdy(date),  # or another format depending on the input
    datetime_str = paste(date_parsed, start),
    Time_Date_Stamp = ymd_hm(datetime_str),
    Project_Name = "Grassland Birds 2025"
  ) %>%
  select(-date_parsed, -datetime_str) %>%
  select(Project_Name, latitude, longitude, 
         Time_Date_Stamp, Collector,
          pt, period:notes)

table(gl2025$minute)
gl2025_fix_min <- gl2025 %>%
  mutate(notes = case_when(minute %in% 0 & notes %in% "" ~ 
                             paste0(notes, "[minute=0 on datasheet; correctd to minute=1]"),
                           minute %in% 0 & !notes %in% "" ~ 
                             paste0(notes, " [minute=0 on datasheet; correctd to minute=1]"),
                           TRUE ~ notes),
         minute = case_when(minute %in% 0 ~ 1,
                            TRUE ~ minute))
table(gl2025_fix_min$minute)

write.csv(gl2025_fix_min,
          "~/Documents/research/duke_2025/data/20250520_grassland_birds_2025_Allen/Grassland_Birds_2025_data.csv",
          row.names = F)

# 2024 fields
# Project_Name	latitude	longitude	Time_Date_Stamp	Collector	pt	period	recorder	
# year	date	temp.s	temp.e	wind.s	wind.e	sky.s	sky.e	field	ID	point	start	species	
# num	minute	distance	az	detect	sex	voc	breed	notes
