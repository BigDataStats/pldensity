# Library ========================================================
library("pldensity") # Our package
library("tidyverse") # Data manipulation
library("progress") # Progress of the cross-validation
library("leaflet") # Geo maps
library("sp") # Spatial data type

# 0. Show how to crossvalidate
t0 <- lubridate::ymd_hms("2017-02-03 20:00:00")
t0seq <- t0 + lubridate:::weeks(1:5)

# 1. Read Data ===================================================
# rides <- read_csv("D:/Datasets/Rides_A.csv")
rides <- read_csv("C:/Users/mbg877/Google Drive/P1_Ride_Austin/00_Data/Clean_Database/Rides_A.csv")

austin5w <- data.frame()
for (i in 1:5) {
  period <- rides %>%
    filter(started_on >= t0seq[i] & started_on < (t0seq[i] + lubridate::hours(1))) %>% 
    # select(started_on, start_location_lat, start_location_long) %>% 
    mutate(week = i)
  austin5w <- rbind(austin5w, period)
}



