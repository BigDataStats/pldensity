# 0. Libraries ===================================
library(tidyverse)
library(leaflet)
library(sp)
library(pldensity)

# 1. Read Data =====================================
rides <- read_csv("D:/Datasets/Rides_A.csv")
rides <- rides %>%
  mutate(datehour = lubridate::ymd_h(
    paste(paste(
      lubridate::year(started_on), 
      lubridate::month(started_on), 
      lubridate::day(started_on), sep = "-"), 
      lubridate::hour(started_on))
    )
  ) %>% 
  filter(
    started_on >= lubridate::ymd("2016-07-01"),
    started_on <=lubridate::ymd("2016-12-31")
  )
rides_count <- rides %>%
  group_by(datehour) %>%
  count()

plot(tail(rides_count$datehour, 500), tail(rides_count$n, 500), type = "l",
     main = "Ride count every hour in December '16")
abline(h = mean(rides_count$n), col = "red")
