library(tidyverse)
library(lubridate)
library(rnoaa)
library(lutz)
library(hms)

source("functions.R")

# list all the files 
# add them to a data frame and create identifier columns
data_files <-
  tibble(full_path = list.files("data", recursive = T, full.names = T)) %>%
  filter(str_detect(full_path, ".csv")) %>%
  mutate(file_name = str_extract(full_path, "(\\d|_)+.csv")) %>%
  distinct(file_name, .keep_all = T) %>%
  separate(file_name, c("device_id", "centitime", "unknown"), sep = "_", convert = T) %>%
  drop_na(device_id) %>%
  mutate(date_time = map(full_path, read_date_time)) %>%
  unnest(date_time) 

session_annotations <-
  read_csv("data/annotations/master.csv") %>%
  group_by(date, type, session, description) %>%
  filter(n() > 2) %>%
  drop_na(start_time_utc, end_time_utc) %>%
  summarise(start = min(start_time_utc, na.rm = T), 
            end   = max(end_time_utc, na.rm = T))

data_files  <-
  data_files %>%
  inner_join(session_annotations, by = "date") %>%
  arrange(start) %>%
  distinct(date, .keep_all = T) %>%
  mutate(first_guaranteed_on = 100 * as.numeric(as_hms(start) - as_hms(time)),
         first_guaranteed_on = if_else(first_guaranteed_on < 0,
                                       24*60*100*60 + first_guaranteed_on, first_guaranteed_on))



long_lat_data <-
  data_files %>%
  mutate(data = map2(full_path,
                     first_guaranteed_on,
                     ~ read_csv(.x, skip = .y, guess_max = 1000, n_max = 5000, col_names = F) %>%
                       select(c(17, 18, 26)) %>%
                       rename(latitude = X17,
                              longitude = X18,
                              gnss_hdop = X26) %>%
                       # drop na latitute
                       drop_na(latitude) %>%
                       # HDOP
                       filter(gnss_hdop == min(gnss_hdop)) %>%
                       slice(1))) %>%
  unnest(cols = c(data)) 

long_lat_data 
  
isd_data <-
  long_lat_data %>%
  mutate(isd_station_id = map2(latitude, longitude, ~ isd_stations_search(.x, .y, 200) %>%
                                                      slice(1)))


new_isd <- 
  isd_data %>% 
  rename(start_time = start,
         end_time = end) %>%
  unnest(cols = c(isd_station_id)) %>%
  mutate(date_time = ymd_hms(paste0(date, time), tz = "UTC")) %>%
  mutate(time_zone = tz_lookup_coords(latitude, longitude),
         new_time = map2(date_time, time_zone, force_tz))  %>% 
  unnest(cols = c(new_time)) %>%
  mutate(new_time = round_date(new_time , unit = "hour"),
         local_hour = hour(new_time)) %>%
  mutate(isd = pmap(list(usaf, wban, year(date)), possibly(~ isd(..1, ..2, ..3) %>%
                                                             select(date, time, temperature, elevation), NA)))
  

temperature <- 
  new_isd %>%
  group_by(full_path) %>%
  mutate(type = map_chr(isd, typeof)) %>% 
  filter(type == "list") %>%
  unnest(cols = c(isd), names_sep = "_") %>% 
  filter(date == ymd(isd_date)) %>%
  mutate(isd_time = parse_number(isd_time)/100,
         time_dist = abs(local_hour - isd_time)) %>% 
  filter(time_dist == min(time_dist)) %>%
  mutate(temperature = parse_number(isd_temperature)/10) %>%
  ungroup() %>%
  select(date, temperature)

write_csv(temperature, "temperature_measurements.csv")

