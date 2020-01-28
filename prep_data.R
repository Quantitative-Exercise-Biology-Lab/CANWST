#################################################################################
# Libraries
#################################################################################
library(tidyverse)
library(snakecase)
library(hms)
library(tictoc)
library(lubridate)
library(broom)
library(minpack.lm)

#################################################################################
# User Defined Functions
#################################################################################

source("functions.R")

# duration limit should be minutes
# moving avg window should be in seconds
get_mean_max_vel_data <- function(device_ids, date_to_use, duration_limit, moving_avg_window){
  # Start timing 
  tic()
  # convert duration from minutes to 10th of a second
  duration_limit <- duration_limit * 60 * 10 + 1
  # convert moving avg window to 10th of a second
  moving_avg_window <- moving_avg_window * 10
  
  
  # read in old saved data
  old_max_mean_vel_data <- 
    read_rds("max_mean_vel_data_acc.rds")
  
  # select only columns we need
  obs_we_already_have <-
    old_max_mean_vel_data %>%
    select(full_path, type, session, description, duration)
  
  # list all the files 
  # add them to a data frame and create identifier columns
  files <-
    tibble(full_path = list.files("data", recursive = T, full.names = T)) %>%
    filter(str_detect(full_path, ".csv")) %>%
    mutate(file_name = str_extract(full_path, "(\\d|_)+.csv")) %>%
    distinct(file_name, .keep_all = T) %>%
    separate(file_name, c("device_id", "centitime", "unknown"), sep = "_", convert = T) %>%
    drop_na(device_id)
  
  
  # read in session annotations
  # because we can't currently link on device id and player
  # take the conservative approach of taking the min and max of a session
  session_annotations <-
    read_csv("data/annotations/master.csv") %>%
    group_by(date, type, session, description) %>%
    filter(n() > 2) %>%
    drop_na(start_time_utc, end_time_utc) %>%
    summarise(start = min(start_time_utc, na.rm = T), 
              end   = max(end_time_utc, na.rm = T))
  
  
  files_w_dt <- 
    files %>%
    filter(device_id %in% device_ids) %>%
    mutate(date_time = map(full_path, read_date_time)) %>%
    unnest(date_time) %>%
    filter(date == date_to_use)
  
  
  all_sessions <- 
    files_w_dt %>%
    filter(date %in% session_annotations$date) %>%
    mutate(data = map(full_path, read_gps_data)) %>%
    left_join(session_annotations, by = "date") %>% 
    drop_na(start) %>%
    # possibly suppresses errors so that it doesn't waste a run
    # but stores an NA in the column if there was an error
    mutate(data = pmap(list(data, time, start, end), possibly(split_data_into_sessions, NA)))
  
  # watch out when you have more than one player
  # all_sessions %>%
  #   anti_join(obs_we_already_have, by = c("full_path", "type", "session", "description")) %>%
  #   unnest(cols = c(data)) %>%
  #   mutate(group = paste(date, session, description, device_id, sep = "_")) %>%
  #   select(data, group) %>%
  #   mutate(data = map(data, ~ filter(., row_number() %% 2000 == 0))) %>%
  #   unnest(cols = c(data)) %>%
  #   ggplot(aes(x = time_stamp, y = smoothed_velocity)) +
  #   geom_point() +
  #   facet_wrap(~ group, scales = "free_x")
  # 
  max_mean_vel_data <-
    all_sessions %>%
    unnest(cols = c(data)) %>%
    mutate(duration   = list(seq(1, duration_limit, by = moving_avg_window)),
           length_of_activity = map_dbl(data, nrow)) %>%
    unnest(duration) %>%
    filter(duration < length_of_activity) %>%
    anti_join(obs_we_already_have, by = c("full_path", "type", "session", "description", "duration")) %>%
    mutate(max_mean_vel = map2(data, duration, max_mean_vel))
  
  
  
  new_max_mean_vel_data <-
    max_mean_vel_data %>%
    unnest(cols = c(max_mean_vel)) %>%
    select(-data)
  
  
  
  new_max_mean_vel_data <-
    new_max_mean_vel_data %>%
    bind_rows(old_max_mean_vel_data) 
  
  new_max_mean_vel_data %>%
    write_rds("max_mean_vel_data_acc.rds")
  
  toc()
  
}

device_ids <- 
  tibble(full_path = list.files("data", recursive = T, full.names = T)) %>%
  filter(str_detect(full_path, ".csv")) %>%
  mutate(file_name = str_extract(full_path, "(\\d|_)+.csv")) %>%
  distinct(file_name, .keep_all = T) %>%
  separate(file_name, c("device_id", "centitime", "unknown"), sep = "_", convert = T) %>%
  drop_na(device_id) %>%
  count(device_id, sort = T) %>%
  pull(device_id)

dates_to_run <- 
  read_csv("data/annotations/master.csv") %>%
  select(date) %>%
  distinct() %>%
  pull(date)




tic()
map(dates_to_run, possibly(~ get_mean_max_vel_data(51533, ., 15, 20), NA))
toc()

tic()
map(dates_to_run, possibly(~ get_mean_max_vel_data(47111, ., 15, 20), NA))
toc()

