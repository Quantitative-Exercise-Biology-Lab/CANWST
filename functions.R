library(zoo)

# read gps data in
read_gps_data <- function(file_name){
  gps_data <- 
    # skip to the head row
    read_csv(file_name, skip = 3, guess_max = 1000000) %>%
    # rename columns
    rename_all(to_snake_case) %>%
    # select the columns we need
    select(time_stamp, odometer, smoothed_velocity, gnss_hdop, acceleration_side, acceleration_forward)
  
  return(gps_data)
}

split_data_into_sessions <- function(data, time, start, end){
  session_data <-
    data %>%
    # add the offset of time that the device was started at to the timestamp
    mutate(clock_time = as_hms(time_stamp) + as_hms(time),
           # take the remainder for when the device was turned on for a certain day but the session overlaps with the next day
           clock_time = if_else(clock_time > (24 * 60 * 60), clock_time - (24 * 60 * 60), clock_time)) %>%
    # based on the game start time we manually figured out 
    # filter out any pregame movement
    filter(clock_time > as_hms(start), clock_time < as_hms(end)) %>%
    mutate(odometer_is_not_na = if_else(is.na(odometer), 0, 1),
           odometer_na_group  = cumsum(odometer_is_not_na)) %>%
    # remove the first group since this is before the gps is turned on
    filter(odometer_na_group != 0) %>%
    group_by(odometer_na_group) %>%
    mutate(length_of_odometer_interval = n()) %>%
    drop_na(odometer) %>%
    mutate(gps_error_group = if_else(length_of_odometer_interval > 20 | gnss_hdop > 2, 1, 0),
           gps_error_group = cumsum(gps_error_group)) %>%
    filter(length_of_odometer_interval <= 20 & gnss_hdop <= 2) %>%
    uncount(length_of_odometer_interval) %>%
    nest(data = c(-gps_error_group))
  
  return(session_data)
}

# read date and time from gps file
read_date_time <- function(file_name){
  date_time <- 
    read_csv(file_name, n_max = 1, col_names = F, col_types = cols(
      X1 = col_character()
    )) %>%
    mutate(date = parse_date_time(X1, orders = "mdy T"),
           time = paste(hour(date), minute(date), second(date), sep = ":"),
           date = date(date)) %>%
    select(date, time)
  
  return(date_time)
}


mean_no_na <- function(x){
  mean(x, na.rm = T)
}

sd_no_na <- function(x){
  sd(x, na.rm = T)
}

# create accelaration summary
# and grab the start and stop frame for the max mean vel interval
pull_interval_data <- function(data, r, dur){
  filter(data, row_number() >= r, row_number() < r + dur) %>%
    summarise(mean_sum_abs_acc = mean(abs(acceleration_side) + abs(acceleration_forward)),
              r_start = r,
              r_stop = r + dur)
}

max_mean_vel <- function(gps_data, dur){
  
  mmv <-
    gps_data %>%
    mutate(max_mean_vel = rollmean(smoothed_velocity, k = dur, fill = NA, align = "left"),
           r = row_number()) %>%
    filter(max_mean_vel == max(max_mean_vel, na.rm = T)) %>%
    select(max_mean_vel, r) %>%
    distinct(max_mean_vel, .keep_all = T) %>%
    mutate(interval_data = map(r, ~ pull_interval_data(gps_data, ., dur))) %>%
    unnest_wider(interval_data) 
  
  return(mmv)
}

fit_extended_model <- function(data, w_prime_guess, vel_crit_guess, max_vel_guess){
  nlsLM(max_mean_vel ~ 
          vel_crit + 
          w_prime * (1 - exp(-1 * duration * (max_vel - vel_crit)/w_prime))/duration - A * log(duration/18000) * (duration <= 18000)
        ,
        data = data,
        start = list(w_prime  = w_prime_guess, 
                     vel_crit = vel_crit_guess, 
                     max_vel  = max_vel_guess ,
                     A = 0.2
        ),
        lower = c(0, 0, 0, 0))
}

# extended model without the extended duration
fit_extended_model <- function(data, w_prime_guess, vel_crit_guess, max_vel_guess){
  nlsLM(max_mean_vel ~ 
          vel_crit + 
          w_prime * (1 - exp(-1 * duration * (max_vel - vel_crit)/w_prime))/duration,
        data = data,
        start = list(w_prime  = w_prime_guess, 
                     vel_crit = vel_crit_guess, 
                     max_vel  = max_vel_guess ),
        lower = c(0, 0, 0))
}

fit_five_p <- function(data, vel_crit_guess, max_vel_guess){
  nlsLM(max_mean_vel ~
          vel_crit + 
          (max_vel - vel_crit) / ((1 + exp(-a *(log(duration) - b) ) )^f),
        data = data,
        start = list(a = -2.71, 
                     b = 2.195, 
                     vel_crit = vel_crit_guess,
                     max_vel  = max_vel_guess, 
                     f = 0.2),
        control = nls.lm.control(maxiter = 300))
  
}

fit_three_p <- function(data, w_prime_guess, vel_crit_guess, max_vel_guess){
  nlsLM(max_mean_vel ~ 
          vel_crit + 
          w_prime / (duration + w_prime/(max_vel - vel_crit)), 
        data = data, 
        start = list(w_prime  = w_prime_guess, 
                     vel_crit = vel_crit_guess, 
                     max_vel  = max_vel_guess),
        lower = c(0, 0, 0),
        control = nls.lm.control(maxiter = 300))
}  



fit_two_p <- function(data, w_prime_guess, vel_crit_guess){
  nlsLM(max_mean_vel ~
          w_prime / duration +
          vel_crit,
        data = data,
        start = list(w_prime = w_prime_guess,
                     vel_crit = vel_crit_guess),
        control = nls.lm.control(maxiter = 300))
}
