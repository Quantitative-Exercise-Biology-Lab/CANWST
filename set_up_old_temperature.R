library(tidyverse)
library(lubridate)
library(broom)

data_3_players <- read_csv("old_methodology_data/3_players_max_vel_v2.csv")

data_ming_games <- read_csv("old_methodology_data/all_ming_games.csv")

data_18_19 <- read_csv("old_methodology_data/2018_2019_data.csv")

data_3_players <-
  data_3_players %>%
  unite(col = date, c("year", "month", "day"), sep = "-") %>%
  mutate(date = ymd(date))

data_ming_games <- 
  data_ming_games %>%
  select(date, device_id, length_of_interval, max_vel)

data_3_players <-
  data_3_players %>%
  select(date, device_id, length_of_interval, max_vel)


data_18_19 <-
  data_18_19 %>%
  unite(col = date, c("year", "month", "day"), sep = "-") %>%
  mutate(date = ymd(date))  %>%
  select(date, device_id, length_of_interval, max_vel)


data <- 
  bind_rows(data_ming_games, data_3_players) %>%
  bind_rows(data_18_19)



# Remove if they played less than 30 minutes
# and if they have any maxes as less than 0.1
# AND REMOVE DISCONTINUITIES
data <-
  data %>%
  group_by(date, device_id) %>%
  mutate(max_length = max(length_of_interval),
         num_zeros = sum(max_vel < 0.1)) %>%
  ungroup() %>%
  filter(max_length == 1800) %>%
  filter(num_zeros == 0) %>%
  filter(length_of_interval <= 240) %>%
  group_by(date, device_id) %>%
  mutate(discontinuity = if_else(max_vel - lead(max_vel) > 1 & length_of_interval > 5, 1, 0),
         group_possible_error = max(discontinuity, na.rm = T)) %>%
  filter(group_possible_error == 0) %>%
  ungroup() %>%
  select(- max_length, - num_zeros, - discontinuity, - group_possible_error)


data_date <- data 

data_to_model <-
  data_date %>%
  rename(max_vel_interval = max_vel) %>%
  nest(data = c(-device_id, -date)) %>%
  mutate(max_int = 240) %>%
  mutate(data = map2(data, max_int, ~ group_by(.x, length_of_interval) %>%
                       mutate(max_vel_rank = rank(-max_vel_interval, ties.method = "random")) %>%
                       filter(max_vel_rank == 1) %>% 
                       select(-max_vel_rank) %>%
                       ungroup() %>%
                       filter(length_of_interval <= .y)))


data_to_model <-
  data_to_model %>%
  mutate(data = map(data, ~ mutate(., 
                                   distance = max_vel_interval * length_of_interval,
                                   inv_time = 1 / length_of_interval,
                                   w        = length_of_interval * (max_vel_interval - min(max_vel_interval)))),
         w_prime_guess  = map_dbl(data, ~ filter(., length_of_interval > 60, 
                                                 length_of_interval < 120) %>%
                                    summarise(w = mean(w)) %>%
                                    pull(w)),
         vel_crit_guess = map_dbl(data, ~ summarise(., 
                                                    max_vel = min(max_vel_interval)) %>% 
                                    pull(max_vel)),
         max_vel_guess  = map_dbl(data, ~ summarise(., 
                                                    max_vel = max(max_vel_interval)) %>% 
                                    pull(max_vel)))

# data_modelled <-
#   data_to_model %>%
#   filter(w_prime_guess != 0) %>%
#   mutate(extended_model = pmap(list(data, w_prime_guess, vel_crit_guess, max_vel_guess), 
#                                possibly(
#                                ~ nlsLM(formula = max_vel_interval ~ 
#                                                  vel_crit + 
#                                                  w_prime * (1 - exp(-1 * length_of_interval * (max_vel - vel_crit)/w_prime))/length_of_interval -
#                                                  A * log(length_of_interval/1800) * (length_of_interval <= 1800),
#                                        data    = ..1,
#                                        start   = list(w_prime  = ..2, 
#                                                       vel_crit = ..3, 
#                                                       max_vel  = ..4, 
#                                                       A = 0.2),
#                                        lower   = c(0, 0, 0, -Inf),
#                                        upper   = c(Inf, 11, 11, Inf)), NA)))


data_modelled <-
  data_to_model %>%
  filter(w_prime_guess != 0) %>%
  mutate(extended_model = pmap(list(data, w_prime_guess, vel_crit_guess, max_vel_guess), 
                               possibly(
                                 ~ nlsLM(formula = max_vel_interval ~ 
                                           vel_crit + 
                                           w_prime * (1 - exp(-1 * length_of_interval * (max_vel - vel_crit)/w_prime))/length_of_interval,
                                         data    = ..1,
                                         start   = list(w_prime  = ..2, 
                                                        vel_crit = ..3, 
                                                        max_vel  = ..4),
                                         lower   = c(0, 0, 0),
                                         upper   = c(Inf, 11, 11)), NA)))

data_modelled <-
  data_modelled %>%  
  mutate(type = map_chr(extended_model, typeof),
         model_fit = extended_model) %>% 
  filter(type == "list") %>%
  mutate(tidy_mod = map(model_fit, tidy),
         fitted   = map(model_fit, fitted)) #%>%
#select(device_id, data, model_name, model_fit, fitted, tidy_mod, sigma) 



data_modelled %>% 
  select(date, device_id, data, fitted, tidy_mod) %>% 
  unnest(cols = c(tidy_mod)) %>%
  filter(term %in% c("vel_crit", "w_prime")) %>%
  select(-std.error, -statistic, -p.value) %>%
  spread(term, estimate) %>%
  arrange(vel_crit) %>%
  filter(vel_crit < 1.5) %>%
  mutate(vel_crit = round(vel_crit, 2),
         w_prime = round(w_prime, 2)) %>%
  unnest(cols = c(data, fitted)) %>%
  ggplot() +
  geom_point(aes(x = length_of_interval, y = max_vel_interval), size = 0.5) +
  geom_line(aes(x = length_of_interval, y = fitted)) +
  geom_text(aes(label = paste0("CV = ", vel_crit), x = Inf, y = Inf), 
            hjust = 1.25,
            vjust = 1.25) +
  geom_text(aes(label = paste0("W' = ", w_prime), x = Inf, y = Inf), 
            hjust = 1.25,
            vjust = 4.25) +
  facet_wrap(~ date + device_id) 



data_modelled <-
  data_modelled %>%
  mutate(cv_estimate = map_dbl(tidy_mod, . %>% filter(term == "vel_crit") %>% pull(estimate))) %>%
  filter(cv_estimate >= 1.5)



game_results <- 
  read_csv("old_methodology_data/wins losses and temprature.csv")

practice <-
  read_csv("old_methodology_data/temperature_measurements.csv") %>%
  filter(temperature < 50)


data_modelled <-
  data_modelled %>%
  filter(device_id %in% c(11240, 11148, 10664, 10658, 4650))

old_data <-
  data_modelled %>%
  select(date, device_id, max_int, tidy_mod) %>%
  unnest(cols = c(tidy_mod)) %>%
  filter(term %in% c("vel_crit")) %>%
  inner_join(game_results, by = "date") %>%
  group_by(device_id, term) %>%
  filter(n() > 4) %>%
  ungroup() %>%
  select(temperature = `Day Average Temp`, crit_vel = estimate, device_id)

#vc <- 
data_modelled %>%
  select(date, device_id, max_int, tidy_mod) %>%
  unnest() %>%
  filter(term %in% c("vel_crit")) %>%
  inner_join(game_results, by = "date") %>%
  group_by(device_id, term) %>%
  filter(n() > 4) %>%
  mutate(`Player ID` = device_id) %>%
  ggplot(aes(x = `Day Average Temp`, y = estimate)) +
  geom_point() +
  facet_grid(~ `Player ID`, scales = "free_y", labeller = label_both) +
  geom_smooth(method = "lm", se = F, colour = "red") +
  theme_bw(base_size = 22) +
  labs(title =  ~ underline("Asymptotic Velocity (AV) Estimate in Games vs. Temperature by Player"),
       y = "AV Estimate (m/sec)",
       x = "Temperature (Celsius)") + 
  theme(plot.title = element_text(hjust = 0.5),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
        axis.text = element_text(size = 20),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0))) #+
ylim(c(0,NA))