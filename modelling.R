library(patchwork)
library(tidyverse)
library(broom)
library(minpack.lm)

source("functions.R")

theme_set(theme_bw() +
            theme(strip.background = element_rect(fill = "white")))


temperature <- 
  read_csv("temperature_measurements.csv")


max_mean_vel_data <-
  read_rds("max_mean_vel_data_acc.rds") %>%
  select(-r, -mean_sum_abs_acc, - r_start, - r_stop)
  




max_mean_vel_data <-
  max_mean_vel_data %>%
  mutate(duration = duration/10)


max_mean_vel_data %>%
  group_by(device_id, duration) %>%
  filter(max_mean_vel == max(max_mean_vel)) %>%
  ungroup() %>%
  count(device_id, type, session, date, description)

max_mean_vel_data %>%
  group_by(device_id, duration) %>%
  filter(max_mean_vel == max(max_mean_vel)) %>%
  ggplot(aes(x = duration, y = max_mean_vel, colour = description)) +
  geom_text(aes(label = description)) +
  facet_wrap(~ device_id)


global_max_mean_vel <- 
  max_mean_vel_data %>%
  group_by(device_id, duration) %>%
  filter(max_mean_vel == max(max_mean_vel)) %>%
  select(device_id, duration, max_mean_vel) %>%
  distinct() %>%
  ungroup() %>%
  nest(data = c(-device_id)) %>%
  mutate(description = "global")

date_max_mean_vel <-
  max_mean_vel_data %>%
  group_by(device_id, duration, date) %>%
  filter(max_mean_vel == max(max_mean_vel)) %>%
  select(date, device_id, duration, max_mean_vel) %>%
  distinct() %>%
  ungroup() %>%
  nest(data = c(-device_id, -date)) %>%
  mutate(description = "date_max")

data_to_model <-
  max_mean_vel_data %>%
  ungroup() %>%
  nest(data = c(duration, max_mean_vel)) %>% 
  bind_rows(global_max_mean_vel) %>%
  bind_rows(date_max_mean_vel) %>%
  mutate(cv_length = list(c(800))) %>%
  unnest(cols = c(cv_length)) %>%
  mutate(data = map2(data, cv_length, ~ filter(.x, duration <= .y)))


data_to_model <-
  data_to_model %>% 
  mutate(data = map(data, ~ mutate(., 
                                   inv_time = 1 / duration,
                                   w        = duration * (max_mean_vel - min(max_mean_vel)))),
         w_prime_guess  = map_dbl(data, ~ filter(., duration > 60, 
                                                 duration < 120) %>%
                                    summarise(w = mean(w)) %>%
                                    distinct() %>%
                                    pull(w)),
         vel_crit_guess = map_dbl(data, ~ summarise(., 
                                                    min_vel = min(max_mean_vel)) %>% 
                                    distinct() %>%
                                    pull(min_vel)),
         max_vel_guess  = map_dbl(data, ~ summarise(., 
                                                    max_vel = max(max_mean_vel)) %>% 
                                    distinct() %>%
                                    pull(max_vel)))

data_modelled <-
  data_to_model %>%
  mutate(two_p          = pmap(list(data, w_prime_guess, vel_crit_guess), 
                               possibly(fit_two_p, NA)),
         three_p        = pmap(list(data, w_prime_guess, vel_crit_guess, max_vel_guess),
                               possibly(fit_three_p, NA)),
         five_p         = pmap(list(data, vel_crit_guess, max_vel_guess), 
                               possibly(fit_five_p, NA)),
         extended_model = pmap(list(data, w_prime_guess, vel_crit_guess, max_vel_guess),
                               possibly(fit_extended_model, NA)))

data_modelled %>%
  pivot_longer(cols = c(two_p, three_p, five_p, extended_model), names_to = "model_name", values_to = "model_fit") %>%
  mutate(type = map_chr(model_fit, typeof)) %>%
  filter(type != "list") %>% View

data_modelled <-
  data_modelled %>%
  pivot_longer(cols = c(two_p, three_p, five_p, extended_model), names_to = "model_name", values_to = "model_fit") %>%
  mutate(type = map_chr(model_fit, typeof)) %>%
  filter(type == "list") %>% 
  mutate(fitted   = map(model_fit, fitted),
         tidy_mod = map(model_fit, tidy),
         sigma    = map_dbl(model_fit, ~ glance(.) %>% pull(sigma)))  


#write_rds(data_modelled, "data_modelled.rds")

data_modelled %>%
  select(-tidy_mod, -model_fit) %>%
  unnest(cols = c(data, fitted)) %>%
  select(device_id, model_name, sigma, fitted, duration, max_mean_vel, session, description, date) %>%
  # filter(description %in% c("global", "11v11coached", "11v11")) %>%
  filter(description %in% c("global")) %>%
  #filter(model_name %in% c("extended_model")) %>%
  ggplot() +
  geom_point(aes(x = duration, y = max_mean_vel)) +
  geom_line(aes(x = duration,  y = fitted, color = model_name)) +
  facet_wrap(~ device_id + session + description + date) +
  expand_limits(y = 0)

data_modelled %>%
  select(device_id, model_name, type, session, description, tidy_mod) %>%
  unnest(cols = tidy_mod) %>%
  filter(term == "vel_crit") %>%
  group_by(device_id, model_name, type, description) %>%
  summarise(mean_vel_crit = mean(estimate)) %>%
  arrange(-mean_vel_crit) 


data_modelled %>%
  select(model_name, type, session, description, tidy_mod) %>%
  unnest(cols = tidy_mod) %>%
  filter(model_name == "five_p", description == "global")


data_modelled %>%
  select(device_id, model_name, type, session, description, model_fit) %>%
  filter(description == "global") %>%
  mutate(g = map(model_fit, glance)) %>%
  unnest(cols = g)

data_modelled %>%
  select(device_id, model_name, type, session, description, model_fit) %>%
  filter(description == "global",
         device_id %in% c(47111, 51533)) %>%
  mutate(g = map(model_fit, glance)) %>%
  unnest(cols = g) %>%
  select(device_id, model_name, AIC, BIC) %>%
  arrange(device_id)



data_modelled %>%
  filter(description == "global",
         device_id %in% c(47111, 51533)) %>%
  select(device_id, model_name, type, session, description, model_fit) %>%
  mutate(crit_vel_30 = map_dbl(model_fit, ~ predict(., tibble(duration = 1800)))) %>%
  ungroup() %>%
  select(device_id, model_name, crit_vel_30) %>%
  arrange(device_id)

annotate_data <- 
  data_modelled %>% 
  filter(model_name == "extended_model",
         device_id %in% c(47111, 51533)) %>%
  filter(description == "global" | (date == "2018-06-05" & description == "11v11free")) %>%
  unnest(cols = c(data, fitted)) %>%
  mutate(duration = duration / 10,
         description = case_when(description == "global" ~ "All Sessions",
                                 description == "11v11free" ~ "One Session"),
         description = factor(description),
         description = fct_rev(description),
         crit_vel = map_dbl(model_fit, ~ predict(., tibble(duration = 1800))),
         crit_vel = round(crit_vel, 2)) %>%
  group_by(description, device_id, crit_vel) %>%
  summarise(rmse = sqrt( sum( (fitted - max_mean_vel)^2))) %>%
  ungroup() %>%
  mutate(rmse = round(rmse, 1),
         device_player = case_when(device_id == 47111 ~ "Player 1",
                               device_id == 51533 ~ "Player 2"),
         label = paste0(device_player, " RMSE = ", rmse)) %>%
  group_by(description) %>%
  arrange(crit_vel) %>%
  mutate(y_plot_cv = if_else(2 == row_number(), crit_vel + 0.3, crit_vel - 0.3),
         label_cv = paste0("AV = ", crit_vel))

data_modelled %>%
  filter(model_name == "extended_model",
         device_id %in% c(47111, 51533)) %>%
  filter(description == "global" | (date == "2018-06-05" & description == "11v11free")) %>%
  unnest(cols = c(data, fitted)) %>%
  mutate(duration = duration,
         description = case_when(description == "global" ~ "All Sessions",
                                 description == "11v11free" ~ "One Session"),
         description = factor(description),
         description = fct_rev(description)) %>%
  group_by(description, device_id) %>%
  mutate(rmse = sqrt( sum( (fitted - max_mean_vel)^2))) %>%
  ggplot() +
  geom_point(aes(x = duration,  y = fitted, color = factor(device_id)), size = 2) +
  facet_wrap(~ description) +
#  geom_text(data = annotate_data, aes(label = label), 
 #           x = Inf, 
  #          y = -Inf ,
   #         hjust   = 1.05,
    #        vjust   = c(- 1.5, -3, -3, -1.5)) +
  geom_hline(aes(yintercept = crit_vel, color = factor(device_id)), data = annotate_data) +
  geom_text(data = annotate_data, aes(label = label_cv, y = y_plot_cv, color = factor(device_id)), 
            x = -Inf,
            hjust   = -0.1, show.legend = F) +
  # font 
  # and moving y axis to 0
  expand_limits(y = 0) +
  scale_x_continuous(expand = expand_scale(mult = c(0, 0.05))) +
  scale_y_continuous(expand = expand_scale(mult = c(0, 0.05))) +
  labs(x = "Duration (sec)",
       y = "Average Velocity (m/sec)",
       title = ) +
  scale_color_manual(name = NULL, values = c("#D7030D", "black"), labels = c("Player 1", "Player 2")) +
  theme(text = element_text(size = 22)) + 
  theme(
    legend.position = c(0.99, .99),
    legend.justification = c("right", "top"),
    legend.box.just = "right",
    legend.margin = margin(0, 6, 6, 6)
  ) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)))

ggsave("plots/one_session_to_all_session.png")


#player_1 <- 
max_mean_vel_data %>%
  filter(duration <= 800) %>%
  filter(device_id %in% c(47111, 51533, 51516, 51528)) %>%
  group_by(device_id, duration) %>%
  filter(max_mean_vel == max(max_mean_vel)) %>% 
  ungroup() %>%
  mutate(device_id = case_when(device_id == 47111 ~ "Player 1",
                                device_id == 51533 ~ "Player 2",
                               device_id == 51516 ~ "Player 3",
                               device_id == 51528 ~ "Player 4")) %>%
  ggplot(aes(x = duration, y = max_mean_vel, fill =  description#shape = description
             ), colour = "black") +
  geom_point(size = 3, pch = 21) +
  facet_wrap(~ device_id, nrow = 2) +
  scale_fill_manual(name = ~ underline("Session Type"), 
                    values =
                       c("#D7030D",
                         "Black",
                         "pink",
                         "#0761fb",
                         "white",
                         "#d4a000",
                         "green"
                         ),
                    labels = c("11v11 Coached", 
                               "11v11 Free", 
                               "Costa Rica", 
                               "Germany", 
                               "Speed", 
                               "USA", 
                               "Warm up")
                     ) +
    scale_x_continuous(expand = expand_scale(mult = c(0, 0.05))) +
    labs(x = "Duration (sec)",
         y = "Average Velocity (m/sec)") +
    theme(text = element_text(size = 30),
          legend.title = element_text(size = 25)) + 
  theme(
      legend.position = c(.99, .99),
      legend.justification = c("right", "top"),
      legend.box.just = "right",
      legend.margin = margin(6, 6, 6, 6)
    ) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
        legend.title.align = 0.5) 
  
ggsave("plots/player_1_2_raw.png", height = 11.91, width = 13.58, units = "in")





new_data <- 
  data_modelled %>% 
  filter(device_id %in% c(47111, 51533, 50975, 50947)) %>%
  left_join(temperature, by = "date") %>%
  mutate(crit_vel = map_dbl(model_fit, ~ predict(., tibble(duration = 1800))),
         crit_vel = round(crit_vel, 2)) %>%
  filter(model_name == "extended_model",
         description == "date_max") %>%  
  select(temperature, crit_vel, device_id)

new_data %>%
  bind_rows(old_data) %>%
  filter(!(device_id %in% c(11148))) %>%
  filter(crit_vel > 1.5)  %>%
  mutate(device_id = case_when(device_id == 47111 ~ "Player 1",
                               device_id == 51533 ~ "Player 2",
                               # device_id == 51516 ~ "Player 3",
                               # device_id == 51528 ~ "Player 4",
                               device_id == 4650  ~ "Player 3",
                               device_id == 10658 ~ "Player 4",
                               device_id == 10664 ~ "Player 5",
                               device_id == 11240 ~ "Player 6",
                               device_id == 50947 ~ "Player 7",
                               device_id == 50975 ~ "Player 8"),
#         device_id = factor(device_id),
 #        device_id = fct_relevel(device_id, "Player 10", after = 10)
)%>%
  ggplot(aes(x = temperature, y = crit_vel)) +
  geom_point() +
  geom_smooth(method = "lm", se = F, colour = "#D7030D") +
  facet_wrap(~ device_id, nrow = 2) +
  labs(x = "Temperature (Celsius)",
       y = "Asymptotic Velocity (m/sec)") +
  theme(text = element_text(size = 34)) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)))

ggsave("plots/temperature.png", height = 11.91, width = 13.58, units = "in")

# start and end of interval when we rerun the data
# temperature plot
# intermittent efforts just the box for a session for 2 points