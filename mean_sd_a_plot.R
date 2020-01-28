library(ggrepel)

max_mean_vel_data <-
  read_rds("max_mean_vel_data_acc_51516.rds")

max_mean_vel_data <-
  max_mean_vel_data %>%
  mutate(duration = duration / 10)






max_mean_vel_data %>%
  group_by(device_id, duration) %>% 
  filter(mean_vel == max(mean_vel)) %>%
  mutate(sd = map_dbl(interval, ~ sd(.$smoothed_velocity)),
         a = map_dbl(interval, ~ mean(abs(.$acceleration_side) + abs(.$acceleration_forward))),
         a = if_else(duration == 0.1, NA_real_, a)) %>%
  arrange(duration) %>%
  select(-interval) %>%
  distinct() %>% 
  ungroup() %>%
  mutate(v_a_label = if_else(row_number() %in% c(2, 4, 8, 12, 14), paste("V/A =", round(mean_vel/a, 2)), "")) %>% 
  ggplot() +
  geom_point(aes(x = duration, y = mean_vel, colour = "Average"), size = 3) +
  geom_point(aes(x = duration, y = sd, colour = "SD"), size = 3) +
  geom_line(aes(x = duration, y = mean_vel, colour = "Average")) +
  geom_line(aes(x = duration, y = sd, colour = "SD")) +
  #geom_point(aes(x = duration, y = a, colour = "Acceleration")) +
  #geom_text_repel(aes(x = duration, y = sd, label = v_a_label), nudge_y = 0.5, size = 8) +
  expand_limits(y = 0) +
  scale_x_continuous(expand = expand_scale(mult = c(0, 0.05))) +
  scale_y_continuous(expand = expand_scale(mult = c(0, 0.05))) +
  labs(x = "Duration (sec)",
       y = "Velocity (m/sec)") +
  theme(text = element_text(size = 38)) +
  scale_color_manual(NULL, values = c("#D7030D", "#0761fb")) + 
  theme(
    legend.position = c(.99, .99),
    legend.justification = c("right", "top"),
    legend.box.just = "right",
    legend.margin = margin(0, 6, 6, 6)
  ) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)))

ggsave("plots/intermittent_effects.png", height = 7.02, width = 12.92, units = "in")




max_mean_vel_data %>%
  group_by(device_id, duration) %>% 
  filter(mean_vel == max(mean_vel)) %>%
  mutate(sd = map_dbl(interval, ~ sd(.$smoothed_velocity)),
         a = map_dbl(interval, ~ mean(abs(.$acceleration_side) + abs(.$acceleration_forward))),
         a = if_else(duration == 0.1, NA_real_, a)) %>%
  arrange(duration) %>%
  filter(duration == 280.1) %>%
  unnest(cols = c(interval)) %>%
  mutate(second = row_number()) %>%
  ggplot() +
  geom_line(aes(x = second/10, y = smoothed_velocity), size = 1) +
  expand_limits(y = 0) +
  scale_x_continuous(expand = expand_scale(mult = c(0, 0.05))) +
  scale_y_continuous(expand = expand_scale(mult = c(0, 0.05))) +
  labs(x = "Duration (sec)",
       y = "Velocity (m/sec)") +
  theme(text = element_text(size = 34)) +
  theme(strip.text = element_blank()) +
  theme(
    legend.position = c(.95, .95),
    legend.justification = c("right", "top"),
    legend.box.just = "right",
    legend.margin = margin(6, 6, 6, 6)
  ) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)))



max_mean_vel_data %>%
  group_by(device_id, duration) %>% 
  filter(mean_vel == max(mean_vel)) %>%
  mutate(sd = map_dbl(interval, ~ sd(.$smoothed_velocity)),
         a = map_dbl(interval, ~ mean(abs(.$acceleration_side) + abs(.$acceleration_forward))),
         a = if_else(duration == 0.1, NA_real_, a)) %>%
  arrange(duration) %>%
  filter(duration == 80.1) %>%
  unnest(cols = c(interval)) %>%
  mutate(second = row_number()) %>%
  ggplot() +
  geom_line(aes(x = second/10, y = smoothed_velocity), size = 1) +
  expand_limits(y = 0) +
  scale_x_continuous(expand = expand_scale(mult = c(0, 0.05))) +
  scale_y_continuous(expand = expand_scale(mult = c(0, 0.05))) +
  labs(x = "Duration (sec)",
       y = "Velocity (m/sec)") +
  theme(text = element_text(size = 34)) +
  theme(strip.text = element_blank()) +
  theme(
    legend.position = c(.95, .95),
    legend.justification = c("right", "top"),
    legend.box.just = "right",
    legend.margin = margin(6, 6, 6, 6)
  ) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)))


max_mean_vel_data %>%
  group_by(device_id, duration) %>% 
  filter(mean_vel == max(mean_vel)) %>%
  mutate(sd = map_dbl(interval, ~ sd(.$smoothed_velocity)),
         a = map_dbl(interval, ~ mean(abs(.$acceleration_side) + abs(.$acceleration_forward))),
         a = if_else(duration == 0.1, NA_real_, a)) %>%
  arrange(duration) %>%
  filter(duration == 80.1 | duration == 280.1) %>%
  unnest(cols = c(interval)) %>%
  mutate(second = row_number()) %>%
  ggplot() +
  geom_line(aes(x = second/10, y = smoothed_velocity, group = duration), size = 1) +
  expand_limits(y = 0) +
  facet_wrap(~ duration, scales = "free_x") +
  scale_x_continuous(expand = expand_scale(mult = c(0, 0.05))) +
  scale_y_continuous(expand = expand_scale(mult = c(0, 0.05))) +
  labs(x = "Duration (sec)",
       y = "Velocity (m/sec)") +
  theme(text = element_text(size = 34)) +
  theme(strip.text = element_blank()) +
  theme(
    legend.position = c(.95, .95),
    legend.justification = c("right", "top"),
    legend.box.just = "right",
    legend.margin = margin(6, 6, 6, 6)
  ) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)))
