library(tidyverse)
source("R/load_data.R")
theme_set(theme_bw())

raw_data <- load_raw_data()
fold_change_data <- load_fold_change_data()
data <- bind_rows(raw_data, fold_change_data)

# Plot
data %>%
  filter(Day == 15) %>% # Note: day 15 is also the start of infection 2
  mutate(Infection = 2) %>%
  bind_rows(., data) %>%
  select(-E) %>% # New model doesn't contain E
  pivot_longer(-c(Infection, Day, Type)) %>%
  ggplot(aes(x = Day, y = value, colour = factor(Infection))) +
  geom_point() +
  geom_smooth(method = "loess", formula = "y ~ x") +
  labs(y = NULL, colour = "Infection") +
  facet_grid(rows = vars(name), cols = vars(Type), scales = "free_y")
ggsave("R/data_visualisation.png", width = 9.25, height = 5.85)
