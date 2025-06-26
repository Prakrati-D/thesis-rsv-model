library(tidyverse)
library(jsonlite) # for loading parameters
source("R/simulation.R")
source("R/load_data.R")
theme_set(theme_bw())

plot_fold_change <- function(pred, data, DT) {
  data_hline <- tibble(name = c("D", "E", "logV"), value = c(DT, NA, NA))
  ggplot() +
    geom_hline(data = datca_hline, mapping = aes(yintercept = value)) +
    geom_line(data = pivot_longer(pred, -c(Day, Infection)), mapping = aes(Day, value, colour = factor(Infection))) +
    geom_point(data = pivot_longer(data, -c(Day, Infection)), mapping = aes(Day, value, colour = factor(Infection))) +
    labs(y = "Fold change", colour = "Infection") +
    facet_grid(rows = vars(name), scales = "free_y")
}

plot_raw <- function(pred, data, DT, title = NULL) {
  data_hline <- tibble(name = c("D", "E", "logV"), value = c(DT, NA, NA))
  ggplot() +
    geom_hline(data = data_hline, mapping = aes(yintercept = value)) +
    geom_line(data = pivot_longer(pred, -c(Day, Infection, Type)), mapping = aes(Day, value, colour = factor(Infection))) +
    geom_point(data = pivot_longer(data, -c(Day, Infection, Type)), mapping = aes(Day, value, colour = factor(Infection))) +
    labs(title = title, y = NULL, colour = "Infection") +
    facet_grid(rows = vars(name), cols = vars(Type), scales = "free_y")
}

data <- load_median_raw_data()
parms <- fromJSON("parameter-values/json/params.json") %>%
  enframe() %>%
  unnest(value) %>%
  deframe()

control <- simulate(parms, V01 = 0, V02 = 0, D0 = parms[["D0"]], E0 = parms[["E0"]]) %>%
  mutate(Type = "Control")
treatment <- simulate(parms, V01 = parms[["V01"]], V02 = parms[["V02"]], D0 = parms[["D0"]], E0 = parms[["E0"]]) %>%
  mutate(Type = "Treatment")
pred <- bind_rows(control, treatment) %>%
  select(-c(threshold_reached))

plot_raw(pred, data, parms[["DT"]])
ggsave("R/plot_fit.png", width = 9.25, height = 5.85)
