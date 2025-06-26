library(tidyverse)
library(jsonlite) # for loading parameters
library(GA) # for Genetic Algorithm
source("R/simulation.R")
source("R/load_data.R")
source("R/plot_fit.R")

# Load the experimental data we are going to use for the fitting ----
data <- load_median_raw_data() %>%
  rename(logV_obs = logV, E_obs = E, D_obs = D)

# Load working parameter set ----
parms <- fromJSON("parameter-values/json/params.json") %>%
  enframe() %>%
  unnest(value) %>%
  deframe()
labels <- names(parms)

# Define cost function ----
cost <- function(parms, data) {
  # GA passes in unnamed vectors but our ODE function needs parameters names
  names(parms) <- labels
  control <- simulate(parms, V01 = 0, V02 = 0, D0 = parms[["D0"]], E0 = parms[["E0"]]) %>%
    mutate(Type = "Control")
  treatment <- simulate(parms, V01 = parms[["V01"]], V02 = parms[["V02"]], D0 = parms[["D0"]], E0 = parms[["E0"]]) %>%
    mutate(Type = "Treatment")
  bind_rows(control, treatment) %>%
    inner_join(data, by = c("Day", "Infection", "Type")) %>%
    mutate(
      logV_se = (logV - logV_obs)**2,
      E_se = (E - E_obs)**2,
      D_se = (D - D_obs)**2
    ) %>%
    summarise(sse = sum(logV_se + E_se + D_se, na.rm = TRUE)) %>%
    pull(sse)
}

# Determine upper and lower bounds ----
# We know that the threshold is not reached on the first infection so we can use
# this maximum D observation as a lower bound
width_factor <- 1.2
lower <- (1 / width_factor) * parms
upper <- width_factor * parms
lower["DT"] <- data %>%
  filter(Type == "Treatment", Infection == 1) %>%
  select(D_obs) %>%
  max()

# Perform Genetic Algorithm ----
# Maximise fitness <=> minimise cost
GA <- ga(
  type = "real-valued",
  fitness = function(parms) -cost(parms, data),
  lower = lower,
  upper = upper,
  maxiter = 100,
  parallel = 8
)

# Save the resulting parameters ----
new_parms <- GA@solution[1, ]
# filename <- paste(round(-1 * GA@fitnessValue, 0), GA@maxiter, "parameter_values.csv", sep = "_")
# write_csv(enframe(new_parms), file.path(dir, filename))
