library(tidyverse)
library(deSolve)

model <- function(t, state, parms) {
  with(as.list(c(state, parms)), {
    dVdt <- kV * (1 + aDV * D) * V - dV * V
    dDdt <- dBV * V * (Btot - D) + dBE * E * (Btot - D) - kB * D
    dEdt <- ifelse(threshold_reached != FALSE, kE + kED * D - dE * E, kE - dE * E)
    nudge <- ifelse(!threshold_reached & (D >= DT), 1, 0)
    return(list(c(dVdt, dDdt, dEdt, nudge)))
  })
}

simulate <- function(parms, V01, V02, D0, E0) {
  # First infection on day 5
  times1 <- seq(from = 5, to = 21, by = 1 / 24)
  yinit1 <- c(V = V01, D = D0, E = E0, threshold_reached = FALSE)
  infec1 <- ode(yinit1, times1, model, parms, method = "ode45") %>%
    as.data.frame() %>%
    tibble() %>%
    mutate(Infection = 1)
  # Second infection on day 15
  times2 <- seq(from = 15, to = 28, by = 1 / 24)
  yinit2 <- infec1 %>%
    filter(time == 15) %>%
    select(-c(time, Infection)) %>%
    as.numeric()
  names(yinit2) <- names(yinit1)
  yinit2["V"] <- yinit2["V"] + V02
  infec2 <- ode(yinit2, times2, model, parms, method = "ode45") %>%
    as.data.frame() %>%
    tibble() %>%
    mutate(Infection = 2)
  bind_rows(infec1, infec2) %>%
    mutate(logV = log(V), V = NULL) %>%
    rename(Day = time)
}
