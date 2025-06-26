library(tidyverse)
library(jsonlite) # for loading parameters
library(lhs) # Latin Hypercube Sampling
library(sensitivity) # PRCC function
library(boot) # for bootstrap replicates
library(latex2exp) # for latex parameters in plot labels
library(cowplot) # for subplots
source("R/simulation.R")
theme_set(theme_bw(base_size = 12))

set.seed(123) # Reproducibility for Latin Hypercube Sample

PRCC <- function(parms, infection_number, width_factor = 0.05, n = 500) {
  # Create a dataframe X with n rows and (no. of parameters) columns.
  # Each row is a random sample of parameter values between lower and
  # upper bounds using Latin Hypercube Sampling.
  lower <- (1 - width_factor) * parms
  upper <- (1 + width_factor) * parms
  m <- length(parms)
  lmat <- matrix(rep(lower, n), nrow = n, ncol = m, byrow = TRUE)
  umat <- matrix(rep(upper, n), nrow = n, ncol = m, byrow = TRUE)
  X <- data.frame(lmat + (umat - lmat) * maximinLHS(n, m))
  colnames(X) <- names(parms)

  # Loop through each row of X and simulate from that set of random parameters.
  # Set the maximum barrier damage D(t) value for the specific `infection_number`
  # in the vector y (containing the responses corresponding to the design of
  # experiments).
  y <- numeric(n)
  for (i in 1:n) {
    print(i) # just to keep track of progress to see if it gets stuck
    y[i] <- simulate(parms = X[i, ], V01 = X[i, "V01"], V02 = X[i, "V02"], D0 = X[i, "D0"], E0 = X[i, "E0"]) %>%
      filter(Infection == infection_number) %>%
      select(D) %>%
      max()
  }

  # Perform the PRCC analysis using the `sensitivity` package
  # see: https://search.r-project.org/CRAN/refmans/sensitivity/html/pcc.html
  alpha <- 0.05
  out <- pcc(
    X = X,
    y = y,
    rank = TRUE, # Partial RANK Correlation Coefficients (PRCC)
    nboot = 100, # Number of bootstrap replicates
    conf = 1 - (alpha / m) # Bonferroni correction to bootstrap confidence intervals
  )

  return(out)
}

parms <- fromJSON("parameter-values/json/params.json") %>%
  enframe() %>%
  unnest(value) %>%
  deframe()

# May take a couple of minutes to run the following two lines
out1 <- PRCC(parms, infection_number = 1)
out2 <- PRCC(parms, infection_number = 2)

# Plot the results
texlabels <- c(
  TeX("$\\alpha_{DV}$"),
  TeX("$\\B_{total}$"),
  TeX("$\\D_0$"),
  TeX("$\\delta_{BE}$"),
  TeX("$\\delta_{BV}$"),
  TeX("$\\delta_{E}$"),
  TeX("$\\D_T$"),
  TeX("$\\delta_{V}$"),
  TeX("$\\E_0$"),
  TeX("$\\kappa_{B}$"),
  TeX("$\\kappa_{E}$"),
  TeX("$\\kappa_{ED}$"),
  TeX("$\\kappa_{V}$"),
  TeX("$\\V_{0I}$"),
  TeX("$\\V_{0II}$")
)

p1 <- out1$PRCC %>%
  rownames_to_column(var = "Parameters") %>%
  rename(PRCC = original) %>%
  ggplot(aes(x = Parameters, y = PRCC)) +
  geom_hline(yintercept = 0) +
  geom_pointrange(aes(ymax = `max. c.i.`, ymin = `min. c.i.`), size = 0.7) +
  coord_cartesian(ylim = c(-1, 1)) +
  scale_x_discrete(labels = texlabels) +
  labs(x = NULL) +
  theme(axis.text.x = element_text(size = 14, color = "black")) +
  ggtitle("Max D(t) for first infection")

p2 <- out2$PRCC %>%
  rownames_to_column(var = "Parameters") %>%
  rename(PRCC = original) %>%
  ggplot(aes(x = Parameters, y = PRCC)) +
  geom_hline(yintercept = 0) +
  geom_pointrange(aes(ymax = `max. c.i.`, ymin = `min. c.i.`), size = 0.7) +
  coord_cartesian(ylim = c(-1, 1)) +
  scale_x_discrete(labels = texlabels) +
  theme(axis.text.x = element_text(size = 14, color = "black")) +
  ggtitle("Max D(t) for second infection")

plot_grid(p1, p2, labels = c("a)", "b)"), nrow = 2, label_size = 12)
ggsave(filename = "R/PRCC.png", width = 9.19, height = 5.75)
