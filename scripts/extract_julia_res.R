library(tidyverse)
library(bayestestR)
library(coda)

path <- "/Volumes/argon_home/Adv_computing/AdvComp_final_proj/data/"

# beta_draws <- read_csv(paste0(path, "normal_approx_full_136.888.csv"))
# beta_draws <- read_csv(paste0(path, "normal_approx_full_FD_1694.0204.csv"))
# beta_draws <- read_csv(paste0(path, "normal_approx_dnc_5.7441.csv"))
# beta_draws <- read_csv(paste0(path, "normal_approx_dnc_FD_82.5457.csv"))

summary <- describe_posterior(as.data.frame(beta_draws))
hpd <- hdi(as.data.frame(beta_draws))

results <- tibble(var = summary[,1],
                  coef = summary[,2],
                  central_lower = summary[, 4],
                  central_upper= summary[, 5],
                  hpd_lower = hpd[,3],
                  hdp_upper = hpd[,4])

coda::effectiveSize(as.matrix(beta_draws))
