library(tidyverse)
library(bayestestR)
library(coda)

nice <- function(x, digs=3){
  format(round(x, digs), nsmall = digs)
}

path <- "/Volumes/argon_home/Adv_computing/AdvComp_final_proj/data/"
estimates <- tibble()
time <- tibble()
ess <- tibble()

## R glm
load(paste0(path, "glm.RData"))
temp <- out$result_table
est <- tibble(Vars = c("(Intercept)", paste0("X", 1:4)),
              estiamte = temp[,1],
              lower_CI =  temp[,1] - qnorm(0.975)*temp[,2],
              upper_CI =  temp[,1] + qnorm(0.975)*temp[,2],
              nice = paste0(nice(estiamte), " (", nice(lower_CI), ", ", nice(upper_CI), ")")) %>% 
  mutate(method = "glm_R") %>% 
  dplyr::select(Vars, estiamte, lower_CI, upper_CI, nice, method)

estimates <- bind_rows(estimates, est)
time <- bind_rows(time, out$comp_time %>% mutate(method = "glm_R"))


## All other R
files <- list.files(path)
mods <- files[grepl("Rdata", files) & !grepl("glm", files)]

for (i in mods){
  
  load(paste0(path, i))
  label <- str_remove(i, ".Rdata")
  temp <- out$result_table
  est <- tibble(Vars = c("(Intercept)", paste0("X", 1:4)),
                estiamte = temp$coef,
                lower_CI = temp$central_lower,
                upper_CI = temp$central_upper,
                nice = paste0(nice(estiamte), " (", nice(lower_CI), ", ", nice(upper_CI), ")")) %>% 
    mutate(method = label) %>% 
    dplyr::select(Vars, estiamte, nice, lower_CI, upper_CI, method)
  estimates <- bind_rows(estimates, est)
  time <- bind_rows(time, out$comp_time %>% mutate(method = label))
  temp_ess <-tibble(Vars = c("(Intercept)", paste0("X", 1:4)),
                    ESS = out$ESS,
                    method = label)
  ess <- bind_rows(ess, temp_ess)
  
}

## Julia glm
mods <- files[grepl("csv", files) & grepl("glm", files)]
temp <- read.csv(paste0(path, mods))
est <- tibble(Vars = c("(Intercept)", paste0("X", 1:4)),
              estiamte = temp$Estimate,
              lower_CI = temp$lowerCI,
              upper_CI = temp$upperCI,
              nice = paste0(nice(estiamte), " (", nice(lower_CI), ", ", nice(upper_CI), ")")) %>% 
  mutate(method = "glm_Julia") %>% 
  dplyr::select(Vars, estiamte, nice, lower_CI, upper_CI, method)
estimates <- bind_rows(estimates, est)

temp_time <- as.numeric(str_remove(str_split_1(mods, "_")[length(str_split_1(mods, "_"))], ".csv"))/60
time <- bind_rows(time, time %>% filter(row_number()==1) %>% 
                    mutate_at(vars(everything()), ~NA) %>% 
                    mutate(mean = temp_time) %>% 
                    mutate(method = "glm_Julia"))

## All other Julia
mods <- files[grepl("csv", files) & !grepl("glm", files) & !(grepl("sim_data", files) | grepl("acc", files))]

for (i in mods){
  
  temp <- read.csv(paste0(path, i))
  label <- ifelse(length(str_split_1(i, "_"))==3,
                  paste0(str_split_1(i, "_")[1:2], collapse ="_"),
                  ifelse(length(str_split_1(i, "_"))==4, 
                         paste0(str_split_1(i, "_")[1:3], collapse ="_"),
                  paste0(str_split_1(i, "_")[1:4], collapse ="_")))
  
  summary <- describe_posterior(as.data.frame(beta_draws))
  hpd <- hdi(as.data.frame(beta_draws))
  temp <- tibble(var = summary[,1],
                 coef = summary[,2],
                 central_lower = summary[, 4],
                 central_upper= summary[, 5],
                 hpd_lower = hpd[,3],
                 hdp_upper = hpd[,4])
  
  est <- tibble(Vars = c("(Intercept)", paste0("X", 1:4)),
                estiamte = temp$coef,
                lower_CI = temp$central_lower,
                upper_CI = temp$central_upper,
                nice = paste0(nice(estiamte), " (", nice(lower_CI), ", ", nice(upper_CI), ")")) %>% 
    mutate(method = label) %>% 
    dplyr::select(Vars, estiamte, nice, lower_CI, upper_CI, method)
  estimates <- bind_rows(estimates, est)
  temp_time <- as.numeric(str_remove(str_split_1(i, "_")[length(str_split_1(i, "_"))], ".csv"))/60
  time <- bind_rows(time, time %>% filter(row_number()==1) %>% 
                      mutate_at(vars(everything()), ~NA) %>% 
                      mutate(mean = temp_time) %>% 
                      mutate(method = label))
  
  temp_ess <-tibble(Vars = c("(Intercept)", paste0("X", 1:4)),
                    ESS =  coda::effectiveSize(as.matrix(beta_draws)),
                    method = label)
  
  ess <- bind_rows(ess, temp_ess)
  
}

