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
mods <- files[grepl("csv", files) & !grepl("glm", files) & !(grepl("sim_data", files) | grepl("acc", files) | grepl("ess", files))]

for (i in mods){
  
  beta_draws <- read.csv(paste0(path, i))
  label <- ifelse(length(str_split_1(i, "_"))==3,
                  paste0(str_split_1(i, "_")[1:2], collapse ="_"),
                  ifelse(length(str_split_1(i, "_"))==4, 
                         paste0(str_split_1(i, "_")[1:3], collapse ="_"),
                  paste0(str_split_1(i, "_")[1:4], collapse ="_")))
  label <- paste0(label, "_Julia")
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
  

## Create key
key <- time %>% distinct(method) %>% 
  mutate(main = c("GLM", rep("MH", 4),
                  rep("Normal Approx.", 4),
                  "GLM", rep("MH", 2),
                  rep("Normal Approx.", 4)),
         program = c("R", rep(c("Rcpp", "R"), 4),
                     rep("Julia", 7)),
         opt_type = c(rep("", 5), rep("Finite Diff.", 4),
                      rep("", 3), rep(c("Auto Diff.", "Finite Diff."), 2)),
         DNC = c("", rep(c("D&C", "D&C","",""), 2),
                 "", "D&C", "", c("D&C", "D&C","",""))) %>% 
  mutate(nice = paste0(main, "-", program, ifelse(DNC == "", "", paste0(" (", DNC, ")")),
                       ifelse(opt_type == "", "", paste0(" (", opt_type, ")"))))


## ESS per sec
ess  %>%
  inner_join(time %>%  select(mean, method)) %>% 
  mutate(ESSS = ESS/(mean*60)) %>% 
  inner_join(key %>% mutate(DNC = ifelse(DNC == "", "non-D&C", DNC),
                            DNC = factor(DNC, levels = c("non-D&C", "D&C")),
                            main = factor(main, level = c("MH", "Normal Approx.", "GLM")))) %>% 
  mutate(group = fct_reorder(nice, ESSS),
         var = factor(Vars)) %>% 
  ggplot(aes(x=group, y = ESSS, fill = program)) + 
  geom_bar(stat = "identity") +
  theme(legend.title=element_blank()) +
  scale_fill_hue(c = 40) +
  facet_grid(main + DNC~var, scales = "free_y")+
  coord_flip() +
  geom_text(aes(label = nice(ESSS,2),hjust = -0.1))+
  ylab("Effective Sample Size per Second")+
  xlab("")+
  ylim(c(0,45000))


  
## Time comparison plot
time %>% inner_join(key) %>% 
  mutate(group = fct_reorder(nice, mean)) %>% 
  ggplot(aes(x=group, y = mean, fill = program)) + 
  geom_bar(stat = "identity") +
  theme(legend.title=element_blank()) +
  scale_fill_hue(c = 40) +
  coord_flip() +
  geom_text(aes(label = nice(mean),hjust = -0.1))+
  ylab("Computation Time (minutes)")+
  xlab("")

time %>%  inner_join(key %>% mutate(DNC = ifelse(DNC == "", "non-D&C", DNC),
                                    DNC = factor(DNC, levels = c("non-D&C", "D&C")),
                                    main = factor(main, level = c("MH", "Normal Approx.", "GLM")))) %>% 
  mutate(group = fct_reorder(nice, mean)) %>% 
  ggplot(aes(x=group, y = mean, fill = program)) + 
  geom_bar(stat = "identity") +
  theme(legend.title=element_blank()) +
  scale_fill_hue(c = 40) +
  facet_grid(main+DNC~., scales = "free_y") +
  coord_flip() +
  geom_text(aes(label = nice(mean),hjust = -0.1))+
  ylab("Computation Time (minutes)")+
  ylim(c(0, 510))+
  xlab("")


## Effect estimates comaprison plot
estimates %>% select(-nice) %>%  inner_join(key) %>% 
  inner_join(tibble(Vars = unique(estimates$Vars),
                    true = c(-3, 3.8, 1.1, 2.3, -0.2),
                    sym = paste0("beta[", 0:4, "]"))) %>% 
  inner_join(key) %>% 
  mutate(group = fct_reorder(nice, estiamte)) %>% 
  ggplot(aes(y = estiamte, x = group))+
  geom_point()+
  facet_wrap(~sym, scales = "free_x", labeller=label_parsed)+
  geom_errorbar(aes(ymin=lower_CI, ymax=upper_CI), colour="black", width=.1) + 
  coord_flip()+
  geom_hline(aes(yintercept = true),linetype="dashed", color = "red") +
  ylab("Estimate (95% CI)")+
  xlab("")+
  ggtitle(expression(beta[0]))+
  theme(plot.title = element_text(hjust = 0.5))
