

###################
#### Libraries ####
###################

library(dplyr)

#######################
#### Main Function ####
#######################

run_glm <- function(){
  
  ###################
  #### Load Data ####
  ###################
  
  mod_dat <- data.table::fread("/Users/atlan/Adv_computing/AdvComp_final_proj/data/sim_data.csv",
                               colClasses = "numeric",
                               showProgress = T) 

  mod <- glm(y ~ ., data = mod_dat, 
             family = "binomial")
  
  return(mod)
}

glm_res <- run_glm()

time <- microbenchmark::microbenchmark(run_glm(),
                                       unit = "s",
                                       times = 10)

elapsed <- tibble(summary(time)) %>% 
  dplyr::select(-expr, -neval) %>% 
  mutate(across(everything(), ~./60)) #convert to minutes

out <- list(result_table = summary(glm_res)$coefficients,
            comp_time = elapsed)

save(out,
     file = "/Users/atlan/Adv_computing/AdvComp_final_proj/data/glm.Rdata")


