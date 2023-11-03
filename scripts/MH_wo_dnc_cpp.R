# This script implements the normal approximation approach with out implementing 
# the divide and conquer algorithm

###################
#### Libraries ####
###################

# devtools::install_git("https://github.com/a1arakkal/DNC")
library(dplyr)
library(bayestestR)
library(DNC)
library(MASS)
library(coda)

#######################
#### Main Function ####
#######################

run_MH <- function(NN = 1e4){
  
  ###################
  #### Load Data ####
  ###################
  
  mod_dat <- data.table::fread("/Users/atlan/Adv_computing/AdvComp_final_proj/data/sim_data.csv",
                               colClasses = "numeric",
                               showProgress = T) 
  
  mod_dat <- as.matrix(mod_dat)
  
  ########################################################
  #### Set prior values / set up results df ##############
  ########################################################
  
  mu <- rep(0, ncol(mod_dat)) # need to add 1 for intercept and subtract one for y so it cancels out
  sd_temp <- sqrt(diag(var(mod_dat[, -1])))
  scale <- sd_temp # scale continuous vars
  sigma <- diag(c(40^2, 3^2 * scale))  # flat priors on beta
  
  ##################################################
  #### Add vector of 1s for intercept ##############
  ##################################################
  
  mod_dat <- cbind(mod_dat[,1], 1, mod_dat[, 2:ncol(mod_dat)])
  colnames(mod_dat)[1] <- "y" 
  colnames(mod_dat)[2] <- "intercept" 
  
  ######################################################################
  #### Proposal variance/covaraince matrix and output sample matrix ####
  ######################################################################
  
  sd_prop <- 0.25
  glm_fit <- glm(y ~ -1 + ., as.data.frame(mod_dat), family = "binomial")
  proposal_cov <- diag(diag(vcov(glm_fit)))*sd_prop
  
  beta_draws_mh <- matrix(0.0, NN, ncol(mod_dat) -1)
  beta_draws_mh[1, ] <- coef(glm_fit)
  
  ################
  #### Run MH ####
  ################
  
  out <- MH(MH_draws = beta_draws_mh, 
            proposal_cov = proposal_cov,
            x = mod_dat[, -1],
            y = mod_dat[, 1],
            mu = mu,
            sigma = sigma,
            ratio = 1) # 1 means no DNC
  
  #remove burn-in
  beta_draws_mh <- beta_draws_mh[-(1:1000), ]
  
  return(list(beta_draws_mh  = beta_draws_mh,
              acc_count = out,
              acc_prob = out/(NN-1)))
}

set.seed(123123)
out <- run_MH(NN = 1e4)
beta_draws <- out$beta_draws_mh
acc_prob <- out$acc_prob

time <- microbenchmark::microbenchmark(run_norm(NN = 1e4),
                                       unit = "s",
                                       times = 10)

elapsed <- tibble(summary(time)) %>% 
  dplyr::select(-expr, -neval) %>% 
  mutate(across(everything(), ~./60)) #convert to minutes

#################
#### Summary ####
#################


summary <- describe_posterior(as.data.frame(beta_draws))
hpd <- hdi(as.data.frame(beta_draws))

results <- tibble(var = summary[,1],
                  coef = summary[,2],
                  central_lower = summary[, 4],
                  central_upper= summary[, 5],
                  hpd_lower = hpd[,3],
                  hdp_upper = hpd[,4])

##############
#### Save ####
##############

out <- list(result_table = results,
            comp_time = elapsed,
            acc_prob = acc_prob,
            ESS = coda::effectiveSize(full_data_draws))

save(out, beta_draws,
     file = "/Users/atlan/Adv_computing/AdvComp_final_proj/data/MH_full_cpp.Rdata")
