# This script implements the normal approximation approach with out implementing 
# the divide and conquer algorithm

###################
#### Libraries ####
###################

library(dplyr)
library(bayestestR)
library(coda)

###################
#### Functions ####
###################

log_post_fun <- function(param, x, y, sigma, mu) {
  
  sum(y*x%*%param - log(1 + exp(x%*%param))) -
    drop(0.5*t((param - mu)) %*% solve(sigma) %*% (param - mu))
  
}

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
  
  acc_count <- 0
  sd_prop <- 0.25
  glm_fit <- glm(y ~ -1 + ., as.data.frame(mod_dat), family = "binomial")
  proposal_cov <- diag(diag(vcov(glm_fit)))*sd_prop
  
  beta_draws_mh <- matrix(0, NN, ncol(mod_dat) -1)
  beta_draws_mh[1, ] <- coef(glm_fit)
  
  ################
  #### Run MH ####
  ################
  
  for(i in 2:NN) {
    
    beta_draws_mh[i,] <- beta_draws_mh[i - 1,]
    
    ## Draw beta
    beta_proposal <- MASS::mvrnorm(1, mu = beta_draws_mh[i,], Sigma = proposal_cov)
    
    acc_prob <- exp(log_post_fun(beta_proposal, 
                                 x = mod_dat[, -1],
                                 y = mod_dat[, 1], 
                                 sigma = sigma,
                                 mu = mu) - 
                      log_post_fun(beta_draws_mh[i, ],
                                   x = mod_dat[, -1],
                                   y = mod_dat[, 1], 
                                   sigma = sigma,
                                   mu = mu))
    
    if (runif(1) < acc_prob) {
      
      beta_draws_mh[i,] <- beta_proposal
      acc_count <- acc_count + 1
      
    }
    
  }
  #remove burn-in
  beta_draws_mh <- beta_draws_mh[-(1:1000), ]
  
  return(list(beta_draws_mh  = beta_draws_mh,
              acc_count = acc_count,
              acc_prob = acc_count/(NN-1)))
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
     file = "/Users/atlan/Adv_computing/AdvComp_final_proj/data/MH_full.Rdata")
