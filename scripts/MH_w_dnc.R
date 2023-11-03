# This script implements the normal approximation approach and implements
# the divide and conquer algorithm

###################
#### Libraries ####
###################

library(dplyr)
library(bayestestR)
library(parallel)
library(coda)

###################
#### Functions ####
###################

log_post_fun <- function(param, x, y, ratio, sigma, mu) {
  (ratio)*sum(y*x%*%param - log(1 + exp(x%*%param))) -
    drop(0.5*t((param - mu)) %*% solve(sigma) %*% (param - mu))
}

#######################
#### Main function ####
#######################

inner_draws <- function(k_dat, NN = 1e4, total_n) {
  
  k_dat <- as.matrix(k_dat)
  
  # Set prior values
  mu <- rep(0, ncol(k_dat)) # need to add 1 for intercept and subtract one for y so it cancels out
  sd_temp <- sqrt(diag(var(k_dat[, -1])))
  scale <- sd_temp # scale continuous vars
  sigma <- diag(c(40^2, 3^2 * scale))  # flat priors on beta

  k_dat <- cbind(k_dat[,1], 1, k_dat[, 2:ncol(k_dat)])
  colnames(k_dat)[1] <- "y" 
  colnames(k_dat)[2] <- "intercept" 
  
  acc_count <- 0
  sd_prop <- 0.02
  glm_fit <- glm(y ~ -1 + ., as.data.frame(k_dat), family = "binomial")
  proposal_cov <- diag(diag(vcov(glm_fit)))*sd_prop
  
  beta_draws_mh <- matrix(0.0, NN, ncol(k_dat) -1)
  beta_draws_mh[1, ] <- coef(glm_fit)
  
  ################
  #### Run MH ####
  ################
  
  for(i in 2:NN) {
    
    beta_draws_mh[i,] <- beta_draws_mh[i - 1,]
    
    ## Draw beta
    beta_proposal <- MASS::mvrnorm(1, mu = beta_draws_mh[i,], Sigma = proposal_cov)
    
    acc_prob <- exp(log_post_fun(beta_proposal, 
                                 x = k_dat[, -1],
                                 y = k_dat[, 1], 
                                 ratio = total_n/nrow(k_dat), #n/m_j
                                 sigma = sigma,
                                 mu = mu) - 
                      log_post_fun(beta_draws_mh[i, ],
                                   x = k_dat[, -1],
                                   y = k_dat[, 1], 
                                   ratio = total_n/nrow(k_dat), #n/m_j
                                   sigma = sigma,
                                   mu = mu))
    
    if (runif(1) < acc_prob) {
      
      beta_draws_mh[i,] <- beta_proposal
      acc_count <- acc_count + 1
      
    }
    
  }
  
  return(list(beta_draws_mh  = beta_draws_mh,
              acc_prob = acc_count/(NN-1)))

}

############################
#### Parallel function #####
############################

par_fun <- function(K, NN = 1e4, burn_in_n = 1000){
  
  if(NN <= burn_in_n){
    stop("Discarding more than sampled")
  }
  
  ###################
  #### Load Data ####
  ###################
  
  mod_dat <- data.table::fread("/Users/atlan/Adv_computing/AdvComp_final_proj/data/sim_data.csv",
                               colClasses = "numeric",
                               showProgress = T) 
  
  mod_dat <- as.data.frame(mod_dat)
  total_n <- nrow(mod_dat)
  
  #########################################################
  #### Partition data into batches and set up clusters ####
  #########################################################
  
  split_dat <- split(mod_dat, rep(1:K, length.out = nrow(mod_dat), each = ceiling(nrow(mod_dat)/K)))
  
  cl <- makeCluster(min(detectCores(), K))
  clusterExport(
    cl,
    c("inner_draws", 
      "log_post_fun",
      "NN",
      "total_n"),
    envir = environment()
  )
  clusterEvalQ(cl, {
    library(MASS)
  })
  
  results_temp <- parLapply(
    cl,
    split_dat,
    function(x){inner_draws(k_dat = x, NN = NN, total_n = total_n)}
  )

  stopCluster(cl)
  gc()
  
  acc_probs <- sapply(1:length(results_temp), function(i) results_temp[[i]]$acc_prob)
  results <- lapply(1:length(results_temp), function(i) results_temp[[i]]$beta_draws_mh)
  
  ########################
  #### Remove burn-in ####
  ########################

  remove_burnin <- function(x, burnin) {
    x[-(1:burnin),]
  }
  
  results <- lapply(results, remove_burnin, burn_in_n)
  
  ########################
  #### Recenter Draws ####
  ########################
  
  subset_mean <- t(sapply(1:K, function(i) colMeans(results[[i]]))) # rows indicate subset
  global_mean <- colMeans(bind_rows(lapply(results, data.frame)))
  recenter <- t(sapply(1:K, function(i) subset_mean[i, ] - global_mean)) # rows indicate subset
  results_recentered <- lapply(1:K,function(i) results[[i]] - matrix(recenter[i, ],
                                                                     nrow = nrow(results[[i]]),
                                                                     ncol = ncol(results[[i]]),
                                                                     byrow = T))
  full_data_draws <- do.call(rbind, results_recentered)

  return(list(full_data_draws = full_data_draws,
              acc_probs = acc_probs))
}

set.seed(123123)

out <- par_fun(K=25)
full_data_draws <- out$full_data_draws
acc_probs <- out$acc_probs

gc()

time <- microbenchmark::microbenchmark(par_fun(K=25),
                                       unit = "s",
                                       times = 10)

elapsed <- tibble(summary(time)) %>% 
  dplyr::select(-expr, -neval) %>% 
  mutate(across(everything(), ~./60)) #convert to minutes

#################
#### Summary ####
#################

summary <- describe_posterior(as.data.frame(full_data_draws))
hpd <- hdi(as.data.frame(full_data_draws))

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
            acc_probs = acc_probs,
            ESS = coda::effectiveSize(full_data_draws))

save(out, full_data_draws,
     file = "/Users/atlan/Adv_computing/AdvComp_final_proj/data/MH_dnc.Rdata")
