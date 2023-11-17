# This script implements the normal approximation approach with out implementing 
# the divide and conquer algorithm

###################
#### Libraries ####
###################

# devtools::install_git("https://github.com/a1arakkal/DNC")
library(dplyr)
library(bayestestR)
library(DNC)
library(coda)

#######################
#### Main Function ####
#######################

run_norm <- function(NN = 1e4){
  
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
  
  ###########################
  #### Run Normal Approx ####
  ###########################
  
  Opt <- optim(par = rep(0,ncol(mod_dat[, -1])),
               fn = log_post_fun, 
               method = "BFGS",
               control = list(fnscale= -1,
                              maxit = 1e6),
               hessian = T,
               x = mod_dat[, -1], 
               y = mod_dat[, 1],
               mu = mu,
               sigma = sigma)
  
  params <- Opt$par
  names(params) <- colnames(mod_dat[, -1])
  SigNew <-  chol_inv(-Opt$hessian)
  ## Draw from multivariate normal
  beta_draws <- MASS::mvrnorm(NN, mu = params, Sigma = SigNew)
  #remove burn-in
  beta_draws <- beta_draws[-(1:1000), ]
  return(beta_draws)
}

set.seed(123123)
beta_draws <- run_norm(NN = 1e4)

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
            ESS = coda::effectiveSize(beta_draws))

save(out, beta_draws,
     file = "/Users/atlan/Adv_computing/AdvComp_final_proj/data/normal_approx_full_cpp.Rdata")
