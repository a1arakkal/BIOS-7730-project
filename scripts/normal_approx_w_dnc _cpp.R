# This script implements the normal approximation approach and implements
# the divide and conquer algorithm

###################
#### Libraries ####
###################

# devtools::install_git("https://github.com/a1arakkal/DNC")
library(dplyr)
library(bayestestR)
library(parallel)
library(DNC)
library(coda)

########################
#### Optim function ####
########################

optim_fun <- function(init, x, y, ratio, sigma, mu) {
  optim(par = init,
        fn = log_post_fun_dnc,
        method = "BFGS",
        control = list(fnscale= -1,
                       maxit = 1e6),
        hessian = T,
        sigma = sigma, 
        mu = mu,
        x = x,
        y = y,
        ratio = ratio)
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
  
  Opt <- optim_fun(init = rep(0, ncol(k_dat[, -1])),
                   x = k_dat[, -1], 
                   y = k_dat[, 1],
                   ratio = total_n/nrow(k_dat), #n/m_j
                   sigma = sigma,
                   mu = mu)
  
  params <- Opt$par
  names(params) <- colnames(k_dat[, -1])
  SigNew <- chol_inv(-Opt$hessian)
  
  ## Draw from multivariate normal
  MASS::mvrnorm(NN, mu =  params, Sigma = SigNew)
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
    c("inner_draws", "NN",
      "optim_fun", "total_n"),
    envir = environment()
  )
  clusterEvalQ(cl, {
    library(MASS)
    library(DNC)
  })
  
  results <- parLapply(
    cl,
    split_dat,
    function(x){inner_draws(k_dat = x, NN = NN, total_n = total_n)}
  )

  stopCluster(cl)
  gc()
  
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

  return(full_data_draws)
}

set.seed(123123)
full_data_draws <- par_fun(K=25)

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
            ESS = coda::effectiveSize(full_data_draws))

save(out, full_data_draws,
     file = "/Users/atlan/Adv_computing/AdvComp_final_proj/data/normal_approx_dnc_cpp.Rdata")
