# This script creates the simulated logistic regression data that is used in
# the analysis.

simdata <- function(n, beta, rho = 0.1) {
  sigma <- matrix(rho, length(beta)-1, length(beta)-1) 
  diag(sigma) <- 1
  X <- mvtnorm::rmvnorm(n, sigma = sigma)
  eta <- cbind(1,X) %*% beta
  probs <- 1 / (exp(-eta) + 1) 
  data.frame(
    y = as.integer(rbinom(n, 1, probs)),
    X )
}


beta_true = c(-3, 3.8, 1.1, 2.3, -0.2)
out <- "/Users/atlan/Adv_computing/AdvComp_final_proj/data/sim_data.csv"
file.create(out) # Create output csv file
out_name <-  c("y", paste0("X", 1:(length(beta_true)-1)))
data.table::fwrite(as.data.frame(t(out_name)), file = out, col.names = F,
                   append = T, showProgress = F)

set.seed(123123)
for (i in 1:10){
  temp <- simdata(n = 1e6, beta = beta_true, rho = 0.1)
  data.table::fwrite(temp, file = out, col.names = F,
                     append = T, showProgress = F)
}
  


