

using CSV, DataFrames, Statistics, Distributions, LinearAlgebra, Optim, NLSolversBase, BenchmarkTools, Random, Distributed, MLUtils, GLM

Threads.nthreads()

###################
#### Functions ####
###################

function log_post_fun(param, x, y, sigma, mu)
    p = 1 ./ (1 .+ exp.(-x * param))
    loglike = sum(logpdf.(Binomial.(1, p), y))
    logprior = logpdf(MvNormal(vec(mu), sigma), param)[1,1]
    return loglike + logprior
end
       
#######################
#### Main Function ####
#######################

function run_MH(NN = 10000, burn_in = 1000)
 
  if (burn_in > NN)
    return "Error: throwing out more than sampled"
  end

  # Load Data   
  mod_dat = CSV.read("/Users/atlan/Adv_computing/AdvComp_final_proj/data/sim_data.csv",
                     DataFrame, header = 1)
  
  # Set prior values
  mu = zeros(size(mod_dat)[2], 1)
  temp_sigma = Diagonal([40^2; 3^2 .* var.(eachcol(mod_dat[:, 1:end .!= 1]))])  # flat priors on beta
  sigma = convert(Matrix,  temp_sigma)

  # Add vector of 1s for intercept 
  y = mod_dat[:, 1]
  x = Matrix([ones(size(mod_dat)[1]) mod_dat[:, 2:size(mod_dat)[2]]])

  ## Run MH
  acc_count = 0
  sd_prop = 0.02

  fm = @formula(y ~ X1 + X2 + X3 + X4)
  logit = glm(fm, DataFrame(mod_dat), Binomial(), LogitLink())
  temp = GLM.coeftable(logit)
  proposal_cov = convert(Matrix, Diagonal(temp.cols[2].^2)*sd_prop)
 
  beta_draws_mh = zeros(NN, size(mod_dat)[2])
  beta_draws_mh[1,:] = temp.cols[1]

  for i in 2:NN
    
    beta_draws_mh[i,:] = beta_draws_mh[i-1,:]

    ## Draw beta
    beta_proposal = rand(MvNormal(beta_draws_mh[i,:], proposal_cov), 1)

    acc_prob = exp(log_post_fun(beta_proposal, x, y, sigma, mu) - log_post_fun(beta_draws_mh[i,:], x, y, sigma, mu))
    
    if (rand(Uniform(0,1), 1)[1] < acc_prob)
      beta_draws_mh[i,:] = beta_proposal
      acc_count += 1
    end

  end  

  acc_prob = acc_count/NN

  #remove burn-in
  beta_draws_mh = beta_draws_mh[(burn_in+1):end, :]

  return beta_draws_mh, acc_prob

end

time_run = @elapsed beta_draws, acc_prob  = run_MH(100, 10)
Random.seed!(123123)
time_run = @elapsed beta_draws, acc_prob = run_MH(10000, 1000)
time_run = round(time_run, digits = 4)
CSV.write("/Users/atlan/Adv_computing/AdvComp_final_proj/data/MH_full_$time_run.csv", 
          DataFrame(beta_draws, :auto))
CSV.write("/Users/atlan/Adv_computing/AdvComp_final_proj/data/MH_full_acc_prob_$time_run.csv", 
          DataFrame(reshape([acc_prob], length([acc_prob]), 1), :auto))

