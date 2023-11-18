

using CSV, DataFrames, Statistics, Distributions, LinearAlgebra, Optim, NLSolversBase, BenchmarkTools, Random, Distributed, MLUtils, GLM

Threads.nthreads()

###################
#### Functions ####
###################

function log_post_fun(param, x, y, ratio, sigma, mu)
    p = 1 ./ (1 .+ exp.(-x * param))
    loglike = sum(logpdf.(Binomial.(1, p), y))
    logprior = logpdf(MvNormal(vec(mu), sigma), param)[1,1]
    return (ratio*loglike) + logprior
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
  x = [ones(size(mod_dat)[1]) mod_dat[:, 2:size(mod_dat)[2]]]
  ratio = total_n/size(mod_dat)[1]

  ## Run MH
  acc_count = 0
  sd_prop = 0.02

  fm = @formula(x1 ~ x2 + x3 + x4+ x5)
  logit = glm(fm, DataFrame(mod_dat), Binomial(), LogitLink())
  temp = GLM.coeftable(logit)
  proposal_cov = convert(Matrix, Diagonal(temp.cols[2].^2)*sd_prop)
 
  beta_draws_mh = zeros(NN, size(mod_dat)[2])
  beta_draws_mh[1,:] = temp.cols[1]

  for i in 2:NN
    
    beta_draws_mh[i,:] = beta_draws_mh[i-1,:]

    ## Draw beta
    beta_proposal = rand(MvNormal(beta_draws_mh[i,:], proposal_cov), 1)

    acc_prob = exp(log_post_fun(beta_proposal, x, y, ratio, sigma, mu) - log_post_fun(beta_draws_mh[i,:], x, y, ratio, sigma, mu))
    
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

@time beta_draws, acc_prob  = run_MH(100, 10)
Random.seed!(123123)
time = @elapsed beta_draws, acc_prob = run_MH(10000, 1000)
time = round(time, digits = 4)
CSV.write("/Users/atlan/Adv_computing/AdvComp_final_proj/data/MH_full_$time.csv", 
          DataFrame(Matrix(beta_draws)))
CSV.write("/Users/atlan/Adv_computing/AdvComp_final_proj/data/MH_full_acc_prob_$time.csv", 
          DataFrame(reshape([acc_prob], length([acc_prob]), 1)))

