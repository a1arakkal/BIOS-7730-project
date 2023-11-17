

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
#### Main function ####
#######################

function inner_draws(k_dat, NN, total_n, burn_in)
  
  # Set prior values
  mu = zeros(size(k_dat)[2], 1)
  temp_sigma = Diagonal([40^2; 3^2 .* var.(eachcol(k_dat[:, 1:end .!= 1]))])  # flat priors on beta
  sigma = convert(Matrix,  temp_sigma)

  # Add vector of 1s for intercept 
  y = k_dat[:, 1]
  x = [ones(size(k_dat)[1]) k_dat[:, 2:size(k_dat)[2]]]
  ratio = total_n/size(k_dat)[1]

  ## Run MH
  acc_count = 0
  sd_prop = 0.02

  fm = @formula(x1 ~ x2 + x3 + x4+ x5)
  logit = glm(fm, DataFrame(k_dat), Binomial(), LogitLink())
  temp = GLM.coeftable(logit)
  proposal_cov = convert(Matrix, Diagonal(temp.cols[2].^2)*sd_prop)
 
  beta_draws_mh = zeros(NN, size(k_dat)[2])
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

#######################
#### Main Function ####
#######################

function par_fun(NN = 10000, burn_in = 1000)
 
  if (burn_in > NN)
    return "Error: throwing out more than sampled"
  end

  ###################
  #### Load Data ####
  ###################
  
  mod_dat = CSV.read("/Users/atlan/Adv_computing/AdvComp_final_proj/data/sim_data.csv",
                     DataFrame, header = 1)

  total_n = size(mod_dat)[1]

  split_dat = chunk(Matrix(mod_dat), 25, dims=1)

  res = Array{Float64}(undef,NN-burn_in,size(mod_dat)[2],25)
  acc_count = Array{Float64}(undef,25,1)

  Threads.@threads for i in 1:25
    draws, acc_prob = inner_draws(split_dat[i], NN, total_n, burn_in)
    (@view acc_count[i]) .= acc_prob
    (@view res[:,:,i]) .= draws
  end
  
  ########################
  #### Recenter Draws ####
  ########################
  
  subset_mean = [mean(res[:,j,i]) for i in 1:size(res)[3], j in 1:size(res)[2]]
  weights =  [size(res[:,:, i])[1] for i in 1:size(res)[3]]
  weight_mat = reshape(weights, 1, :)/sum(weights)
  global_mean = vec(*(weight_mat, subset_mean))'
  
  for i in 1:size(res)[3]
    res[:,:,i] = res[:,:,i] .- subset_mean[i,:]' .+ global_mean
  end

  full_data_draws = reduce(vcat,(eachslice(res, dims=3)))

  return full_data_draws, acc_count

end

@time beta_draws, acc_prob  = par_fun(100, 10)
Random.seed!(123123)
time = @elapsed beta_draws, acc_prob = par_fun(10000, 1000)
time = round(time, digits = 4)
CSV.write("/Users/atlan/Adv_computing/AdvComp_final_proj/data/MH_dnc_$time.csv", 
          DataFrame(Matrix(beta_draws)))
CSV.write("/Users/atlan/Adv_computing/AdvComp_final_proj/data/MH_dnc_acc_prob_$time.csv", 
          DataFrame(Matrix(acc_prob)))

