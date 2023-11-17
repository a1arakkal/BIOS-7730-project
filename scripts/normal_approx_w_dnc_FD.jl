

using CSV, DataFrames, Statistics, Distributions, LinearAlgebra, Optim, NLSolversBase, BenchmarkTools, Random, Distributed, MLUtils

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

function inner_draws(k_dat, NN, total_n)
  
  # Set prior values
  mu = zeros(size(k_dat)[2], 1)
  temp_sigma = Diagonal([40^2; 3^2 .* var.(eachcol(k_dat[:, 1:end .!= 1]))])  # flat priors on beta
  sigma = convert(Matrix,  temp_sigma)

  # Add vector of 1s for intercept 
  y = k_dat[:, 1]
  x = [ones(size(k_dat)[1]) k_dat[:, 2:size(k_dat)[2]]]

  ## Run Normal Approx 
  ratio = total_n/size(k_dat)[1]
  obj = TwiceDifferentiable(theta -> -1.0*log_post_fun(theta, x, y, ratio, sigma, mu), vec(zeros(size(x)[2], 1)))

  Opt = optimize(obj, vec(zeros(size(x)[2], 1)), BFGS(),
                 Optim.Options(x_tol = 1e-8,
                              iterations = 1000000))

  params = Optim.minimizer(Opt)
  SigNew = inv(NLSolversBase.hessian!(obj, params))

  ## Draw from multivariate normal
  beta_draws = rand(MvNormal(params, Symmetric(SigNew)), NN)'

  #remove burn-in
  beta_draws = beta_draws[1001:end, :]

  return beta_draws

end

#######################
#### Main Function ####
#######################

function par_fun(NN = 10000)
 
  ###################
  #### Load Data ####
  ###################
  
  mod_dat = CSV.read("/Users/atlan/Adv_computing/AdvComp_final_proj/data/sim_data.csv",
                     DataFrame, header = 1)

  total_n = size(mod_dat)[1]

  split_dat = chunk(Matrix(mod_dat), 25, dims=1)

  res = Array{Float64}(undef,NN-1000,size(mod_dat)[2],25)

  Threads.@threads for i in 1:25
    (@view res[:,:,i]) .= inner_draws(split_dat[i], NN, total_n)
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

  return full_data_draws

end

beta_draws = par_fun(10000)
Random.seed!(123123)
time = @elapsed beta_draws = par_fun(10000)
time = round(time, digits = 4)
CSV.write("/Users/atlan/Adv_computing/AdvComp_final_proj/data/normal_approx_dnc_FD_$time.csv", 
          DataFrame(Matrix(beta_draws)))

