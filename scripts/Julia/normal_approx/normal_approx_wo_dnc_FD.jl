

using CSV, DataFrames, Statistics, Distributions, LinearAlgebra, Optim, NLSolversBase, BenchmarkTools, Random, MLUtils

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

function run_norm(NN = 10000)
 
  ###################
  #### Load Data ####
  ###################
  
  mod_dat = CSV.read("/Users/atlan/Adv_computing/AdvComp_final_proj/data/sim_data.csv",
                     DataFrame, header = 1)

  ########################################################
  #### Set prior values / set up results df ##############
  ########################################################

  mu = zeros(size(mod_dat)[2], 1)
  temp_sigma = Diagonal([40^2; 3^2 .* var.(eachcol(mod_dat[:, 1:end .!= 1]))])  # flat priors on beta
  sigma = convert(Matrix,  temp_sigma)
  
  ##################################################
  #### Add vector of 1s for intercept ##############
  ##################################################

  insert!(mod_dat, 2, vec(ones(size(mod_dat)[1], 1)), :intercept)
  mod_dat_mat = convert(Matrix,  mod_dat)
  y = mod_dat_mat[:, 1]
  x = mod_dat_mat[:, 1:end .!= 1]

  ###########################
  #### Run Normal Approx ####
  ###########################

  obj = TwiceDifferentiable(theta -> -1.0*log_post_fun(theta, x, y, sigma, mu), vec(zeros(size(mod_dat)[2]-1, 1)))

  Opt = optimize(obj, vec(zeros(size(mod_dat)[2]-1, 1)), BFGS(),
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

beta_draws = run_norm(10000)
Random.seed!(123123)
time = @elapsed beta_draws = run_norm(10000)
time = round(time, digits = 4)
CSV.write("/Users/atlan/Adv_computing/AdvComp_final_proj/data/normal_approx_full_FD_$time.csv", 
          DataFrame(Matrix(beta_draws)))