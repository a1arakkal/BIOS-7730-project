# Advanced Computing Final Project

This repository contains code for the final project for BIOS 7730.

The repository is structed as follows:


# `jobs`
  - contains the job text files used for HPC job submisson 

# `scripts`
 - ## `Julia`
   -  ### `MH`
      - `MH_w_dnc.jl`: Julia script implementing Metropolis–Hastings algorithm along with Divide and Conquer (D&C)-recentering approach.
      - `MH_wo_dnc.jl`: Julia script implementing Metropolis–Hastings algorithm along without D&C.
   -  ### `normal_approx`
      - `normal_approx_w_dnc.jl`: Julia script implementing normal appoximation approach along with D&C-recentering approach. Here foward automatic differentiation is used in the optimization procedures.
      - `normal_approx_wo_dnc.jl`: Julia script implementing normal appoximation approach along without D&C. Here foward automatic differentiation is used in the optimization procedures.
      - `normal_approx_w_dnc_FD.jl`: Julia script implementing normal appoximation approach along with D&C-recentering approach. Here finite difference is used in the optimization procedures.
      - `normal_approx_wo_dnc_FD.jl`: Julia script implementing normal appoximation approach along without D&C. Here finite difference is used in the optimization procedures.
   -  ### `glm.jl`
      - Julia script to run simple logistic regression model on the full data.
 - ## `R`
   -  ### `MH`
      - `MH_w_dnc.R`: R script implementing Metropolis–Hastings algorithm along with D&C-recentering approach.
      - `MH_wo_dnc.R`: R script implementing Metropolis–Hastings algorithm along without D&C.
   -  ### `normal_approx`
      - `normal_approx_w_dnc.R`: R script implementing normal appoximation approach along with D&C-recentering approach. Here finite difference is used in the optimization procedures.
      - `normal_approx_wo_dnc.R`: R script implementing normal appoximation approach along without D&C. Here finite difference is used in the optimization procedures.
    -  ### `glm.R`
       - R script to run simple logistic regression model on the full data.
    - ### `sim_data.R`
      - R script to generate the simulated logistic regression data used in the various analyses.
    -  ### `extract_results.R`
        - R script to extract results from all methods employed in the various analyses. Also contains scripts to plot comparisons of results. 
    
 - ## `Rcpp`
   - Note all Rcpp functions used in the various analyses can be found in the DNC package in the following GitHub repository: https://github.com/a1arakkal/DNC/
   -  ### `MH`
      - `MH_w_dnc_cpp.R`: R script with Rcpp integration implementing Metropolis–Hastings algorithm along with D&C-recentering approach.
      - `MH_wo_dnc_cpp.R`: R script with Rcpp integration implementing Metropolis–Hastings algorithm along without D&C.
   -  ### `normal_approx`
      - `normal_approx_w_dnc_cpp.R`: R script with Rcpp integration implementing normal appoximation approach along with D&C-recentering approach. Here finite difference is used in the optimization procedures.
      - `normal_approx_wo_dnc_cpp.R`: R script with Rcpp integration implementing normal appoximation approach along without D&C. Here finite difference is used in the optimization procedures.

