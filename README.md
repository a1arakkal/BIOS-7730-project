# Advanced Computing Final Project

This repository contains code for the final project for BIOS 7730.

The repository is structed as follows:

# `jobs`
  - contains the job text files used for HPC job submisson 

# `scripts`
 - ## `Julia`
   -  ### `MH`
     - `MH_w_dnc.jl`: Julia script implementing Metropolis–Hastings algorithm along with Divide and Conquer (D&C)-recentering approach
     - `MH_wo_dnc.jl`: Julia script implementing Metropolis–Hastings algorithm along without D&C
   -  ### `normal_approx`
 - ## `R`
   -  ### `MH`
     - `MH_w_dnc.jl`: R script implementing Metropolis–Hastings algorithm along with Divide and Conquer (D&C)-recentering approach
     - `MH_wo_dnc.jl`: R script implementing Metropolis–Hastings algorithm along without D&C
   -  ### `normal_approx`
 - ## `Rcpp`
   - Note all Rcpp functions used in the various analyses can be found in the DNC package in the following GitHub repository: https://github.com/a1arakkal/DNC/tree/main
   -  ### `MH`
     - `MH_w_dnc.jl`: R script implementing Metropolis–Hastings algorithm along with Divide and Conquer (D&C)-recentering approach
     - `MH_wo_dnc.jl`: R script implementing Metropolis–Hastings algorithm along without D&C
   -  ### `normal_approx`
 - ## `sim_data.R`
   - Conatins the R-script to generate the simulated logistic regression data used in the various analyses.
