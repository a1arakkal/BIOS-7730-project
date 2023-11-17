using DataFrames, GLM, CSV, BenchmarkTools

function run_glm()
    mod_dat = CSV.read("/Users/atlan/Adv_computing/AdvComp_final_proj/data/sim_data.csv",
                   DataFrame, header = 1)
                   
    fm = @formula(y ~ X1 + X2 + X3 + X4)

    logit = glm(fm, mod_dat, Binomial(), LogitLink())
    temp = GLM.coeftable(logit)

    out = DataFrame(variable = temp.rownms,
                    Estimate = temp.cols[1],
                    StdError = temp.cols[2],
                    z_val = temp.cols[3],
                    lowerCI = temp.cols[5],
                    upperCI = temp.cols[6])
    return out
end

run_glm()
run_time = @elapsed res = run_glm()
run_time = round(run_time, digits = 4)
CSV.write("/Users/atlan/Adv_computing/AdvComp_final_proj/data/glm_$run_time.csv", res)