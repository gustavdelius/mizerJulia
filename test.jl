include("params.jl")
include("rates.jl")

R"""
library(mizer)
params <- NS_params
params <- setPredKernel(params, pred_kernel = getPredKernel(params))
rates_mizer <- getRates(params)
"""
@rget params;
typeof(params)

n = params.initial_n;
n_pp = params.initial_n_pp;
effort = params.initial_effort;

r = allocate_rates(params, n, n_pp)

r = get_rates(params, n, n_pp, effort)

import BenchmarkTools: @benchmark, @btime
@btime get_rates(params, n, n_pp, effort);
# 283.681 μs (45 allocations: 171.52 KiB)
@btime get_rates!(r, params, n, n_pp, effort);
# 269.983 μs (31 allocations: 71.88 KiB)

r = get_rates(params, n, n_pp, effort);
rm = @rget rates_mizer;

isapprox(r.encounter, rm[:encounter], rtol = 1e-6)
isapprox(r.feeding_level, rm[:feeding_level], rtol = 1e-6)
isapprox(r.e_growth, rm[:e_growth], rtol = 1e-6)
isapprox(r.pred_mort, rm[:pred_mort], rtol = 1e-6)
isapprox(r.mort, rm[:mort], rtol = 1e-6)
isapprox(r.rdd, rm[:rdd], rtol = 1e-6)
isapprox(r.resource_mort, rm[:resource_mort], rtol = 1e-6)

isapprox.(r.resource_mort, rm[:resource_mort], rtol = 1e-3)