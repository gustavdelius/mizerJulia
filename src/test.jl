using Dates, DataFrames, RCall, LoopVectorization, Tullio, LinearAlgebra
include("params.jl")
include("rates.jl")
include("project.jl")


R"""
library(mizer)
params <- NS_params
params <- setPredKernel(params, pred_kernel = getPredKernel(params))
rates_mizer <- getRates(params)
""";
@rget params;
typeof(params)

n = copy(params.initial_n);
n_pp = copy(params.initial_n_pp);
effort = copy(params.initial_effort);
r = get_rates(params, n, n_pp, effort);

import BenchmarkTools: @benchmark, @btime
@btime get_rates(params, n, n_pp, effort);
# 213.398 μs (14 allocations: 99.64 KiB) in office
@btime get_rates!(r, params, n, n_pp, effort);
# 204.844 μs (0 allocations: 0 bytes) in office

r = get_rates(params, n, n_pp, effort);
r_mizer = @rget rates_mizer;

isapprox(r.encounter, r_mizer[:encounter])
isapprox(r.one_minus_feeding_level, 1.0 .- r_mizer[:feeding_level])
isapprox(r.e_growth, r_mizer[:e_growth])
isapprox(r.pred_mort, r_mizer[:pred_mort])
isapprox(r.mort, r_mizer[:mort])
isapprox(r.rdd, r_mizer[:rdd])
isapprox(r.resource_mort, r_mizer[:resource_mort])

# Benchmark individual rate functions
@btime get_encounter!(r.encounter, r.pred_rate, r.e, params, n, n_pp);
# 70.988 μs (0 allocations: 0 bytes)

@btime get_one_minus_feeding_level!(r.one_minus_feeding_level, params, r.encounter);
# 2.646 μs (0 allocations: 0 bytes)

@btime get_e_repro_and_growth!(r.e, params, r.encounter, r.one_minus_feeding_level);
# 1.741 μs (0 allocations: 0 bytes)

@btime get_e_repro!(r.e_repro, params, r.e);
# 1.518 μs (0 allocations: 0 bytes)

@btime get_e_growth!(r.e_growth, r.e_repro, r.e);
# 1.530 μs (0 allocations: 0 bytes)

@btime get_pred_rate!(r.pred_rate, n, params.pred_rate_kernel, r.one_minus_feeding_level);
# 67.230 μs (1 allocation: 240 bytes)

@btime get_pred_mort!(r.pred_mort, params, n, n_pp, r.pred_rate);
# Calculate fishing mortality
@btime get_f_mort!(r.f_mort, params, effort)
# Calculate total mortality \mu_i(w)
@btime get_mort!(r.mort, params, r.f_mort, r.pred_mort);

## Reproduction ----
@btime get_rdi!(r.rdi, params, n, r.e_repro)
# 1.961 μs (0 allocations: 0 bytes)

@btime get_rdd!(r.rdd, r.rdi, params.R_max)

## Resource ----
# Calculate mortality on the resource spectrum
@btime get_resource_mort!(r.resource_mort, params, r.pred_rate)

@btime resource_dynamics!(n_pp, params, r, 0.1);
# 2.328 μs (0 allocations: 0 bytes)

# Test project

R"""
sim <- project(params, t_max = 100)
n_final <- finalN(sim)
n_pp_final <- finalNResource(sim)
""";
@rget n_final;
@rget n_pp_final;

effort = params.initial_effort;
n, n_pp = project(params, effort = effort);

isapprox(n, n_final)
isapprox(n_pp, n_pp_final)

import BenchmarkTools: @btime
@btime project(params, effort = $effort, t_max = 100);
# 264.589 ms (4022 allocations: 511.80 KiB) in office

using Profile
import ProfileView: @profview
@profview project(params, effort = effort, t_max = 1);
@profview project(params, effort = effort, t_max = 1000);
