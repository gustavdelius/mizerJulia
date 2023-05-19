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
# 277.203 μs (35 allocations: 140.28 KiB)
@btime get_rates!(r, params, n, n_pp, effort);
# 260.167 μs (20 allocations: 31.14 KiB)
# 379.200 μs (16 allocations: 432 bytes) on slow notebook

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
# 89.700 μs (3 allocations: 336 bytes)

@btime get_feeding_level!(r.one_minus_feeding_level, params, r.encounter);
# 2.000 μs (2 allocations: 224 bytes)

@btime get_e_repro_and_growth!(r.e, params, r.encounter, r.one_minus_feeding_level);
# 4.157 μs (8 allocations: 560 bytes)

@btime get_e_repro!(r.e_repro, params, r.e);
# 2.011 μs (2 allocations: 224 bytes)

@btime get_e_growth!(r.e_growth, r.e_repro, r.e);
# 2.067 μs (3 allocations: 336 bytes)

@btime get_pred_rate!(r.pred_rate, n, params.pred_rate_kernel, r.one_minus_feeding_level);
# 87.700 μs (3 allocations: 432 bytes)

@btime get_pred_mort!(r.pred_mort, params, n, n_pp, r.pred_rate);
# Calculate fishing mortality
@btime get_f_mort!(r.f_mort, params, effort)
# Calculate total mortality \mu_i(w)
@btime get_mort!(r.mort, params, r.f_mort, r.pred_mort);

## Reproduction ----
# R_di
@btime get_rdi!(r.rdi, params, n, r.e_repro)
# R_dd,
@btime get_rdd!(r.rdd, r.rdi, params.species_params.R_max)

## Resource ----
# Calculate mortality on the resource spectrum
@btime get_resource_mort!(r.resource_mort, params, r.pred_rate)

n_pp = copy(params.initial_n_pp)
temp1 = similar(n_pp)
temp2 = similar(n_pp)

@btime resource_dynamics!(n_pp, temp1, temp2, params, r, 0.1);

# Test project

R"""
sim <- project(params, t_max = 100)
n_final <- finalN(sim)
n_pp_final <- finalNResource(sim)
""";
@rget n_final;
@rget n_pp_final;

include("project.jl")
effort = params.initial_effort;
n, n_pp = project(params, effort = effort);

isapprox(n, n_final)
isapprox(n_pp, n_pp_final)

import BenchmarkTools: @btime
@btime project(params, effort = $effort);

using Profile
using PProf

# Collect an allocation profile
Profile.Allocs.clear()
Profile.Allocs.@profile project(params, effort = effort, num_steps = 10000);

# Export pprof allocation profile and open interactive profiling web interface.
PProf.Allocs.pprof()

# Collect a profile
Profile.clear()
@profile project(params, effort = effort);

# Export pprof profile and open interactive profiling web interface.
pprof()

import ProfileView: @profview
@profview project(params, effort = effort);
using Profile
@profile project(params, effort = effort);
import Juno
Juno.profiler()