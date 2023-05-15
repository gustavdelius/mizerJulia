include("params.jl")
include("rates.jl")

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

r = get_rates(params, n, n_pp, effort);
r_mizer = @rget rates_mizer;

isapprox(r.encounter, r_mizer[:encounter])
isapprox(r.feeding_level, r_mizer[:feeding_level])
isapprox(r.e_growth, r_mizer[:e_growth])
isapprox(r.pred_mort, r_mizer[:pred_mort])
isapprox(r.mort, r_mizer[:mort])
isapprox(r.rdd, r_mizer[:rdd])
isapprox(r.resource_mort, r_mizer[:resource_mort])

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
Profile.Allocs.@profile project(params, effort = effort);

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