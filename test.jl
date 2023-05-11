include("rates.jl")

R"""
library(mizer)
params <- NS_params
"""
@rget params;
typeof(params)

n = params.initial_n;
n_pp = params.initial_n_pp;
effort = params.initial_effort;
r = julia_rates(params, n, n_pp, effort)

import BenchmarkTools: @benchmark
@benchmark julia_rates(params, n, n_pp, effort)

@profview julia_rates(params, n, n_pp, effort)
