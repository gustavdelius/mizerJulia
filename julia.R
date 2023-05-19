library(mizer)
library(JuliaCall)

julia <- julia_setup()

julia_eval("sqrt(2)")

julia_source("params.jl")
julia_source("rates.jl")

juliaRates <- function(params, n, n_pp, n_other, t = 0, effort, rates_fns, ...) {
    julia_call("get_rates", params, n, n_pp, effort)
}
params <- NS_params
params <- setRateFunction(params, "Rates", "juliaRates")
params_slow <- setPredKernel(NS_params, getPredKernel(NS_params))

install.packages("microbenchmark")
library(microbenchmark)
microbenchmark(getRates(NS_params), getRates(params_slow))
rates_fns <- lapply(NS_params@rates_funcs, get)
microbenchmark(mizerRates(NS_params, n = NS_params@initial_n, 
  n_pp = NS_params@initial_n_pp, effort = NS_params@initial_effort, 
  rates_fns = rates_fns), times = 1000)

# This does not work yet because get_rates() returns a struct, not a list.
r <- getRates(params)
r_old <- getRates(NS_params)
waldo::compare(r$rdd, r_old$rdd)
waldo::compare(r$encounter, r_old$encounter)
waldo::compare(r$feeding_level, r_old$feeding_level)

microbenchmark(getRates(NS_params), getRates(params_slow), getRates(params))

microbenchmark(project(NS_params))

microbenchmark(project_simple(NS_params, steps = 1000), times = 20)
