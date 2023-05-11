library(mizer)
library(JuliaCall)

julia <- julia_setup()

julia_eval("sqrt(2)")

julia_source("rates.jl")

julia_assign("params", NS_params)
julia_eval("params.initial_effort")

julia_console()


juliaEncounter <- function(params, n, n_pp, n_other, t,  ...) {
    julia_call("get_encounter", params, n, n_pp)
}
params <- NS_params
params <- setRateFunction(params, "Encounter", "juliaEncounter")
waldo::compare(getEncounter(NS_params), getEncounter(params))

juliaFeedingLevel <- function(params, n, n_pp, n_other, t, encounter, ...) {
    julia_call("get_feeding_level", params, encounter)
}

params <- NS_params
params <- setRateFunction(params, "FeedingLevel", "juliaFeedingLevel")
waldo::compare(getFeedingLevel(NS_params), getFeedingLevel(params))

juliaRates <- function(params, n, n_pp, n_other, t = 0, effort, rates_fns, ...) {
    julia_call("julia_rates", params, n, n_pp, effort)
}
params <- NS_params
params <- setRateFunction(params, "Rates", "juliaRates")
params_slow <- setPredKernel(NS_params, getPredKernel(NS_params))

r <- getRates(params)
r_old <- getRates(NS_params)
waldo::compare(r$rdd, r_old$rdd)
waldo::compare(r$encounter, r_old$encounter)
waldo::compare(r$feeding_level, r_old$feeding_level)

library(microbenchmark)
microbenchmark(getRates(NS_params), getRates(params_slow), getRates(params))

