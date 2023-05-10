library(mizer)
library(JuliaCall)

julia <- julia_setup()

julia_eval("sqrt(2)")

julia_source("rates.jl")

params <- setPredKernel(NS_params, pred_kernel = getPredKernel(NS_params))

julia_assign("params", NS_params)
julia_eval("params.initial_effort")

julia_console()

juliaFeedingLevel <- function(params, n, n_pp, n_other, t, encounter, ...) {
    julia_call("feeding_level", params, n, n_pp, n_other, t, encounter)
}

params <- NS_params
params <- setRateFunction(params, "FeedingLevel", "juliaFeedingLevel")

sim <- project(params, t_max = 10)

waldo::compare(getFeedingLevel(NS_params), getFeedingLevel(params))
