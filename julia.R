library(mizer)
library(JuliaCall)

julia <- julia_setup()

julia_eval("sqrt(2)")

julia_source("rates.jl")

params <- setPredKernel(NS_params, pred_kernel = getPredKernel(NS_params))

julia_assign("params", NS_params)
julia_eval("params.initial_effort")

julia_console()


juliaEncounter <- function(params, n, n_pp, n_other, t,  ...) {
    julia_call("encounter", params, n, n_pp, t = t)
}
params <- NS_params
params <- setRateFunction(params, "Encounter", "juliaEncounter")
waldo::compare(getEncounter(NS_params), getEncounter(params))

juliaFeedingLevel <- function(params, n, n_pp, n_other, t, encounter, ...) {
    julia_call("feeding_level", params, n, n_pp, n_other, t, encounter)
}
params <- NS_params
params <- setRateFunction(params, "FeedingLevel", "juliaFeedingLevel")
waldo::compare(getFeedingLevel(NS_params), getFeedingLevel(params))

juliaEReproAndGrowth <- function(params, n, n_pp, n_other, t, encounter,
                                 feeding_level, ...) {
    julia_call("e_repro_and_growth", params, n, n_pp, t = t, encounter = encounter,
               feeding_level = feeding_level)
}
params <- NS_params
params <- setRateFunction(params, "EReproAndGrowth", "juliaEReproAndGrowth")
waldo::compare(getEReproAndGrowth(NS_params), getEReproAndGrowth(params))

sim <- project(params, t_max = 10)

