using Dates
using DataFrames
using RCall

struct Params
    # metadata::Dict{String,Any}
    # mizer_version::Any
    # extensions::String
    # time_created::DateTime
    # time_modified::DateTime
    w::Array{Float64, 1}
    dw::Array{Float64, 1}
    w_full::Array{Float64, 1}
    dw_full::Array{Float64, 1}
    w_min_idx::Array{Int, 1}
    maturity::Array{Float64, 2}
    psi::Array{Float64, 2}
    initial_n::Array{Float64, 2}
    intake_max::Array{Float64, 2}
    search_vol::Array{Float64, 2}
    metab::Array{Float64, 2}
    pred_kernel::Array{Float64, 3}
    # ft_pred_kernel_e::Array{Float64,2}
    # ft_pred_kernel_p::Array{Float64,2}
    mu_b::Array{Float64, 2}
    rr_pp::Array{Float64, 1}
    cc_pp::Array{Float64, 1}
    # resource_dynamics::String
    resource_params::Dict{String, Any}
    # other_dynamics::Dict{String,Any}
    # other_params::Dict{String,Any}
    # other_encounter::Dict{String,Any}
    # other_mort::Dict{String,Any}
    # rates_funcs::Dict{String,String}
    # sc::Array{Float64,1}
    initial_n_pp::Array{Float64, 1}
    # initial_n_other::Dict{String,Any}
    species_params::DataFrame
    interaction::Array{Float64, 2}
    gear_params::DataFrame
    selectivity::Array{Float64, 3}
    catchability::Array{Float64, 2}
    initial_effort::Array{Float64, 1}
    # A::Array{Float64,1}
    # linecolour::Array{String,1
    # linetype::Array{String,1
    # ft_mask::Array{Float64,2}
end

import RCall.rcopy
function rcopy(::Type{Params}, robj::Ptr{S4Sxp})
    # Check that the R object is a MizerParams object
    if !convert(Bool, rcall(:is, robj, "MizerParams")[])
        throw(ArgumentError("RObject is not a MizerParams object"))
    end
    if isna(robj[:pred_kernel], 1)
        pred_kernel = rcall(:getPredKernel, robj)
    else
        pred_kernel = robj[:pred_kernel]
    end
    Params(
           # rcopy(Dict{Any,Any}, robj[:metadata]),
           # rcopy(Any, robj[:mizer_version]),
           # rcopy(String, robj[:extensions]),
           # rcopy(DateTime, robj[:time_created]),
           # rcopy(DateTime, robj[:time_modified]),
           rcopy(Array{Float64, 1}, robj[:w]),
           rcopy(Array{Float64, 1}, robj[:dw]),
           rcopy(Array{Float64, 1}, robj[:w_full]),
           rcopy(Array{Float64, 1}, robj[:dw_full]),
           rcopy(Array{Int, 1}, robj[:w_min_idx]),
           rcopy(Array{Float64, 2}, robj[:maturity]),
           rcopy(Array{Float64, 2}, robj[:psi]),
           rcopy(Array{Float64, 2}, robj[:initial_n]),
           rcopy(Array{Float64, 2}, robj[:intake_max]),
           rcopy(Array{Float64, 2}, robj[:search_vol]),
           rcopy(Array{Float64, 2}, robj[:metab]),
           rcopy(Array{Float64, 3}, pred_kernel),
           # rcopy(Array{Float64,2}, robj[:ft_pred_kernel_e]),
           # rcopy(Array{Float64,2}, robj[:ft_pred_kernel_p]),
           rcopy(Array{Float64, 2}, robj[:mu_b]),
           rcopy(Array{Float64, 1}, robj[:rr_pp]),
           rcopy(Array{Float64, 1}, robj[:cc_pp]),
           # rcopy(String, robj[:resource_dynamics]),
           rcopy(Dict{String, Any}, robj[:resource_params]),
           # rcopy(Dict{Any,Any}, robj[:other_dynamics]),
           # rcopy(Dict{Any,Any}, robj[:other_params]),
           # rcopy(Dict{Any,Any}, robj[:other_encounter]),
           # rcopy(Dict{Any,Any}, robj[:other_mort]),
           # rcopy(Dict{Any,Any}, robj[:rates_funcs]),
           # rcopy(Array{Float64,1}, robj[:sc]),
           rcopy(Array{Float64, 1}, robj[:initial_n_pp]),
           # rcopy(Dict{String,Any}, robj[:initial_n_other]),
           rcopy(DataFrame, robj[:species_params]),
           rcopy(Array{Float64, 2}, robj[:interaction]),
           rcopy(DataFrame, robj[:gear_params]),
           rcopy(Array{Float64, 3}, robj[:selectivity]),
           rcopy(Array{Float64, 2}, robj[:catchability]),
           rcopy(Array{Float64, 1}, robj[:initial_effort])
           # rcopy(Array{Float64,1}, robj[:A]),
           # rcopy(String, robj[:linecolour]),
           # rcopy(String, robj[:linetype]),
           # rcopy(Array{Float64,3}, robj[:ft_mask])
           )
end

import RCall: RClass, rcopytype
rcopytype(::Type{RClass{:MizerParams}}, s::Ptr{S4Sxp}) = Params

using RCall
R"""
library(mizer)
params <- NS_params
"""

@rget params;
typeof(params)

"""
    juliaRates(params, n, n_pp, n_other; t=0, effort, rates_fns, kwargs...)

Get all rates needed to project the standard mizer model.

Calls other rate functions in sequence and collects the results in a dictionary.

# Arguments
- `params`: A `MizerParams` object.
- `n`: A matrix of species abundances (species x size).
- `n_pp`: A vector of the resource abundance by size.
- `n_other`: A list of abundances for other dynamical components of the ecosystem.
- `t`: The time for which to do the calculation (Not used by standard mizer rate functions but useful for extensions with time-dependent parameters). Default is `0`.
- `effort`: The effort for each fishing gear.
- `rates_fns`: Named dictionary of the functions to call to calculate the rates. Note that this dictionary holds the functions themselves, not their names.
- `kwargs`: Unused.

# Returns
- Dictionary of rates with the following components:
- `encounter` from `mizer_encounter()`
- `feeding_level` from `mizer_feeding_level()`
- `e` from `mizer_e_repro_and_growth()`
- `e_repro` from `mizer_e_repro()`
- `e_growth` from `mizer_e_growth()`
- `pred_rate` from `mizer_pred_rate()`
- `pred_mort` from `mizer_pred_mort()`
- `f_mort` from `mizer_f_mort()`
- `mort` from `mizer_mort()`
- `rdi` from `mizer_rdi()`
- `rdd` from `beverton_holt_rdd()`
- `resource_mort` from `mizer_resource_mort()`
"""
function julia_rates(params::Params, n, n_pp; effort, t = 0)
    r = Dict{String, Any}()

    ## Growth ----
    # Calculate rate E_{e,i}(w) of encountered food
    r["encounter"] = encounter(params, n, n_pp, t)
    # Calculate feeding level f_i(w)
    r["feeding_level"] = feeding_level(params, n, n_pp,
                                       encounter = r["encounter"], t = t)
    # Calculate the energy available for reproduction and growth
    r["e"] = e_repro_and_growth(params, n, n_pp, encounter = r["encounter"],
                                feeding_level = r["feeding_level"], t = t)
    # # Calculate the energy for reproduction
    # r["e_repro"] = juliaERepro(params, n, n_pp, e = r["e"], t = t, kwargs...)
    # # Calculate the growth rate g_i(w)
    # r["e_growth"] = juliaEGrowth(params, n, n_pp, e_repro = r["e_repro"],
    #                              e = r["e"], t = t, kwargs...)

    # ## Mortality ----
    # # Calculate the predation rate
    # r["pred_rate"] = juliaPredRate(params, n, n_pp, n_other,
    #                                feeding_level = r["feeding_level"], t = t, kwargs...)
    # # Calculate predation mortality on fish \mu_{p,i}(w)
    # r["pred_mort"] = juliaPredMort(params, n, n_pp, n_other,
    #                                pred_rate = r["pred_rate"], t = t, kwargs...)
    # # Calculate fishing mortality
    # r["f_mort"] = juliaFMort(params, n, n_pp, n_other, effort = effort, t = t,
    #                          e_growth = r["e_growth"], pred_mort = r["pred_mort"],
    #                          kwargs...)
    # # Calculate total mortality \mu_i(w)
    # r["mort"] = juliaMort(params, n, n_pp, n_other, f_mort = r["f_mort"],
    #                       pred_mort = r["pred_mort"], t = t, kwargs...)

    # ## Reproduction ----
    # # R_di
    # r["rdi"] = juliaRDI(params, n, n_pp, n_other, e_growth = r["e_growth"],
    #                     mort = r["mort"],
    #                     e_repro = r["e_repro"], t = t, kwargs...)
    # # R_dd
    # r["rdd"] = juliaRDD(rdi = r["rdi"], species_params = params.species_params, kwargs...)

    # ## Resource ----
    # # Calculate mortality on the resource spectrum
    # r["resource_mort"] = juliaResourceMort(params, n, n_pp, n_other,
    #                                        pred_rate = r["pred_rate"], t = t, kwargs...)
    return r
end

function encounter(params::Params, n, n_pp; t = 0)
    E = zeros(size(n))
    num_sp, num_w = size(n)
    for i in 1:num_sp
        for wi in 1:num_w
            for j in 1:num_sp
                for wj in 1:num_w
                    E[i, wi] += params.interaction[i, j] * n[j, wj] *
                                params.pred_kernel[i, wi, wj] * params.w[wj] *
                                params.dw[wj]
                end
            end
            for wj in eachindex(n_pp)
                E[i, wi] += params.species_params.interaction_resource[i] * n_pp[wj] *
                            params.pred_kernel[i, wi, wj] * params.w_full[wj] *
                            params.dw_full[wj]
            end
            E[i, wi] *= params.search_vol[i, wi]
        end
    end

    return E
end

function feeding_level(params::Params, n, n_pp; t, encounter)
    encounter ./ (encounter .+ params.intake_max)
end

function e_repro_and_growth(params::Params, n, n_pp; t, encounter, feeding_level)
    (1.0 .- feeding_level) .* encounter .* params.species_params.alpha .- params.metab
end

n = params.initial_n;
n_pp = params.initial_n_pp;
effort = params.initial_effort;
r = julia_rates(params, n, n_pp, effort = effort, t = 0)
