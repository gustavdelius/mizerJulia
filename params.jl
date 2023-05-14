using Dates
using DataFrames
using RCall
using LoopVectorization
using Tullio

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
    pred_rate_kernel::Array{Float64, 3}
    encounter_kernel::Array{Float64, 3}
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

function encounter_kernel(pred_kernel, w_full, dw_full, search_vol)
    @tullio K[i, w, wp] := pred_kernel[i, w, wp] * w_full[wp] * dw_full[wp] *
                           search_vol[i, w]
end

function pred_rate_kernel(pred_kernel, dw, search_vol)
    @tullio K[i, w, wp] := pred_kernel[i, w, wp] * dw[w] * search_vol[i, w]
end

import RCall.rcopy
function rcopy(::Type{Params}, robj::Ptr{S4Sxp})
    # Check that the R object is a MizerParams object
    if !convert(Bool, rcall(:is, robj, "MizerParams")[])
        throw(ArgumentError("RObject is not a MizerParams object"))
    end
    # Calculate pred_kernel if not yet available
    if isna(robj[:pred_kernel], 1)
        pred_kernel = rcall(:getPredKernel, robj)
    else
        pred_kernel = robj[:pred_kernel]
    end
    # Copy some arrays that will be used below to calculate the kernels
    # needed for pred_rate and encounter rate.
    dw = rcopy(Array{Float64, 1}, robj[:dw])
    w_full = rcopy(Array{Float64, 1}, robj[:w_full])
    dw_full = rcopy(Array{Float64, 1}, robj[:dw_full])
    search_vol = rcopy(Array{Float64, 2}, robj[:search_vol])
    pred_kernel = rcopy(Array{Float64, 3}, pred_kernel)
    Params(
           # rcopy(Dict{Any,Any}, robj[:metadata]),
           # rcopy(Any, robj[:mizer_version]),
           # rcopy(String, robj[:extensions]),
           # rcopy(DateTime, robj[:time_created]),
           # rcopy(DateTime, robj[:time_modified]),
           rcopy(Array{Float64, 1}, robj[:w]),
           dw,
           w_full,
           dw_full,
           rcopy(Array{Int, 1}, robj[:w_min_idx]),
           rcopy(Array{Float64, 2}, robj[:maturity]),
           rcopy(Array{Float64, 2}, robj[:psi]),
           rcopy(Array{Float64, 2}, robj[:initial_n]),
           rcopy(Array{Float64, 2}, robj[:intake_max]),
           search_vol,
           rcopy(Array{Float64, 2}, robj[:metab]),
           # We do not store the pred_kernel but instead the kernels that are
           # actually needed to calculate the pred_rate and the encounter_rate.
           pred_rate_kernel(pred_kernel, dw, search_vol),
           encounter_kernel(pred_kernel, w_full, dw_full, search_vol),
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
