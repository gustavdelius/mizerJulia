struct Rates
    encounter::Matrix{Float64}
    feeding_level::Matrix{Float64}
    e::Matrix{Float64}
    e_repro::Matrix{Float64}
    e_growth::Matrix{Float64}
    pred_rate::Matrix{Float64}
    pred_mort::Matrix{Float64}
    f_mort::Matrix{Float64}
    mort::Matrix{Float64}
    rdi::Vector{Float64}
    rdd::Vector{Float64}
    resource_mort::Vector{Float64}
end

function allocate_rates(n, n_pp)
    Rates(similar(n), # encounter
          similar(n), # feeding_level
          similar(n), # e
          similar(n), # e_repro
          similar(n), # e_growth
          zeros(size(n, 1), length(n_pp)), # pred_rate
          similar(n), # pred_mort
          similar(n), # f_mort
          similar(n), # mort
          zeros(size(n, 1)), # rdi
          zeros(size(n, 1)), # rdd
          similar(n_pp)) # resource_mort
end

function get_rates(params::Params, n, n_pp, effort)
    r = allocate_rates(n, n_pp)
    get_rates!(r, params::Params, n, n_pp, effort)
end

function get_rates!(r::Rates, params::Params, n, n_pp, effort)

    ## Growth ----
    # Calculate rate E_{e,i}(w) of encountered food
    # r.pred_rate and r.e only passed so that the function has preallocated
    # arrays to hold temporary results
    get_encounter!(r.encounter, r.pred_rate, r.e, params, n, n_pp)
    # Calculate feeding level f_i(w)
    get_feeding_level!(r.feeding_level, params, r.encounter)
    # Calculate the energy available for reproduction and growth
    get_e_repro_and_growth!(r.e, params, r.encounter, r.feeding_level)
    # Calculate the energy for reproduction
    get_e_repro!(r.e_repro, params, r.e)
    # Calculate the growth rate g_i(w)
    get_e_growth!(r.e_growth, r.e_repro, r.e)

    ## Mortality ----
    # Calculate the predation rate
    get_pred_rate!(r.pred_rate, n, params.pred_rate_kernel, r.feeding_level)
    # Calculate predation mortality on fish \mu_{p,i}(w)
    get_pred_mort!(r.pred_mort, params, n, n_pp, r.pred_rate)
    # Calculate fishing mortality
    get_f_mort!(r.f_mort, params, effort)
    # Calculate total mortality \mu_i(w)
    get_mort!(r.mort, params, r.f_mort, r.pred_mort)

    ## Reproduction ----
    # R_di
    get_rdi!(r.rdi, params, n, r.e_repro)
    # R_dd,
    get_rdd!(r.rdd, r.rdi, params.species_params.R_max)

    ## Resource ----
    # Calculate mortality on the resource spectrum
    get_resource_mort!(r.resource_mort, params, r.pred_rate)
    return r
end

function padded_add!(B, A)
    start = size(B, 2) - size(A, 2) + 1
    B[:, start:size(B, 2)] .+= A
    return B
end

function prey!(P, Q, interaction, n, interaction_resource, n_pp)
    mul!(P, interaction_resource, n_pp')
    mul!(Q, interaction, n)
    P .= padded_add!(P, Q)
end

function get_encounter!(E, P, Q, params, n, n_pp)
    P = prey!(P, Q, params.interaction, n, params.species_params.interaction_resource, n_pp)
    K = params.encounter_kernel
    @tullio E[i, w] = P[i, wp] * K[i, w, wp]
end

function get_feeding_level!(feeding_level, params::Params, encounter)
    feeding_level .= encounter ./ (encounter .+ params.intake_max)
end

function get_e_repro_and_growth!(e, params::Params, encounter, feeding_level)
    e .= (1.0 .- feeding_level) .* encounter .* params.species_params.alpha .- params.metab
end

function get_e_repro!(e_repro, params::Params, e)
    e_repro .= max.(e .* params.psi, 0)
end

function get_e_growth!(e_growth, e_repro, e)
    e_growth .= max.(e .- e_repro, 0)
end

function get_pred_rate!(pred_rate, n, pred_rate_kernel, feeding_level)
    @tullio pred_rate[i, wj] = pred_rate_kernel[i, wi, wj] * (1.0 .- feeding_level)[i, wi] *
                               n[i, wi]
end

function get_pred_mort!(pred_mort, params::Params, n, n_pp, pred_rate)
    num_sp, num_w = size(n)
    num_w_pp = length(n_pp)
    idx = (num_w_pp + 1 - num_w):num_w_pp
    mul!(pred_mort, transpose(params.interaction), @view pred_rate[:, idx])
end

function get_f_mort!(f_mort, params::Params, effort)
    @tullio f_mort[i, w] = params.selectivity[g, i, w] * params.catchability[g, i] *
                           effort[g]
end

function get_mort!(mort, params::Params, f_mort, pred_mort)
    mort .= f_mort .+ pred_mort .+ params.mu_b
end

function get_rdi!(rdi, params::Params, n, e_repro)
    rdi .= 0.5 .* (((e_repro .* n) * params.dw) .* params.species_params.erepro) ./
           params.w[params.w_min_idx]
end

function get_rdd!(rdd, rdi, R_max)
    rdd .= rdi ./ (1.0 .+ rdi ./ R_max)
end

function get_resource_mort!(resource_mort, params::Params, pred_rate)
    # @tullio resource_mort[w] = params.species_params.interaction_resource[i] * pred_rate[i, w]
    resource_mort .= (params.species_params.interaction_resource' * pred_rate)'
end
