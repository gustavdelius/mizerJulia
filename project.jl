function project(params::Params; effort, num_steps = 1000, dt = 0.1)
    n = copy(params.initial_n)
    n_pp = copy(params.initial_n_pp)
    b = similar(n)
    temp1 = similar(n_pp)
    temp2 = similar(n_pp)
    r = allocate_rates(n, n_pp)
    num_sp, num_w = size(n)
    w_min_idx_array_ref = (params.w_min_idx .- 1) * num_sp + (1:num_sp)
    n_egg = view(n, w_min_idx_array_ref)
    t = 0.0
    for i in 1:num_steps
        get_rates!(r, params, n, n_pp, effort)
        resource_dynamics!(n_pp, temp1, temp2, params, r, dt)
        @tullio b[i, w] = 1.0 + (r.e_growth[i, w] / params.dw[w] + r.mort[i, w]) * dt
        b_egg = view(b, w_min_idx_array_ref)
        dw_egg = view(params.dw, params.w_min_idx)
        n_egg .= (n_egg .+ r.rdd .* dt ./ dw_egg) ./ b_egg
        for i in 1:num_sp
            for j in (params.w_min_idx[i] + 1):num_w
                n[i, j] = (n[i, j] + r.e_growth[i, j-1] * dt / params.dw[j] * n[i, j-1]) / b[i, j]
            end
        end
        t += dt
    end
    return n, n_pp
end

function resource_dynamics!(n_pp::Vector, temp1, temp2, params::Params, r::Rates, dt::Float64)
    temp1 .= r.resource_mort .+ params.rr_pp
    temp2 .= params.rr_pp .* params.cc_pp ./ temp1
    temp2 .= temp2 + (n_pp .- temp2) .* exp.(-temp1 .* dt)
    n_pp[isfinite.(temp2)] .= temp2[isfinite.(temp2)]
end