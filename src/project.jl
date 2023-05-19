function project(params::Params; effort, t_max = 100, t_save = 1, dt = 0.1)
    n = copy(params.initial_n)
    n_pp = copy(params.initial_n_pp)
    # TODO: Implement saving of results in a Sim object and move allocations from
    # project_simple to here.
    num_steps = Int(t_max / dt)
    project_simple(params, n = n, n_pp = n_pp, effort = effort, 
                   num_steps = num_steps, dt = dt)
end

function project_simple(params::Params; n, n_pp, effort, num_steps = 1000, dt = 0.1)
    b = similar(n)
    r = allocate_rates(n, n_pp)
    num_sp, num_w = size(n)
    w_min_idx_array_ref = (params.w_min_idx .- 1) * num_sp + (1:num_sp)
    n_egg = view(n, w_min_idx_array_ref)
    t = 0.0
    for i in 1:num_steps
        get_rates!(r, params, n, n_pp, effort)
        resource_dynamics!(n_pp, params, r, dt)
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

function resource_dynamics!(n_pp::Vector, params::Params, r::Rates,
                            dt::Float64)
    for i in eachindex(n_pp)
        temp1 = r.resource_mort[i] + params.rr_pp[i]
        temp2 = params.rr_pp[i] * params.cc_pp[i] / temp1
        temp2 = temp2 + (n_pp[i] - temp2) * exp(-temp1 * dt)
        if isfinite(temp2)
            n_pp[i] = temp2
        end
    end
    nothing
end