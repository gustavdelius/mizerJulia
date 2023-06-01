struct Sim
    params::Params
    n::Array{Float64, 3}
    n_pp::Matrix{Float64}
    # TODO: Allow effort to be a function of time
    # effort::Matrix{Float64}
end

function allocate_sim(params::Params, n, n_pp, num_t)
    num_sp, num_w = size(n)
    n_t = zeros(num_sp, num_w, num_t + 1)
    @views n_t[:, :, 1] = copy(n)
    n_pp_t = zeros(length(n_pp), num_t + 1)
    @views n_pp_t[:, 1] = copy(n_pp)
    Sim(deepcopy(params), n_t, n_pp_t)
end