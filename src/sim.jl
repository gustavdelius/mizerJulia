struct Sim
    params::params
    n::Matrix{Float64}
    n_pp::Vector{Float64}
    effort::Matrix{Float64}
end

function allocate_sim(params::Params, n, n_pp, effort)
    num_sp, num_w = size(n)
    num_t = size(effort, 2)
    n_t = zeros(num_sp, num_w, num_t)
    @views n_t[:, :, 1] = copy(n)
    n_pp_t = zeros(length(n_pp), num_t)
    @views n_pp_t[:, 1] = copy(n_pp)
    Sim(deepcopy(params), n_t, n_pp_t, effort)
end