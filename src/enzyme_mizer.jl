# Differentiating mizer
using Mizer, Enzyme, BenchmarkTools
Enzyme.API.runtimeActivity!(true)
params = small_params;
#params = NS_params;
n = copy(params.initial_n);
n_pp = copy(params.initial_n_pp);
effort = copy(params.initial_effort);
r = get_rates(params, n, n_pp, effort);

# Differentiate get_encounter!()
∂z_∂n = zero(n);
∂z_∂n_pp = zero(n_pp);
E = similar(n);
∂z_∂E = ones(size(E));
P = zeros(size(n, 1), length(n_pp));
Q = similar(n);
∂z_∂P = zero(P);
∂z_∂Q = zero(Q);
get_encounter!(E, P, Q, params, n, n_pp)
autodiff(Reverse, get_encounter!, Const,
         Duplicated(E, ∂z_∂E), Duplicated(P, ∂z_∂P),
         Duplicated(Q, ∂z_∂Q), Const(params),
         Duplicated(n, ∂z_∂n), Const(n_pp))
@btime Enzyme.autodiff(Reverse, get_encounter!, Const,
                       Duplicated(E, ∂z_∂E), Duplicated(P, ∂z_∂P),
                       Duplicated(Q, ∂z_∂Q), Const(params),
                       Duplicated(n, ∂z_∂n), Const(n_pp))



# Differentiate get_pred_rate!()

∂z_∂n = zero(n);
pred_rate = similar(r.pred_rate);
∂z_∂pred_rate = ones(size(pred_rate));

autodiff(Reverse, Mizer.get_pred_rate!, Const,
         Duplicated(pred_rate, ∂z_∂pred_rate), Duplicated(n, ∂z_∂n),
         Const(params.pred_rate_kernel), Const(r.one_minus_feeding_level))

@btime autodiff(Reverse, Mizer.get_pred_rate!, Const,
               Duplicated(pred_rate, ∂z_∂pred_rate), Duplicated(n, ∂z_∂n),
               Const(params.pred_rate_kernel), Const(r.one_minus_feeding_level))

# Differentiate get_growth!
function get_growth!(g, params::Params, n, n_pp, effort)
    r = get_rates(params, n, n_pp, effort)
    g .= r.e_growth
end

g = similar(n);
∂z_∂g = ones(size(g));
∂z_∂n = zero(n);

autodiff(Reverse, get_growth!, Const,
         Duplicated(g, ∂z_∂g), Const(params),
         Duplicated(n, ∂z_∂n), Const(n_pp), Const(effort))

# Differentiate get_mort!
function get_mort!(m, params::Params, n, n_pp, effort)
    r = get_rates(params, n, n_pp, effort)
    m .= r.mort
end

m = similar(n);
∂z_∂m = ones(size(m));
∂z_∂n = zero(n);

autodiff(Reverse, get_mort!, Const,
         Duplicated(m, ∂z_∂m), Const(params),
         Duplicated(n, ∂z_∂n), Const(n_pp), Const(effort))

# Differentiate get_pred_mort!
function get_resource_mort!(m, params::Params, n, n_pp, effort)
    r = get_rates(params, n, n_pp, effort)
    m .= r.resource_mort
end

m = similar(r.resource_mort);
∂z_∂m = ones(size(m));
∂z_∂n = zero(n);

autodiff(Reverse, get_resource_mort!, Const,
         Duplicated(m, ∂z_∂m), Const(params),
         Duplicated(n, ∂z_∂n), Const(n_pp), Const(effort))
