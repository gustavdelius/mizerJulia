using Enzyme, Random

# Example from autodiff documentation
a = 4.2
bb = [2.2, 3.3];
∂f_∂bb = zero(bb);
c = 55;
d = 9;

f(a, bb, c, d) = a * √(bb[1]^2 + bb[2]^2) + c^2 * d^2
∂f_∂a, _, _, ∂f_∂d = autodiff(Reverse, f, Active(a), Duplicated(bb, ∂f_∂bb), c,
                              Active(d))[1]
∂f_∂bb

# Example from "Implementing pullbacks" https://enzyme.mit.edu/julia/stable/pullbacks/

function mymul!(R, A, B)
    @assert axes(A,2) == axes(B,1)
    @inbounds @simd for i in eachindex(R)
        R[i] = 0
    end
    @inbounds for j in axes(B, 2), i in axes(A, 1)
        @inbounds @simd for k in axes(A,2)
            R[i,j] += A[i,k] * B[k,j]
        end
    end
    nothing
end

Random.seed!(1234)
A = rand(5, 3)
B = rand(3, 7)

R = zeros(size(A,1), size(B,2))
∂z_∂R = rand(size(R)...)  # Some gradient/tangent passed to us
∂z_∂R0 = copyto!(similar(∂z_∂R), ∂z_∂R)  # exact copy for comparison

∂z_∂A = zero(A)
∂z_∂B = zero(B)

Enzyme.autodiff(Reverse, mymul!, Const, Duplicated(R, ∂z_∂R), Duplicated(A, ∂z_∂A), Duplicated(B, ∂z_∂B))

∂z_∂A
∂z_∂B

R ≈ A * B &&
    ∂z_∂A ≈ ∂z_∂R0 * B' &&  # equivalent to Zygote.pullback(*, A, B)[2](∂z_∂R)[1]
    ∂z_∂B ≈ A' * ∂z_∂R0       # equivalent to Zygote.pullback(*, A, B)[2](∂z_∂R)[2]

using BenchmarkTools
@btime Enzyme.autodiff(Reverse, mymul!, Const, Duplicated(R, ∂z_∂R), Duplicated(A, ∂z_∂A), Duplicated(B, ∂z_∂B))


# Enzyme for adjoint tutorial: Stommel three-box ocean model
# from https://enzyme.mit.edu/julia/stable/generated/box/

struct ModelParameters

    # handy to have constants
    day::Float64
    year::Float64

    # Information related to the boxes
    boxlength::Vector{Float64}      ## Vector with north-south size of each box  [cm]
    boxdepth::Vector{Float64}       ## "          " the depth of each box  [cm]
    boxwidth::Float64               ## "          " the width of each box  [cm]
    boxarea::Vector{Float64}        ## "          " the area of each box   [cm^2]
    boxvol::Vector{Float64}         ## "          " the volume of each box   [cm^3]

    delta::Float64                  ## Constant ratio depth(box1) / (depth(box1) + depth(box3))

    # Parameters that appear in the box model equations
    u0::Float64
    alpha::Float64
    beta::Float64
    gamma::Float64

    # Coefficient for the Robert filter smoother
    rf_coeff::Float64

    # Freshwater forcing
    FW::Vector{Float64}

    # Restoring atmospheric temperatures and salinities
    Tstar::Vector{Float64}
    Sstar::Vector{Float64}
end

function setup()
    blength = [5000.0e5; 1000.0e5; 5000.0e5]
    bdepth = [1.0e5; 5.0e5; 4.0e5]

    delta = bdepth[1] / (bdepth[1] + bdepth[3])

    bwidth = 4000.0 * 1e5  ## box width, centimeters

    # box areas
    barea = [blength[1] * bwidth;
             blength[2] * bwidth;
             blength[3] * bwidth]

    # box volumes
    bvolume = [barea[1] * bdepth[1];
               barea[2] * bdepth[2];
               barea[3] * bdepth[3]]

    # parameters that are used to ensure units are in CGS (cent-gram-sec)

    day = 3600.0 * 24.0
    year = day * 365.0
    Sv = 1e12                       ## one Sverdrup (a unit of ocean transport), 1e6 meters^3/second

    # parameters that appear in box model equations
    u0 = 16.0 * Sv / 0.0004
    alpha = 1668e-7
    beta = 0.7811e-3

    gamma = 1 / (300 * day)

    # robert filter coefficient for the smoother part of the timestep
    robert_filter_coeff = 0.25

    # freshwater forcing
    FW = [(100 / year) * 35.0 * barea[1]; -(100 / year) * 35.0 * barea[1]]

    # restoring atmospheric temperatures
    Tstar = [22.0; 0.0]
    Sstar = [36.0; 34.0]

    structure_with_parameters = ModelParameters(day,
                                                year,
                                                blength,
                                                bdepth,
                                                bwidth,
                                                barea,
                                                bvolume,
                                                delta,
                                                u0,
                                                alpha,
                                                beta,
                                                gamma,
                                                robert_filter_coeff,
                                                FW,
                                                Tstar,
                                                Sstar)

    return structure_with_parameters
end

# function to compute transport
#       Input: rho - the density vector
#       Output: U - transport value

function compute_transport(rho, params)
    U = params.u0 * (rho[2] - (params.delta * rho[1] + (1 - params.delta) * rho[3]))
    return U
end

# function to compute density
#       Input: state = [T1; T2; T3; S1; S2; S3]
#       Output: rho

function compute_density(state, params)
    rho = -params.alpha * state[1:3] + params.beta * state[4:6]
    return rho
end

# lastly, a function that takes one step forward
#       Input: state_now = [T1(t), T2(t), ..., S3(t)]
#              state_old = [T1(t-dt), ..., S3(t-dt)]
#              u = transport(t)
#              dt = time step
#       Output: state_new = [T1(t+dt), ..., S3(t+dt)]

function compute_update(state_now, state_old, u, params, dt)
    dstate_now_dt = zeros(6)
    state_new = zeros(6)

    # first computing the time derivatives of the various temperatures and salinities
    if u > 0
        dstate_now_dt[1] = u * (state_now[3] - state_now[1]) / params.boxvol[1] +
                           params.gamma * (params.Tstar[1] - state_now[1])
        dstate_now_dt[2] = u * (state_now[1] - state_now[2]) / params.boxvol[2] +
                           params.gamma * (params.Tstar[2] - state_now[2])
        dstate_now_dt[3] = u * (state_now[2] - state_now[3]) / params.boxvol[3]

        dstate_now_dt[4] = u * (state_now[6] - state_now[4]) / params.boxvol[1] +
                           params.FW[1] / params.boxvol[1]
        dstate_now_dt[5] = u * (state_now[4] - state_now[5]) / params.boxvol[2] +
                           params.FW[2] / params.boxvol[2]
        dstate_now_dt[6] = u * (state_now[5] - state_now[6]) / params.boxvol[3]

    elseif u <= 0
        dstate_now_dt[1] = u * (state_now[2] - state_now[1]) / params.boxvol[1] +
                           params.gamma * (params.Tstar[1] - state_now[1])
        dstate_now_dt[2] = u * (state_now[3] - state_now[2]) / params.boxvol[2] +
                           params.gamma * (params.Tstar[2] - state_now[2])
        dstate_now_dt[3] = u * (state_now[1] - state_now[3]) / params.boxvol[3]

        dstate_now_dt[4] = u * (state_now[5] - state_now[4]) / params.boxvol[1] +
                           params.FW[1] / params.boxvol[1]
        dstate_now_dt[5] = u * (state_now[6] - state_now[5]) / params.boxvol[2] +
                           params.FW[2] / params.boxvol[2]
        dstate_now_dt[6] = u * (state_now[4] - state_now[6]) / params.boxvol[3]
    end

    # update fldnew using a version of Euler's method
    state_new .= state_old + 2.0 * dt * dstate_now_dt

    return state_new
end

function integrate(state_now, state_old, dt, M, parameters)

    # Because of the adjoint problem we're setting up, we need to store both the states before
    # and after the Robert filter smoother has been applied
    states_before = [state_old]
    states_after = [state_old]

    for t in 1:M
        rho = compute_density(state_now, parameters)
        u = compute_transport(rho, parameters)
        state_new = compute_update(state_now, state_old, u, parameters, dt)

        # Applying the Robert filter smoother (needed for stability)
        state_new_smoothed = state_now +
                             parameters.rf_coeff * (state_new - 2.0 * state_now + state_old)

        push!(states_after, state_new_smoothed)
        push!(states_before, state_new)

        # cycle the "now, new, old" states
        state_old = state_new_smoothed
        state_now = state_new
    end

    return states_after, states_before
end

function one_step_forward(state_now, state_old, out_now, out_old, parameters, dt)
    state_new_smoothed = zeros(6)
    rho = compute_density(state_now, parameters)                             ## compute density
    u = compute_transport(rho, parameters)                                   ## compute transport
    state_new = compute_update(state_now, state_old, u, parameters, dt)      ## compute new state values

    # Robert filter smoother
    state_new_smoothed[:] = state_now +
                            parameters.rf_coeff * (state_new - 2.0 * state_now + state_old)

    out_old[:] = state_new_smoothed
    out_now[:] = state_new

    return nothing
end

parameters = setup()

Tbar = [20.0; 1.0; 1.0]         ## initial temperatures
Sbar = [35.5; 34.5; 34.5]       ## initial salinities

# Running the model one step forward
states_after_smoother, states_before_smoother = integrate(copy([Tbar; Sbar]),
                                                          copy([Tbar; Sbar]),
                                                          10 * parameters.day,
                                                          1,
                                                          parameters)

# Run Enzyme one time on `one_step_forward``
dstate_now = zeros(6)
dstate_old = zeros(6)
out_now = zeros(6);
dout_now = ones(6);
out_old = zeros(6);
dout_old = ones(6);

autodiff(Reverse,
         one_step_forward,
         Duplicated([Tbar; Sbar], dstate_now),
         Duplicated([Tbar; Sbar], dstate_old),
         Duplicated(out_now, dout_now),
         Duplicated(out_old, dout_old),
         parameters,
         Const(10 * parameters.day))

@show out_now, out_old

@show dstate_now

# Differentiating mizer
using Mizer
params = NS_params;
n = params.initial_n;
∂z_∂n = zero(n);
n_pp = params.initial_n_pp;
∂z_∂n_pp = zero(n_pp);
E = similar(n);
∂z_∂E = ones(size(E));
P = zeros(size(n, 1), length(n_pp));
Q = similar(n);
get_encounter!(E, P, Q, params, n, n_pp)

Enzyme.autodiff(Reverse, get_encounter!, Const, 
                Duplicated(E, ∂z_∂E), Const(P), Const(Q), Const(params),
                Duplicated(n, ∂z_∂n), Const(n_pp))


Mizer.padded_add!(P, Q)

∂z_∂P = ones(size(P));
∂z_∂Q = zero(Q);
autodiff(Reverse, Mizer.padded_add!, Const, Duplicated(P, ∂z_∂P), Duplicated(Q, ∂z_∂Q))
∂z_∂Q

autodiff(Reverse, Mizer.prey!, Const, Duplicated(P, ∂z_∂P), Const(Q), 
         Const(params.interaction), Duplicated(n, ∂z_∂n), 
         Const(params.interaction_resource), Const(n_pp))
∂z_∂n

using LinearAlgebra
∂z_∂P = ones(size(P));
∂z_∂n_pp = zero(n_pp);
autodiff(Reverse, mul!, Const, Duplicated(P, ∂z_∂P), Const(params.interaction_resource), 
         Duplicated(n_pp, ∂z_∂n_pp))
∂z_∂n_pp

using FiniteDiff



using ForwardDiff

function mymul_flat(AB, m, n, p)
    A = reshape(AB[1:(m * n)], (m, n))
    B = reshape(AB[(m * n + 1):end], (n, p))
    R = zeros(m, p)
    mymul!(R, A, B)
    return vec(R)
end

m, n, p = size(A, 1), size(A, 2), size(B, 2)
AB = [vec(A); vec(B)]
#∂s_∂AB = ForwardDiff.jacobian(AB -> mymul_flat(AB, m, n, p), AB)

∂s_∂A = reshape(∂s_∂AB[1:(m * n)], (m, n))
∂s_∂B = reshape(∂s_∂AB[(m * n + 1):end], (n, p))

