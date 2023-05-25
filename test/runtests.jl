using Mizer, RCall
using Test

@testset "get_rates" begin
    # We get the NS_params from the R mizer package
    R"""
    library(mizer)
    params <- NS_params
    # We set the predation kernel like this so that mizer does not
    # use fast Fourier transforms
    params <- setPredKernel(params, pred_kernel = getPredKernel(params))
    rates_mizer <- getRates(params)
    """
    @rget params
    r_mizer = @rget rates_mizer
    @test typeof(params) == Mizer.Params
    @test typeof(r_mizer) == Mizer.Rates

    n = params.initial_n
    n_pp = params.initial_n_pp
    effort = params.initial_effort
    r = get_rates(params, n, n_pp, effort)

    @test isapprox(r.encounter, r_mizer[:encounter])
    @test isapprox(r.one_minus_feeding_level, 1.0 .- r_mizer[:feeding_level])
    @test isapprox(r.e_growth, r_mizer[:e_growth])
    @test isapprox(r.pred_mort, r_mizer[:pred_mort])
    @test isapprox(r.mort, r_mizer[:mort])
    @test isapprox(r.rdd, r_mizer[:rdd])
    @test isapprox(r.resource_mort, r_mizer[:resource_mort])
end

@testset "project" begin
    R"""
    library(mizer)
    params <- NS_params
    params <- setPredKernel(params, pred_kernel = getPredKernel(params))
    sim <- project(params, t_max = 100)
    n_final <- finalN(sim)
    n_pp_final <- finalNResource(sim)
    """
    @rget params
    @rget n_final
    @rget n_pp_final

    effort = params.initial_effort
    n, n_pp = project(params, effort = effort)

    @test isapprox(n, n_final)
    @test isapprox(n_pp, n_pp_final)
end
