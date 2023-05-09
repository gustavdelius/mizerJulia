"""
    mizer_rates(params, n, n_pp, n_other; t=0, effort, rates_fns, kwargs...)

Get all rates needed to project the standard mizer model.

Calls other rate functions in sequence and collects the results in a dictionary.

By default, this function returns a dictionary with the following components:
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

However, you can replace any of these rate functions with your own rate function if you wish. See `set_rate_function()` for details.

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
- Dictionary of rates.
"""
function mizer_rates(params, n, n_pp, n_other;
  t=0, effort, rates_fns, kwargs...)
  r = Dict{String,Any}()

  ## Growth ----
  # Calculate rate E_{e,i}(w) of encountered food
  r["encounter"] = juliaEncounter(params, n, n_pp, n_other, t, kwargs...)
  # Calculate feeding level f_i(w)
  r["feeding_level"] = juliaFeedingLevel(params, n, n_pp, n_other,
    encounter=r["encounter"], t=t, kwargs...)
  # Calculate the energy available for reproduction and growth
  r["e"] = juliaEReproAndGrowth(params, n, n_pp, n_other,
    encounter=r["encounter"], feeding_level=r["feeding_level"], t=t, kwargs...)
  # Calculate the energy for reproduction
  r["e_repro"] = juliaERepro(params, n, n_pp, n_other,
    e=r["e"], t=t, kwargs...)
  # Calculate the growth rate g_i(w)
  r["e_growth"] = juliaEGrowth(params, n, n_pp, n_other,
    e_repro=r["e_repro"], e=r["e"], t=t, kwargs...)

  ## Mortality ----
  # Calculate the predation rate
  r["pred_rate"] = juliaPredRate(params, n, n_pp, n_other,
    feeding_level=r["feeding_level"], t=t, kwargs...)
  # Calculate predation mortality on fish \mu_{p,i}(w)
  r["pred_mort"] = juliaPredMort(params, n, n_pp, n_other,
    pred_rate=r["pred_rate"], t=t, kwargs...)
  # Calculate fishing mortality
  r["f_mort"] = juliaFMort(params, n, n_pp, n_other,
    effort=effort, t=t,
    e_growth=r["e_growth"], pred_mort=r["pred_mort"], kwargs...)
  # Calculate total mortality \mu_i(w)
  r["mort"] = juliaMort(params, n, n_pp, n_other,
    f_mort=r["f_mort"], pred_mort=r["pred_mort"], t=t, kwargs...)

  ## Reproduction ----
  # R_di
  r["rdi"] = juliaRDI(params, n, n_pp, n_other,
    e_growth=r["e_growth"],
    mort=r["mort"],
    e_repro=r["e_repro"], t=t, kwargs...)
  # R_dd
  r["rdd"] = juliaRDD(rdi=r["rdi"], species_params=params.species_params, kwargs...)

  ## Resource ----
  # Calculate mortality on the resource spectrum
  r["resource_mort"] = juliaResourceMort(params, n, n_pp, n_other,
    pred_rate=r["pred_rate"], t=t, kwargs...)
  return r
end