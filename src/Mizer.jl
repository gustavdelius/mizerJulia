module Mizer

# Export Types
export Params, Rates, Sim
# Export functions
export project, get_rates, get_rates!, get_encounter!, get_mort!
# Export constants
export NS_params, small_params

using Dates, DataFrames, RCall, LoopVectorization, Tullio, LinearAlgebra, Serialization
include("params.jl")
include("sim.jl")
include("rates.jl")
include("project.jl")

# Load data
# NS_params from mizer package
NS_params = deserialize(pkgdir(Mizer, "data", "NS_params.jls"))
# A small trait-based model with two species and 36 size classes
small_params = deserialize(pkgdir(Mizer, "data", "small_params.jls"))

end