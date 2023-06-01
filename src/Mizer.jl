module Mizer

# Export Types
export Params, Rates, Sim
# Export functions
export project, get_rates, get_rates!
# Export constants
export NS_params

using Dates, DataFrames, RCall, LoopVectorization, Tullio, LinearAlgebra, Serialization
include("params.jl")
include("sim.jl")
include("rates.jl")
include("project.jl")

datapath = 
NS_params = deserialize(pkgdir(Mizer, "data", "NS_params.jls"))

end