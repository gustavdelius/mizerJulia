using RCall, Mizer

R"""
library(mizer)
params <- newTraitParams(no_sp = 2, no_w = 10) |> steady()
""";
@rget params;

using Serialization
serialize("data/small_params.jls", params)