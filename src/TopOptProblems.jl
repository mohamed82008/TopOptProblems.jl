module TopOptProblems

using JuAFEM
using StaticArrays

abstract type AbstractTopOptProblem end

include("grids.jl")
include("metadata.jl")
include("stiffness_problems_types.jl")
include("matrices_and_vectors.jl")

export RectilinearGrid, PointLoadCantilever, HalfMBB

end # module
