module TopOptProblems

using JuAFEM
using StaticArrays
using InpParser

abstract type AbstractTopOptProblem end

include("utils.jl")
include("grids.jl")
include("metadata.jl")
include("stiffness_problems_types.jl")
include("inp_to_juafem.jl")
include("matrices_and_vectors.jl")

export RectilinearGrid, PointLoadCantilever, HalfMBB

end # module
