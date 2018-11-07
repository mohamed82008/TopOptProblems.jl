module TopOptProblems

using JuAFEM
using StaticArrays
#using Makie
#using GeometryTypes

import JuAFEM: assemble!

abstract type AbstractTopOptProblem end

include("grids.jl")
include("metadata.jl")
include("problem_types.jl")
include("matrices_and_vectors.jl")
include("penalties.jl")
include("assemble.jl")
include("utils.jl")
#include("makie.jl")

export PointLoadCantilever, HalfMBB, LBeam, TieBeam, InpStiffness, StiffnessTopOptProblem, AbstractTopOptProblem, GlobalFEAInfo, ElementFEAInfo, YoungsModulus, assemble, assemble_f!

end # module
