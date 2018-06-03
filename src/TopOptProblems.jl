module TopOptProblems

using JuAFEM
using StaticArrays

abstract type AbstractTopOptProblem end

include("grids.jl")
include("metadata.jl")
include("stiffness_problems_types.jl")
include("inp_to_juafem.jl")
include("assemble_f.jl")
include("matrices_and_vectors.jl")
include("utils.jl")

export PointLoadCantilever, HalfMBB, InpStiffness, StiffnessTopOptProblem, AbstractTopOptProblem, GlobalFEAInfo, ElementFEAInfo, YoungsModulus

end # module
