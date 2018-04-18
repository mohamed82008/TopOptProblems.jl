# TopOptProblems

[![Build Status](https://travis-ci.org/mohamed82008/TopOptProblems.jl.svg?branch=master)](https://travis-ci.org/mohamed82008/TopOptProblems.jl)

[![Coverage Status](https://coveralls.io/repos/mohamed82008/TopOptProblems.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/mohamed82008/TopOptProblems.jl?branch=master)

[![codecov.io](http://codecov.io/github/mohamed82008/TopOptProblems.jl/coverage.svg?branch=master)](http://codecov.io/github/mohamed82008/TopOptProblems.jl?branch=master)


This package allows the definition of topology optimization problem contexts outside the topology optimization code. The problem context includes the mesh, boundary conditions, loads as well as any elements to be fixed as black or white.

```julia
julia> using TopOptProblems

help?> HalfMBB
search: HalfMBB

   |
   |
   v
  O*********************************
  O*                               *
  O*                               *
  O*                               *
  O*********************************
                                   O

  struct HalfMBB{dim, T, N, M} <: StiffnessTopOptProblem{dim, T}
      rect_grid::RectilinearGrid{dim, T, N, M}
      E::T
      ν::T
      ch::ConstraintHandler{DofHandler{dim, N, T, M}, T}
      force::T
      force_dof::Int
      metadata::Metadata
  end

  dim: dimension of the problem

  T: number type for computations and coordinates

  N: number of nodes in a cell of the grid

  M: number of faces in a cell of the grid

  rect_grid: a RectilinearGrid struct

  E: Young's modulus

  ν: Poisson's ration

  ch: a JuAFEM.ConstraintHandler struct

  force: force at the top left of half the MBB (negative is downward)

  force_dof: dof number at which the force is applied

  metadata:: Metadata having various cell-node-dof relationships

  API:

          HalfMBB(nels::NTuple{dim,Int}, sizes::NTuple{dim}, E, ν, force) where {dim, T}

  Example:


  nels = (60,20);
  sizes = (1.0,1.0);
  E = 1.0;
  ν = 0.3;
  force = -1.0;
  problem = HalfMBB(nels, sizes, E, ν, force)

julia>  nels = (60,20);

julia>   sizes = (1.0,1.0);

julia>   E = 1.0;

julia>   ν = 0.3;

julia>   force = -1.0;

julia>   problem = HalfMBB(nels, sizes, E, ν, force)
TopOptProblems.HalfMBB{2,Float64,4,4}(TopOptProblems.RectilinearGrid{2,Float64,4,4}(JuAFEM.Grid{2,4,Float64,4} with 1200 Quadrilateral cells and 1281 nodes, (60, 20), (1.0, 1.0), ([0.0, 0.0], [60.0, 20.0]), Bool[false, false, false, false, true, true, true, true, true, false  …  false, false, false, false, false, false, false, false, false, false], Bool[false, false, false, false, true, true, true, true, true, false  …  false, false, false, false, false, false, false, false, false, false], Bool[false, false, false, false, true, true, true, true, true, false  …  false, false, false, false, false, false, false, false, false, false]), 1.0, 0.3, ConstraintHandler:
  Face BCs:
    Field: u, Components: [1]
  Node BCs:
    Field: u, Components: [2], -1.0, 2444, TopOptProblems.Metadata([1 3 … 2435 2437; 2 4 … 2436 2438; … ; 7 5 … 2557 2559; 8 6 … 2558 2560], Tuple{Int64,Int64}[(1, 1) (1, 2) … (1200, 5) (1200, 6); (0, 0) (0, 0) … (0, 0) (0, 0); (0, 0) (0, 0) … (0, 0) (0, 0); (0, 0) (0, 0) … (0, 0) (0, 0)], Tuple{Int64,Int64}[(1, 1), (1, 2), (2, 2), (3, 2), (4, 2), (5, 2), (6, 2), (7, 2), (8, 2), (9, 2)  …  (1191, 3), (1192, 3), (1193, 3), (1194, 3), (1195, 3), (1196, 3), (1197, 3), (1198, 3), (1199, 3), (1200, 3)], [1 3 … 2559 2561; 2 4 … 2560 2562]))

```
