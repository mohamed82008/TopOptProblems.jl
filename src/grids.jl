abstract type AbstractGrid{dim, T} end

"""
```
struct RectilinearGrid{dim, T, N, M} <: AbstractGrid{dim, T}
    grid::JuAFEM.Grid{dim, N, T, M}
    nels::NTuple{dim, Int}
    sizes::NTuple{dim, T}
    corners::NTuple{2, Vec{dim,T}}
    white_cells::BitVector
    black_cells::BitVector
    constant_cells::BitVector
end
```

`dim`: dimension of the problem

`T`: number type for computations and coordinates

`N`: number of nodes in a cell of the grid

`M`: number of faces in a cell of the grid


`grid`: a JuAFEM.Grid struct

`nels`: number of elements in every dimension

`sizes`: dimensions of each rectilinear cell

`corners`: 2 corner points of the rectilinear grid

`white_cells`: cells fixed to be void during optimization

`black_cells`: cells fixed to have material during optimization

`constant_cells`: cells fixed to be either void or have material during optimization


API:

`function RectilinearGrid(nels::NTuple{dim,Int}, sizes::NTuple{dim,T}) where {dim,T}`


Example:

`rectgrid = RectilinearGrid((60,20), (1.0,1.0))`
"""
struct RectilinearGrid{dim, T, N, M} <: AbstractGrid{dim, T}
    grid::JuAFEM.Grid{dim, N, T, M}
    nels::NTuple{dim, Int}
    sizes::NTuple{dim, T}
    corners::NTuple{2, Vec{dim, T}}
    white_cells::BitVector
    black_cells::BitVector
    constant_cells::BitVector
end
nnodespercell(::RectilinearGrid{dim,T,N,M}) where {dim, T, N, M} = N
nfacespercell(::RectilinearGrid{dim,T,N,M}) where {dim, T, N, M} = M

left(rectgrid::RectilinearGrid, x) = x[1] ≈ rectgrid.corners[1][1]
right(rectgrid::RectilinearGrid, x) = x[1] ≈ rectgrid.corners[2][1]
bottom(rectgrid::RectilinearGrid, x) = x[2] ≈ rectgrid.corners[1][2]
top(rectgrid::RectilinearGrid, x) = x[2] ≈ rectgrid.corners[2][2]
back(rectgrid::RectilinearGrid, x) = x[3] ≈ rectgrid.corners[1][3]
front(rectgrid::RectilinearGrid, x) = x[3] ≈ rectgrid.corners[2][3]
middlex(rectgrid::RectilinearGrid, x) = x[1] ≈ (rectgrid.corners[1][1] + rectgrid.corners[2][1]) / 2
middley(rectgrid::RectilinearGrid, x) = x[2] ≈ (rectgrid.corners[1][2] + rectgrid.corners[2][2]) / 2
middlez(rectgrid::RectilinearGrid, x) = x[3] ≈ (rectgrid.corners[1][3] + rectgrid.corners[2][3]) / 2

function RectilinearGrid(nels::NTuple{dim,Int}, sizes::NTuple{dim,T}) where {dim,T}
    geoshape = dim === 2 ? Quadrilateral : Hexahedron
    corner1 = Vec{dim}(fill(T(0), dim))
    corner2 = Vec{dim}((nels .* sizes))
    grid = generate_grid(geoshape, nels, corner1, corner2);

    N = JuAFEM.nnodes(geoshape)
    M = JuAFEM.nfaces(geoshape)
    ncells = prod(nels)
    return RectilinearGrid{dim, T, N, M}(grid, nels, sizes, (corner1, corner2),BitVector(ncells), BitVector(ncells), BitVector(ncells))    
end
