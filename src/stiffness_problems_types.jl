
abstract type StiffnessTopOptProblem{dim, T} <: AbstractTopOptProblem end

"""
```
///**********************************
///*                                *
///*                                * |
///*                                * |
///++++++++++++++++++++++++++++++++++ v


struct PointLoadCantilever{dim, T, N, M} <: StiffnessTopOptProblem{dim, T}
    rect_grid::RectilinearGrid{dim, T, N, M}
    E::T
    ν::T
    ch::ConstraintHandler{DofHandler{dim, N, T, M}, T}
    force::T
    force_dof::Int
    metadata::Metadata
end
```

`dim`: dimension of the problem

`T`: number type for computations and coordinates

`N`: number of nodes in a cell of the grid

`M`: number of faces in a cell of the grid


`rect_grid`: a RectilinearGrid struct

`E`: Young's modulus

`ν`: Poisson's ration

`ch`: a JuAFEM.ConstraintHandler struct

`force`: force at the center right of the cantilever beam (negative is downward)

`force_dof`: dof number at which the force is applied

`metadata`:: Metadata having various cell-node-dof relationships


API:
```
        PointLoadCantilever(nels::NTuple{dim,Int}, sizes::NTuple{dim}, E, ν, force) where {dim}
```

Example:
```

nels = (60,20);
sizes = (1.0,1.0);
E = 1.0;
ν = 0.3;
force = -1.0;
problem = PointLoadCantilever(nels, sizes, E, ν, force)
```
"""
struct PointLoadCantilever{dim, T, N, M} <: StiffnessTopOptProblem{dim, T}
    rect_grid::RectilinearGrid{dim, T, N, M}
    E::T
    ν::T
    ch::ConstraintHandler{DofHandler{dim, N, T, M}, T}
    force::T
    force_dof::Int
    metadata::Metadata
end

function PointLoadCantilever(nels::NTuple{dim,Int}, sizes::NTuple{dim}, E, ν, force) where {dim}
    iseven(nels[2]) && (length(nels) < 3 || iseven(nels[3])) || throw("Grid does not have an even number of elements along the y and/or z axes.")

    _T = promote_type(eltype(sizes), typeof(E), typeof(ν), typeof(force))
    if _T <: Integer
        T = Float64
    else
        T = _T
    end
    rect_grid = RectilinearGrid(nels, T.(sizes))

    if haskey(rect_grid.grid.facesets, "fixed_all") 
        pop!(rect_grid.grid.facesets, "fixed_all")
    end
    addfaceset!(rect_grid.grid, "fixed_all", x -> left(rect_grid, x));
    
    if haskey(rect_grid.grid.nodesets, "down_force") 
        pop!(rect_grid.grid.nodesets, "down_force")
    end
    addnodeset!(rect_grid.grid, "down_force", x -> right(rect_grid, x) && middley(rect_grid, x));

    # Create displacement field u
    dh = DofHandler(rect_grid.grid)
    push!(dh, :u, dim) # Add a displacement field
    close!(dh)
    
    ch = ConstraintHandler(dh)

    dbc = Dirichlet(:u, getfaceset(rect_grid.grid, "fixed_all"), (x,t) -> zeros(T, dim), collect(1:dim))
    add!(ch, dbc)
    close!(ch)
    t = T(0)
    update!(ch, t)

    metadata = Metadata(dh)
    
    fnode = Tuple(getnodeset(rect_grid.grid, "down_force"))[1]
    force_dof = metadata.node_dofs[2, fnode]

    N = nnodespercell(rect_grid)
    M = nfacespercell(rect_grid)

    return PointLoadCantilever{dim, T, N, M}(rect_grid, E, ν, ch, force, force_dof, metadata)
end

"""
```
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
```

`dim`: dimension of the problem

`T`: number type for computations and coordinates

`N`: number of nodes in a cell of the grid

`M`: number of faces in a cell of the grid

`rect_grid`: a RectilinearGrid struct

`E`: Young's modulus

`ν`: Poisson's ration

`ch`: a JuAFEM.ConstraintHandler struct

`force`: force at the top left of half the MBB (negative is downward)

`force_dof`: dof number at which the force is applied

`metadata`:: Metadata having various cell-node-dof relationships

API:
```
        HalfMBB(nels::NTuple{dim,Int}, sizes::NTuple{dim}, E, ν, force) where {dim}
```

Example:
```

nels = (60,20);
sizes = (1.0,1.0);
E = 1.0;
ν = 0.3;
force = -1.0;
problem = HalfMBB(nels, sizes, E, ν, force)
```
"""
struct HalfMBB{dim, T, N, M} <: StiffnessTopOptProblem{dim, T}
    rect_grid::RectilinearGrid{dim, T, N, M}
    E::T
    ν::T
    ch::ConstraintHandler{DofHandler{dim, N, T, M}, T}
    force::T
    force_dof::Int
    metadata::Metadata
end
function HalfMBB(nels::NTuple{dim,Int}, sizes::NTuple{dim}, E, ν, force) where {dim}
    _T = promote_type(eltype(sizes), typeof(E), typeof(ν), typeof(force))
    if _T <: Integer
        T = Float64
    else
        T = _T
    end
    rect_grid = RectilinearGrid(nels, T.(sizes))

    if haskey(rect_grid.grid.facesets, "fixed_u1")
        pop!(rect_grid.grid.facesets, "fixed_u1")
    end        
    addfaceset!(rect_grid.grid, "fixed_u1", x -> left(rect_grid, x));
    
    if haskey(rect_grid.grid.nodesets, "fixed_u2")
        pop!(rect_grid.grid.nodesets, "fixed_u2")
    end
    addnodeset!(rect_grid.grid, "fixed_u2", x -> bottom(rect_grid, x) && right(rect_grid, x));

    if haskey(rect_grid.grid.nodesets, "down_force")
        pop!(rect_grid.grid.nodesets, "down_force")
    end
    addnodeset!(rect_grid.grid, "down_force", x -> top(rect_grid, x) && left(rect_grid, x));

    # Create displacement field u
    dh = DofHandler(rect_grid.grid)
    push!(dh, :u, dim)
    close!(dh)
    
    ch = ConstraintHandler(dh)
    dbc1 = Dirichlet(:u, getfaceset(rect_grid.grid, "fixed_u1"), (x,t)->T[0], [1])
    add!(ch, dbc1)
    dbc2 = Dirichlet(:u, getnodeset(rect_grid.grid, "fixed_u2"), (x,t)->T[0], [2])
    add!(ch, dbc2)
    close!(ch)

    t = T(0)
    update!(ch, t)

    metadata = Metadata(dh)

    fnode = Tuple(getnodeset(rect_grid.grid, "down_force"))[1]
    force_dof = metadata.node_dofs[2, fnode]

    N = nnodespercell(rect_grid)
    M = nfacespercell(rect_grid)

    return HalfMBB{dim, T, N, M}(rect_grid, E, ν, ch, force, force_dof, metadata)
end