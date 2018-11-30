
"""
Abstract stiffness topology optimization problem. All subtypes must have the following fields:
    ch::ConstraintHandler
    black::BitVector
    white::BitVector
    varind::AbstractVector{Int}
    metadata::Metadata
"""
abstract type StiffnessTopOptProblem{dim, T} <: AbstractTopOptProblem end

whichdevice(p::StiffnessTopOptProblem) = whichdevice(p.ch)
whichdevice(ch::ConstraintHandler) = whichdevice(ch.dh)
whichdevice(dh::DofHandler) = whichdevice(dh.grid)
whichdevice(g::JuAFEM.Grid) = whichdevice(g.cells)

# Fallbacks
getdim(::StiffnessTopOptProblem{dim, T}) where {dim, T} = dim
floattype(::StiffnessTopOptProblem{dim, T}) where {dim, T} = T
getE(p::StiffnessTopOptProblem) = p.E
getν(p::StiffnessTopOptProblem) = p.ν
getgeomorder(p::StiffnessTopOptProblem) = nnodespercell(p) == 9 ? 2 : 1
getdensity(::StiffnessTopOptProblem{dim, T}) where {dim, T} = T(0)
getmetadata(p::StiffnessTopOptProblem) = p.metadata
getdh(p::StiffnessTopOptProblem) = p.ch.dh
getcloaddict(p::StiffnessTopOptProblem{dim, T}) where {dim, T} = Dict{String, Vector{T}}()
getpressuredict(p::StiffnessTopOptProblem{dim, T}) where {dim, T} = Dict{String, T}()
getfacesets(p::StiffnessTopOptProblem{dim, T}) where {dim, T} = Dict{String, Tuple{Int, T}}()
JuAFEM.getncells(problem::StiffnessTopOptProblem) = JuAFEM.getncells(getdh(problem).grid)

"""
Stiffness problem imported from a .inp file.
"""
struct InpStiffness{dim, N, TF, M, TI, GO, TInds<:AbstractVector{TI}, TMeta<:Metadata} <: StiffnessTopOptProblem{dim, TF}
    inp_content::InpParser.InpContent{dim, TF, N, TI}
    geom_order::Type{Val{GO}}
    ch::ConstraintHandler{DofHandler{dim, N, TF, M}, TF}
    black::BitVector
    white::BitVector
    varind::TInds
    metadata::TMeta
end

"""
Imports stiffness problem from a .inp file.
"""
function InpStiffness(filepath_with_ext::AbstractString)
    problem = InpParser.extract_inp(filepath_with_ext)
    return InpStiffness(problem)
end
function InpStiffness(problem::InpParser.InpContent)
    ch = InpParser.inp_to_juafem(problem)
    black, white = find_black_and_white(ch.dh)
    varind = find_varind(black, white)
    metadata = Metadata(ch.dh)
    geom_order = JuAFEM.getorder(ch.dh.field_interpolations[1])
    return InpStiffness(problem, Val{geom_order}, ch, black, white, varind, metadata)
end

getE(p::InpStiffness) = p.inp_content.E
getν(p::InpStiffness) = p..inp_content.ν
nnodespercell(::InpStiffness{dim, N}) where {dim, N} = N
getgeomorder(p::InpStiffness{dim, N, TF, M, TI, GO}) where {dim, N, TF, M, TI, GO} = GO
getdensity(p::InpStiffness) = p.inp_content.density
getpressuredict(p::InpStiffness) = p.inp_content.dloads
getcloaddict(p::InpStiffness) = p.inp_content.cloads
getfacesets(p::InpStiffness) = p.inp_content.facesets

"""
```
///**********************************
///*                                *
///*                                * |
///*                                * |
///********************************** v


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

`force`: force at the center right of the cantilever beam (positive is downward)

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
force = 1.0;
problem = PointLoadCantilever(nels, sizes, E, ν, force)
```
"""
struct PointLoadCantilever{dim, T, N, M, TInds<:AbstractVector{Int}, TG <: RectilinearGrid{dim, T, N, M}, TMeta<:Metadata, CH<:ConstraintHandler{<:DofHandler{dim, N, T, M}, T}} <: StiffnessTopOptProblem{dim, T}
    rect_grid::TG
    E::T
    ν::T
    ch::CH
    force::T
    force_dof::Int
    black::BitVector
    white::BitVector
    varind::TInds
    metadata::TMeta
end

function PointLoadCantilever(::Type{Val{CellType}}, nels::NTuple{dim,Int}, sizes::NTuple{dim}, E, ν, force) where {dim, CellType}
    iseven(nels[2]) && (length(nels) < 3 || iseven(nels[3])) || throw("Grid does not have an even number of elements along the y and/or z axes.")

    _T = promote_type(eltype(sizes), typeof(E), typeof(ν), typeof(force))
    if _T <: Integer
        T = Float64
    else
        T = _T
    end
    if CellType === :Linear || dim === 3
        rect_grid = RectilinearGrid(Val{:Linear}, nels, T.(sizes))
    else
        rect_grid = RectilinearGrid(Val{:Quadratic}, nels, T.(sizes))
    end

    if haskey(rect_grid.grid.facesets, "fixed_all") 
        pop!(rect_grid.grid.facesets, "fixed_all")
    end
    #addfaceset!(rect_grid.grid, "fixed_all", x -> left(rect_grid, x));
    addnodeset!(rect_grid.grid, "fixed_all", x -> left(rect_grid, x));
    
    if haskey(rect_grid.grid.nodesets, "down_force") 
        pop!(rect_grid.grid.nodesets, "down_force")
    end
    addnodeset!(rect_grid.grid, "down_force", x -> right(rect_grid, x) && middley(rect_grid, x));

    # Create displacement field u
    dh = DofHandler(rect_grid.grid)
    if CellType === :Linear || dim === 3
        push!(dh, :u, dim) # Add a displacement field
    else
        ip = Lagrange{2, RefCube, 2}()
        push!(dh, :u, dim, ip) # Add a displacement field        
    end
    close!(dh)
    
    ch = ConstraintHandler(dh)

    #dbc = Dirichlet(:u, getfaceset(rect_grid.grid, "fixed_all"), (x,t) -> zeros(T, dim), collect(1:dim))
    dbc = Dirichlet(:u, getnodeset(rect_grid.grid, "fixed_all"), (x,t) -> zeros(T, dim), collect(1:dim))
    add!(ch, dbc)
    close!(ch)
    t = T(0)
    update!(ch, t)

    metadata = Metadata(dh)
    
    fnode = Tuple(getnodeset(rect_grid.grid, "down_force"))[1]
    node_dofs = metadata.node_dofs
    force_dof = node_dofs[2, fnode]

    N = nnodespercell(rect_grid)
    M = nfacespercell(rect_grid)

    black, white = find_black_and_white(dh)
    varind = find_varind(black, white)
    
    return PointLoadCantilever(rect_grid, E, ν, ch, force, force_dof, black, white, varind, metadata)
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

`force`: force at the top left of half the MBB (positive is downward)

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
struct HalfMBB{dim, T, N, M, TInds<:AbstractVector{Int}, TG<:RectilinearGrid{dim, T, N, M}, TMeta<:Metadata, CH<:ConstraintHandler{<:DofHandler{dim, N, T, M}, T}} <: StiffnessTopOptProblem{dim, T}
    rect_grid::TG
    E::T
    ν::T
    ch::CH
    force::T
    force_dof::Int
    black::BitVector
    white::BitVector
    varind::TInds
    metadata::TMeta
end
function HalfMBB(::Type{Val{CellType}}, nels::NTuple{dim,Int}, sizes::NTuple{dim}, E, ν, force) where {dim, CellType}
    _T = promote_type(eltype(sizes), typeof(E), typeof(ν), typeof(force))
    if _T <: Integer
        T = Float64
    else
        T = _T
    end
    if CellType === :Linear || dim === 3
        rect_grid = RectilinearGrid(Val{:Linear}, nels, T.(sizes))
    else
        rect_grid = RectilinearGrid(Val{:Quadratic}, nels, T.(sizes))
    end

    if haskey(rect_grid.grid.facesets, "fixed_u1")
        pop!(rect_grid.grid.facesets, "fixed_u1")
    end
    #addfaceset!(rect_grid.grid, "fixed_u1", x -> left(rect_grid, x));
    addnodeset!(rect_grid.grid, "fixed_u1", x -> left(rect_grid, x));
    
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
    if CellType === :Linear || dim === 3
        push!(dh, :u, dim)
    else
        ip = Lagrange{2, RefCube, 2}()
        push!(dh, :u, dim, ip)
    end
    close!(dh)
    
    ch = ConstraintHandler(dh)
    #dbc1 = Dirichlet(:u, getfaceset(rect_grid.grid, "fixed_u1"), (x,t)->T[0], [1])
    dbc1 = Dirichlet(:u, getnodeset(rect_grid.grid, "fixed_u1"), (x,t)->T[0], [1])
    add!(ch, dbc1)
    dbc2 = Dirichlet(:u, getnodeset(rect_grid.grid, "fixed_u2"), (x,t)->T[0], [2])
    add!(ch, dbc2)
    close!(ch)

    t = T(0)
    update!(ch, t)

    metadata = Metadata(dh)

    fnode = Tuple(getnodeset(rect_grid.grid, "down_force"))[1]
    node_dofs = metadata.node_dofs
    force_dof = node_dofs[2, fnode]

    N = nnodespercell(rect_grid)
    M = nfacespercell(rect_grid)

    black, white = find_black_and_white(dh)
    varind = find_varind(black, white)

    return HalfMBB(rect_grid, E, ν, ch, force, force_dof, black, white, varind, metadata)
end

nnodespercell(p::Union{PointLoadCantilever, HalfMBB}) = nnodespercell(p.rect_grid)
function getcloaddict(p::Union{PointLoadCantilever{dim, T}, HalfMBB{dim, T}}) where {dim, T}
    f = T[0, -p.force, 0]
    fnode = Tuple(getnodeset(p.rect_grid.grid, "down_force"))[1]
    return Dict{Int, Vector{T}}(fnode => f)
end

"""
```
////////////
............
.          .
.          .
.          . 
.          .                    
.          ......................
.                               .
.                               . 
.                               . |
................................. v
                                force
```
"""
struct LBeam{T, N, M, TInds<:AbstractVector{Int}, TMeta<:Metadata, CH<:ConstraintHandler{<:DofHandler{2, N, T, M}, T}} <: StiffnessTopOptProblem{2, T}
    E::T
    ν::T
    ch::CH
    force::T
    force_dof::Int
    black::BitVector
    white::BitVector
    varind::TInds
    metadata::TMeta
end
function LBeam(::Type{Val{CellType}}, ::Type{T}=Float64; length = 100, height = 100, upperslab = 50, lowerslab = 50, E = 1.0, ν = 0.3, force = 1.0) where {T, CellType}
    # Create displacement field u
    grid = LGrid(Val{CellType}, T, length=length, height=height, upperslab=upperslab, 
        lowerslab=lowerslab)

    dh = DofHandler(grid)
    if CellType === :Linear
        push!(dh, :u, 2)
    else
        ip = Lagrange{2, RefCube, 2}()
        push!(dh, :u, 2, ip)
    end
    close!(dh)
    
    ch = ConstraintHandler(dh)
    dbc = Dirichlet(:u, getfaceset(grid, "top"), (x,t)->T[0, 0], [1, 2])
    add!(ch, dbc)
    close!(ch)

    t = T(0)
    update!(ch, t)

    metadata = Metadata(dh)

    fnode = Tuple(getnodeset(grid, "load"))[1]
    node_dofs = metadata.node_dofs
    force_dof = node_dofs[2, fnode]

    black, white = find_black_and_white(dh)
    varind = find_varind(black, white)

    TInds = typeof(varind)
    TMeta = typeof(metadata)
    return LBeam(E, ν, ch, force, force_dof, black, white, varind, metadata)
end

function boundingbox(grid::JuAFEM.Grid{dim}) where dim
    xmin1 = minimum(n->n.x[1], grid.nodes)
    xmax1 = maximum(n->n.x[1], grid.nodes)
    xmin2 = minimum(n->n.x[2], grid.nodes)
    xmax2 = maximum(n->n.x[2], grid.nodes)
    if dim === 2
        return ((xmin1, xmin2), (xmax1, xmax2))
    else
        xmin3 = minimum(n->n.x[3], grid.nodes)
        xmax3 = maximum(n->n.x[3], grid.nodes)
        return ((xmin1, xmin2, xmin3), (xmax1, xmax2, xmax3))
    end
end

function RectilinearTopology(b, topology = ones(getncells(getdh(b).grid)))
    bb = boundingbox(getdh(b).grid)
    go = getgeomorder(b)
    nels = Int.(round.(bb[2] .- bb[1]))
    dim = length(nels)
    if go === 1
        rectgrid = generate_grid(Quadrilateral, nels, Vec{dim}(bb[1]), Vec{dim}(bb[2]))
    elseif go === 2
        rectgrid = generate_grid(QuadraticQuadrilateral, nels, Vec{dim}(bb[1]), Vec{dim}(bb[2]))
    else
        throw("Unsupported geometry.")
    end
    new_topology = zeros(prod(nels))
    for (i, cell) in enumerate(CellIterator(getdh(b).grid))
        sub = Int.(round.((cell.coords[1]...,))) .+ (1, 1)
        ind = sub2ind(nels, sub...,)
        new_topology[ind] = topology[i]
    end
    return reshape(new_topology, nels)'
end

nnodespercell(p::LBeam{T, N}) where {T, N} = N
getdim(::LBeam) = 2
function getcloaddict(p::LBeam{T}) where {T}
    f = T[0, -p.force]
    fnode = Tuple(getnodeset(getdh(p).grid, "load"))[1]
    return Dict{Int, Vector{T}}(fnode => f)
end

"""
```
                                                               1
                                                               
                                                              OOO
                                                              ...
                                                              . .
                                                           4  . . 
                                30                            . .   
/ . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . <-
/ .                                                                 . <- 2 f 
/ .    3                                                            . <- 
/ . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . <-
                                                              ^^^
                                                              |||
                                                              1 f
```
"""
struct TieBeam{T, N, M, TInds<:AbstractVector{Int}, TMeta<:Metadata, CH<:ConstraintHandler{<:DofHandler{2, N, T, M}, T}} <: StiffnessTopOptProblem{2, T}
    E::T
    ν::T
    force::T
    ch::CH
    black::BitVector
    white::BitVector
    varind::TInds
    metadata::TMeta
end
function TieBeam(::Type{Val{CellType}}, ::Type{T}=Float64, refine = 1, force=T(1); E=T(1), ν=T(0)) where {T, CellType}
    grid = TieBeamGrid(Val{CellType}, T, refine)
    dh = DofHandler(grid)
    if CellType === :Linear
        push!(dh, :u, 2)
    else
        ip = Lagrange{2, RefCube, 2}()
        push!(dh, :u, 2, ip)
    end
    close!(dh)
        
    ch = ConstraintHandler(dh)
    dbc1 = Dirichlet(:u, getfaceset(grid, "leftfixed"), (x,t)->T[0, 0], [1, 2])
    add!(ch, dbc1)
    dbc2 = Dirichlet(:u, getfaceset(grid, "toproller"), (x,t)->T[0], [2])
    add!(ch, dbc2)
    close!(ch)

    t = T(0)
    update!(ch, t)

    metadata = Metadata(dh)

    black, white = find_black_and_white(dh)
    varind = find_varind(black, white)

    return TieBeam(E, ν, force, ch, black, white, varind, metadata)
end

getdim(::TieBeam) = 2
nnodespercell(::TieBeam{T, N}) where {T, N} = N
getpressuredict(p::TieBeam{T}) where {T} = Dict{String, T}("rightload"=>2*p.force, "bottomload"=>-p.force)
getfacesets(p::TieBeam) = getdh(p).grid.facesets
