const RectilinearPointLoad{dim, T, N, M} = Union{PointLoadCantilever{dim, T, N, M}, HalfMBB{dim, T, N, M}, LBeam{T, N, M}}

struct ElementFEAInfo{dim, T, TK<:AbstractMatrix{T}, Tf<:AbstractVector{T}, TKes<:AbstractVector{TK}, Tfes<:AbstractVector{Tf}, Tcload<:AbstractVector{T}, refshape, TCV<:CellValues{dim, T, refshape}, dimless1, TFV<:FaceValues{dimless1, T, refshape}}
    Kes::TKes
    fes::Tfes
    fixedload::Tcload
    cellvolumes::Tcload
    cellvalues::TCV
    facevalues::TFV
    metadata::Metadata
end
#=
function ElementFEAInfo(sp::RectilinearPointLoad{dim, T}, quad_order=2, ::Type{Val{mat_type}}=Val{:Static}) where {dim, T, mat_type}
    Kes, cellvalues, facevalues = make_Kes(sp, quad_order, Val{mat_type})
    fes = [zeros(T, ndofs_per_cell(sp.ch.dh)) for i in 1:getncells(sp.ch.dh.grid)]
    fixedload = make_f(sp)
    cellvolumes = get_cell_volumes(sp, cellvalues)
    return ElementFEAInfo(Kes, fes, fixedload, cellvolumes, cellvalues, facevalues, sp.metadata)
end
=#

function ElementFEAInfo(sp, quad_order=2, ::Type{Val{mat_type}}=Val{:Static}) where {mat_type} 
    Kes, weights, dloads, cellvalues, facevalues = make_Kes_and_fes(sp, quad_order, Val{mat_type})
    fixedload = full(make_cload(sp))
    assemble_f!(fixedload, sp, dloads)
    cellvolumes = get_cell_volumes(sp, cellvalues)
    ElementFEAInfo(Kes, weights, fixedload, cellvolumes, cellvalues, facevalues, sp.metadata)
end

function ElementFEAInfo(Kes::AbstractVector{TK}, fes::AbstractVector{Tf},fixedload::AbstractVector{T}, cellvolumes::AbstractVector{T}, cellvalues::TCV, facevalues::TFV, metadata::Metadata) where {T, TK<:AbstractMatrix{T}, Tf<:AbstractVector{T}, dim, dimless1, refshape, TCV<:CellValues{dim, T, refshape}, TFV<:FaceValues{dimless1, T, refshape}}
    ElementFEAInfo{dim, T, TK, Tf, typeof(Kes), typeof(fes), typeof(fixedload), refshape, TCV, dimless1, TFV}(Kes, fes, fixedload, cellvolumes, cellvalues, facevalues, metadata)
end

function get_cell_volumes(sp::StiffnessTopOptProblem{dim, T}, cellvalues) where {dim, T}
    dh = sp.ch.dh
    cellvolumes = zeros(T, getncells(dh.grid))
    for (i, cell) in enumerate(CellIterator(dh))
        reinit!(cellvalues, cell)
        cellvolumes[i] = sum(JuAFEM.getdetJdV(cellvalues, q_point) for q_point in 1:JuAFEM.getnquadpoints(cellvalues))
    end
    return cellvolumes
end

mutable struct GlobalFEAInfo{T, TK<:AbstractMatrix{T}, Tf<:AbstractVector{T}}
    K::TK
    f::Tf
end

GlobalFEAInfo(::Type{T}) where T = GlobalFEAInfo{T}()

GlobalFEAInfo() = GlobalFEAInfo{Float64}()

GlobalFEAInfo{T}() where T = GlobalFEAInfo{T, SparseMatrixCSC{T, Int}, Vector{T}}(sparse(zeros(T, 0, 0)), zeros(T, 0))

GlobalFEAInfo(sp::StiffnessTopOptProblem) = GlobalFEAInfo(make_empty_K(sp), make_empty_f(sp))

GlobalFEAInfo(K::AbstractMatrix{T}, f::AbstractVector{T}) where T = GlobalFEAInfo{T, typeof(K), typeof(f)}(K, f)

make_empty_K(sp::StiffnessTopOptProblem) = Symmetric(create_sparsity_pattern(sp.ch.dh))

make_empty_f(sp::StiffnessTopOptProblem{dim, T}) where {dim, T} = zeros(T, ndofs(sp.ch.dh))

#=
function make_Kes(sp::RectilinearPointLoad, quad_order=2)
    make_Kes(sp, quad_order, Val{:Static})
end

function make_Kes(sp::RectilinearPointLoad, ::Type{Val{mat_type}}) where mat_type
    make_Kes(sp, 2, Val{mat_type})
end

function make_Kes(sp::RectilinearPointLoad{dim, T}, quad_order, ::Type{Val{mat_type}}) where {dim, T, mat_type}
    E = sp.E
    ν = sp.ν

    λ = E*ν / ((1 + ν) * (1 - 2*ν))
    μ = E / (2*(1 + ν))
    δ(i,j) = i == j ? T(1) : T(0)
    g(i,j,k,l) = λ*δ(i,j)*δ(k,l) + μ*(δ(i,k)*δ(j,l) + δ(i,l)*δ(j,k))
    C = SymmetricTensor{4, dim}(g)

    dh = sp.ch.dh

    # Shape functions and quadrature rule
    interpolation_space = Lagrange{dim, RefCube, 1}()
    quadrature_rule = QuadratureRule{dim, RefCube}(quad_order)
    cellvalues = CellScalarValues(quadrature_rule, interpolation_space)
    facevalues = FaceScalarValues(QuadratureRule{dim-1, RefCube}(quad_order), interpolation_space)

    # Calculate element stiffness matrices
    n_basefuncs = getnbasefunctions(cellvalues)
    
    return _make_Kes(sp, Val{mat_type}, Val{n_basefuncs}, Val{dim*n_basefuncs}, C, interpolation_space, quadrature_rule, cellvalues, facevalues), cellvalues, facevalues
end

function _make_Kes(sp::RectilinearPointLoad{dim, T}, ::Type{Val{mat_type}}, ::Type{Val{n_basefuncs}}, ::Type{Val{ndofs_per_cell}}, C, interpolation_space, quadrature_rule, cellvalues, facevalues) where {dim, T, mat_type, n_basefuncs, ndofs_per_cell}
    dh = sp.ch.dh

    # Calculate element stiffness matrices
    Kesize = ndofs_per_cell
    nel = getncells(dh.grid)

    if mat_type === :Static || mat_type === :SMatrix
        if !(T === BigFloat)
            Kes = Symmetric{T, SMatrix{Kesize, Kesize, T, Kesize^2}}[]
        else
            Kes = Symmetric{T, SizedMatrix{Kesize, Kesize, T, Kesize^2}}[]
        end
        sizehint!(Kes, nel)

        Ke_e = zeros(T, dim, dim)
        Ke_0 = Matrix{T}(Kesize, Kesize)
        celliteratortype = CellIterator{typeof(dh).parameters...}
        _celliterator::celliteratortype = CellIterator(dh)
        for (k, cell) in enumerate(_celliterator)
            Ke_0 .= 0
            reinit!(cellvalues, cell)
            for q_point in 1:getnquadpoints(cellvalues)
                dΩ = getdetJdV(cellvalues, q_point)
                for b in 1:n_basefuncs
                    ∇ϕb = shape_gradient(cellvalues, q_point, b)
                    for d2 in 1:dim
                        for a in 1:n_basefuncs
                            ∇ϕa = shape_gradient(cellvalues, q_point, a)
                            Ke_e .= dotdot(∇ϕa, C, ∇ϕb) * dΩ
                            for d1 in 1:dim
                                #if dim*(b-1) + d2 >= dim*(a-1) + d1
                                Ke_0[dim*(a-1) + d1, dim*(b-1) + d2] += Ke_e[d1,d2]
                                #end
                            end
                        end
                    end
                end
            end
            if !(T === BigFloat)
                push!(Kes, Symmetric(SMatrix{Kesize, Kesize, T, Kesize*Kesize}(Ke_0)))
            else
                push!(Kes, Symmetric(SizedMatrix{Kesize, Kesize, T}(Ke_0)))
            end
        end
    else
        Kes = let Kesize=Kesize, nel=nel
            [Symmetric(zeros(T, Kesize, Kesize), :U) for i = 1:nel]
        end
    
        Ke_e = zeros(T, dim, dim)
        
        celliteratortype = CellIterator{typeof(dh).parameters...}
        _celliterator = CellIterator(dh)
        for (k, cell) in enumerate(_celliterator)
            reinit!(cellvalues, cell)
            for q_point in 1:getnquadpoints(cellvalues)
                dΩ = getdetJdV(cellvalues, q_point)
                for b in 1:n_basefuncs
                    ∇ϕb = shape_gradient(cellvalues, q_point, b)
                    for d2 in 1:dim
                        for a in 1:n_basefuncs
                            ∇ϕa = shape_gradient(cellvalues, q_point, a)
                            Ke_e .= dotdot(∇ϕa, C, ∇ϕb) * dΩ
                            for d1 in 1:dim
                                #if dim*(b-1) + d2 >= dim*(a-1) + d1
                                Kes[k].data[dim*(a-1) + d1, dim*(b-1) + d2] += Ke_e[d1,d2]
                                #end
                            end
                        end
                    end
                end
            end
        end
    end
    return Kes
end

function make_f(sp::RectilinearPointLoad)
    dof_cells = sp.metadata.dof_cells
    fdof = sp.force_dof
    force = sp.force
    f = sparsevec([fdof], [-force], ndofs(sp.ch.dh))
    return f
end
=#

function make_Kes_and_fes(problem, quad_order=2)
    make_Kes_and_fes(problem, quad_order, Val{:Static})
end

function make_Kes_and_fes(problem, ::Type{Val{mat_type}}) where mat_type
    make_Kes_and_fes(problem, 2, Val{mat_type})
end

function make_Kes_and_fes(problem, quad_order, ::Type{Val{mat_type}}) where {mat_type}
    T = floattype(problem)
    dim = getdim(problem)
    geom_order = getgeomorder(problem)
    dh = getdh(problem)
    E = getE(problem)
    ν = getν(problem)
    ρ = getdensity(problem)

    refshape = JuAFEM.getrefshape(dh.field_interpolations[1])

    λ = E*ν / ((1 + ν) * (1 - 2*ν))
    μ = E / (2*(1 + ν))
    δ(i,j) = i == j ? T(1) : T(0)
    g(i,j,k,l) = λ*δ(i,j)*δ(k,l) + μ*(δ(i,k)*δ(j,l) + δ(i,l)*δ(j,k))
    C = SymmetricTensor{4, dim}(g)

    # Shape functions and quadrature rule
    interpolation_space = Lagrange{dim, refshape, geom_order}()
    quadrature_rule = QuadratureRule{dim, refshape}(quad_order)
    cellvalues = CellScalarValues(quadrature_rule, interpolation_space)
    facevalues = FaceScalarValues(QuadratureRule{dim-1, refshape}(quad_order), interpolation_space)

    # Calculate element stiffness matrices
    n_basefuncs = getnbasefunctions(cellvalues)
    
    Kes, weights = _make_Kes_and_weights(dh, Val{mat_type}, Val{n_basefuncs}, Val{dim*n_basefuncs}, C, ρ, quadrature_rule, cellvalues)
    dloads = _make_dloads(weights, problem, facevalues)

    return Kes, weights, dloads, cellvalues, facevalues
end

const g = [0., 9.81, 0.] # N/kg or m/s^2

function _make_Kes_and_weights(dh::DofHandler{dim, N, T}, ::Type{Val{mat_type}}, ::Type{Val{n_basefuncs}}, ::Type{Val{ndofs_per_cell}}, C, ρ, quadrature_rule, cellvalues) where {dim, N, T, mat_type, n_basefuncs, ndofs_per_cell}
    # Calculate element stiffness matrices
    Kesize = ndofs_per_cell
    nel = getncells(dh.grid)
    body_force = ρ .* g # Force per unit volume
    
    if !(T === BigFloat)
        if mat_type === :Static || mat_type === :SMatrix
            MatrixType = SMatrix{Kesize, Kesize, T, Kesize^2}
            VectorType = MVector{Kesize, T}
        elseif mat_type === :MMatrix
            MatrixType = MMatrix{Kesize, Kesize, T, Kesize^2}
            VectorType = MVector{Kesize, T}
        else
            MatrixType = Matrix{T}
            VectorType = Vector{T}
        end
    else
        if mat_type === :Static || mat_type === :SMatrix  || mat_type === :MMatrix
            MatrixType = SizedMatrix{Kesize, Kesize, T, Kesize^2}
            VectorType = SizedVector{Kesize, T}
        else
            MatrixType = Matrix{T}
            VectorType = Vector{T}
        end
    end

    if MatrixType <: StaticArray
        Kes = Symmetric{T, MatrixType}[]
        sizehint!(Kes, nel)
        weights = [zeros(VectorType) for i in 1:nel]
        
        Ke_e = zeros(T, dim, dim)
        fe = zeros(T, Kesize)
        Ke_0 = Matrix{T}(Kesize, Kesize)
        celliteratortype = CellIterator{typeof(dh).parameters...}
        _celliterator::celliteratortype = CellIterator(dh)
        for (k, cell) in enumerate(_celliterator)
            Ke_0 .= 0
            try
                reinit!(cellvalues, cell)
            catch
                @show cell.coords, cell.nodes, cell.current_cellid[]
            end
            fe = weights[k]
            for q_point in 1:getnquadpoints(cellvalues)
                dΩ = getdetJdV(cellvalues, q_point)
                for b in 1:n_basefuncs
                    ∇ϕb = shape_gradient(cellvalues, q_point, b)
                    ϕb = shape_value(cellvalues, q_point, b)
                    for d2 in 1:dim
                        fe[(b-1)*dim + d2] += ϕb * body_force[d2] * dΩ
                        for a in 1:n_basefuncs
                            ∇ϕa = shape_gradient(cellvalues, q_point, a)
                            Ke_e .= dotdot(∇ϕa, C, ∇ϕb) * dΩ
                            for d1 in 1:dim
                                #if dim*(b-1) + d2 >= dim*(a-1) + d1
                                Ke_0[dim*(a-1) + d1, dim*(b-1) + d2] += Ke_e[d1,d2]
                                #end
                            end
                        end
                    end
                end
            end
            if MatrixType <: SizedMatrix # Work around because full constructor errors
                push!(Kes, Symmetric(SizedMatrix{Kesize,Kesize,T}(Ke_0)))
            else
                push!(Kes, Symmetric(MatrixType(Ke_0)))
            end
        end
    else
        Kes = let Kesize=Kesize, nel=nel
            [Symmetric(zeros(T, Kesize, Kesize), :U) for i = 1:nel]
        end
        weights = let Kesize=Kesize, nel=nel
            [zeros(T, Kesize) for i = 1:nel]
        end
    
        Ke_e = zeros(T, dim, dim)
        
        celliteratortype = CellIterator{typeof(dh).parameters...}
        _celliterator = CellIterator(dh)
        for (k, cell) in enumerate(_celliterator)
            reinit!(cellvalues, cell)
            fe = weights[k]
            for q_point in 1:getnquadpoints(cellvalues)
                dΩ = getdetJdV(cellvalues, q_point)
                for b in 1:n_basefuncs
                    ∇ϕb = shape_gradient(cellvalues, q_point, b)
                    ϕb = shape_value(cellvalues, q_point, b)
                    for d2 in 1:dim
                        # Force business
                        fe[(b-1)*dim + d2] += ϕb * body_force[d2] * dΩ
                        for a in 1:n_basefuncs
                            ∇ϕa = shape_gradient(cellvalues, q_point, a)
                            Ke_e .= dotdot(∇ϕa, C, ∇ϕb) * dΩ
                            for d1 in 1:dim
                                #if dim*(b-1) + d2 >= dim*(a-1) + d1
                                Kes[k].data[dim*(a-1) + d1, dim*(b-1) + d2] += Ke_e[d1,d2]
                                #end
                            end
                        end
                    end
                end
            end
        end
    end
    return Kes, weights
end

function _make_dloads(fes, problem, facevalues)
    dim = getdim(problem)
    N = nnodespercell(problem)
    T = floattype(problem)
    dloads = deepcopy(fes)
    for i in 1:length(dloads)
        dloads[i] .= 0
    end
    pressuredict = getpressuredict(problem)
    dh = getdh(problem)
    grid = dh.grid
    boundary_matrix = grid.boundary_matrix
    cell_coords = zeros(JuAFEM.Vec{dim, T}, N)
    n_basefuncs = getnbasefunctions(facevalues)
    for k in keys(pressuredict)
        t = -pressuredict[k] # traction = negative the pressure
        faceset = getfacesets(problem)[k]
        for (cellid, faceid) in faceset
            boundary_matrix[faceid, cellid] || throw("Face $((cellid, faceid)) not on boundary.")
            fe = dloads[cellid]
            getcoordinates!(cell_coords, grid, cellid)
            reinit!(facevalues, cell_coords, faceid)
            for q_point in 1:getnquadpoints(facevalues)
                dΓ = getdetJdV(facevalues, q_point) # Face area
                normal = getnormal(facevalues, q_point) # Nomral vector at quad point
                for i in 1:n_basefuncs
                    ϕ = shape_value(facevalues, q_point, i) # Shape function value
                    for d = 1:dim
                        fe[(i-1)*dim + d] += ϕ * t * normal[d] * dΓ
                    end
                end
            end
        end
    end
    
    return dloads
end

function make_cload(problem)
    T = floattype(problem)
    dim = getdim(problem)
    cloads = getcloaddict(problem)
    dh = getdh(problem)
    metadata = getmetadata(problem)
    node_dofs = metadata.node_dofs
    inds = Int[]
    vals = T[]
    for nodeidx in keys(cloads)
        for (dofidx, force) in enumerate(cloads[nodeidx])
            if force != 0
                dof = node_dofs[(nodeidx-1)*dim+dofidx]
                push!(inds, dof)
                push!(vals, force)
            end
        end
    end
    return sparsevec(inds, vals, ndofs(dh))
end
