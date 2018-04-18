
struct ElementFEAInfo{T, TK<:AbstractMatrix{T}, Tf<:AbstractVector{T}, TKes<:AbstractVector{TK}, Tfes<:AbstractVector{Tf}}
    Kes::TKes
    fes::Tfes
end

struct GlobalFEAInfo{T, TK<:AbstractMatrix{T}, Tf<:AbstractVector{T}}
    K::TK
    f::Tf
end
GlobalFEAInfo(K::AbstractMatrix{T}, f::AbstractVector) where T = GlobalFEAInfo{T, typeof(K), typeof(f)}(K, f)
GlobalFEAInfo(::Type{T}) where T = GlobalFEAInfo{T}()
GlobalFEAInfo() = GlobalFEAInfo{Float64}()
GlobalFEAInfo{T}() where T = GlobalFEAInfo{T, SparseMatrixCSC{T, Int}, Vector{T}}(sparse(zeros(T, 0, 0)), zeros(T, 0))

make_empty_K(sp::StiffnessTopOptProblem) = Symmetric(create_sparsity_pattern(sp.ch.dh))

make_empty_f(sp::StiffnessTopOptProblem{dim, T}) where {dim, T} = zeros(T, ndofs(sp.ch.dh))

function make_Kes(sp::StiffnessTopOptProblem, quad_order=2)
    make_Kes(sp, quad_order, Val{:Static})
end
function make_Kes(sp::StiffnessTopOptProblem, ::Type{Val{mat_type}}) where mat_type
    make_Kes(sp, 2, Val{mat_type})
end
function make_Kes(sp::StiffnessTopOptProblem{dim, T}, quad_order, ::Type{Val{mat_type}}) where {dim, T, mat_type}
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
    
    return _make_Kes(sp, Val{mat_type}, Val{n_basefuncs}, Val{dim*n_basefuncs}, C, interpolation_space, quadrature_rule, cellvalues, facevalues)
end
function _make_Kes(sp::StiffnessTopOptProblem{dim, T}, ::Type{Val{mat_type}}, ::Type{Val{n_basefuncs}}, ::Type{Val{ndofs_per_cell}}, C, interpolation_space, quadrature_rule, cellvalues, facevalues) where {dim, T, mat_type, n_basefuncs, ndofs_per_cell}
    dh = sp.ch.dh

    # Calculate element stiffness matrices
    Kesize = ndofs_per_cell
    nel = getncells(dh.grid)

    if mat_type === :Static || mat_type === :SMatrix
        Kes = Symmetric{T, SMatrix{Kesize, Kesize, T, Kesize^2}}[]
        resize!(Kes, nel)

        Ke_e = zeros(T, dim, dim)
        Ke_0 = Matrix{T}(Kesize, Kesize)
        celliteratortype = CellIterator{typeof(dh).parameters...}
        _celliterator::celliteratortype = CellIterator(dh)
        for (k, cell) in enumerate(_celliterator)
            Ke_0 .= 0
            reinit!(cellvalues, cell)
            for q_point in 1:getnquadpoints(cellvalues)
                for b in 1:n_basefuncs
                    for a in 1:n_basefuncs
                        ∇ϕa = shape_gradient(cellvalues, q_point, a)
                        ∇ϕb = shape_gradient(cellvalues, q_point, b)
                        Ke_e .= dotdot(∇ϕa, C, ∇ϕb) * getdetJdV(cellvalues, q_point)
                        for d2 in 1:dim
                            for d1 in 1:dim
                                #if dim*(b-1) + d2 >= dim*(a-1) + d1
                                Ke_0[dim*(a-1) + d1, dim*(b-1) + d2] += Ke_e[d1,d2]
                                #end
                            end
                        end
                    end
                end
            end
            push!(Kes, Symmetric(SMatrix{Kesize, Kesize, T, Kesize*Kesize}(Ke_0)))
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
                for b in 1:n_basefuncs
                    for a in 1:n_basefuncs
                        ∇ϕa = shape_gradient(cellvalues, q_point, a)
                        ∇ϕb = shape_gradient(cellvalues, q_point, b)
                        Ke_e .= dotdot(∇ϕa, C, ∇ϕb) * getdetJdV(cellvalues, q_point)
                        for d2 in 1:dim
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

function make_fes(sp::Union{PointLoadCantilever, HalfMBB}, ::Type{Val{vec_type}}=Val{:Matrix}) where {vec_type}
    _ndofs_per_cell = ndofs_per_cell(sp.ch.dh)
    nel = getncells(sp.ch.dh.grid)
    # Function barrier
    _make_fes(sp, Val{vec_type}, Val{_ndofs_per_cell}, nel)
end
function _make_fes(sp::Union{PointLoadCantilever{dim, T}, HalfMBB{dim, T}}, ::Type{Val{vec_type}}, ::Type{Val{ndofs_per_cell}}, nel) where {dim, T, vec_type, ndofs_per_cell}
    dof_cells = sp.metadata.dof_cells
    fdof = sp.force_dof
    force = sp.force

    counter = 0
    for i in 1:size(dof_cells, 1)
        cellidx, localidx = dof_cells[i,fdof]
        if (cellidx, localidx) != (0, 0)
            counter += 1
        else
            break
        end
    end

    if vec_type === :SVector || vec_type === :Static
        fes = [zeros(SVector{ndofs_per_cell, T}) for i in 1:nel]
    elseif vec_type === :MVector
        fes = [zeros(MVector{ndofs_per_cell, T}) for i in 1:nel]
    else
        fes = zeros(T, ndofs_per_cell, nel)        
    end

    force_per_cell = force/counter
    for i in 1:size(dof_cells, 1)
        cellidx, localidx = dof_cells[i,fdof]
        if (cellidx, localidx) != (0, 0)
            if vec_type === :SVector || vec_type === :Static || vec_type === :MVector
                fes[cellidx] = ntuple((i)->(i == localidx ? force_per_cell : fes[cellidx][i]), Val{ndofs_per_cell})
            else
                for j in 1:ndofs_per_cell
                    if j == localidx
                        fes[j, cellidx] = force_per_cell
                        break
                    end
                end
            end
        end
    end

    return fes
end
