function inp_to_juafem(filepath_with_ext)
    problem = import_inp(filepath_with_ext)
    return inp_to_juafem(problem)
end
function inp_to_juafem(problem::InpContent)
    _celltype = problem.celltype
    if _celltype == "C3D10"
        celltype = JuAFEM.QuadraticTetrahedron
        geom_order = 2
        refshape = RefTetrahedron
        dim = 3
        # Add more cell types
    end
    cells = celltype.(problem.cells)
    nodes = Node.(problem.node_coords)
    grid = Grid(cells, nodes)

    for k in keys(problem.cellsets)
        grid.cellsets[k] = Set(problem.cellsets[k])
    end
    for k in keys(problem.nodesets)
        grid.nodesets[k] = Set(problem.nodesets[k])
    end
    for k in keys(problem.facesets)
        grid.facesets[k] = Set(problem.facesets[k])
    end
    # Define boundary faces
    grid.boundary_matrix = extract_boundary_matrix(grid);

    dh = DofHandler(grid)
    # Isoparametric
    field_interpolation = Lagrange{dim, refshape, geom_order}()
    push!(dh, :u, dim, field_interpolation) # Add a displacement field
    close!(dh)

    ch = ConstraintHandler(dh)
    for k in keys(problem.nodedbcs)
        vec = problem.nodedbcs[k]
        f(x, t) = [vec[i][2] for i in 1:length(vec)]
        components = [vec[i][1] for i in 1:length(vec)]
        dbc = Dirichlet(:u, getnodeset(grid, k), f, components)
        add!(ch, dbc)
        close!(ch)
        update!(ch, 0.0)
    end

    black, white = find_black_and_white(dh)
    varind = find_varind(black, white)
    metadata = Metadata(dh)

    return InpStiffness(problem, Val{geom_order}, ch, black, white, varind, metadata)
end

@static if VERSION < v"0.7.0-DEV.2563"
    const ht_keyindex2! = Base.ht_keyindex2
else
    import Base.ht_keyindex2!
end

function extract_boundary_matrix(grid::Grid{dim}) where dim
    nfaces = length(JuAFEM.faces(grid.cells[1]))
    ncells = length(grid.cells)
    countedbefore = Dict{NTuple{dim,Int},Bool}()
    boundary_matrix = ones(Bool, nfaces, ncells) # Assume all are boundary faces
    for (ci, cell) in enumerate(getcells(grid))    
        for (fi, face) in enumerate(JuAFEM.faces(cell))
            sface = JuAFEM.sortface(face) # TODO: faces(cell) may as well just return the sorted list
            token = ht_keyindex2!(countedbefore, sface)
            if token > 0 # haskey(countedbefore, sface)
                boundary_matrix[fi, ci] = 0
            else # distribute new dofs
                Base._setindex!(countedbefore, true, sface, -token)# countedbefore[sface] = true,  mark the face as counted
            end
        end
    end
    sparse(boundary_matrix)
end
