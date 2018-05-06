struct Metadata
    cell_dofs::Matrix{Int}
    dof_cells::Vector{Tuple{Int,Int}}
    dof_cells_offset::Vector{Int}
    #node_first_cells::Vector{Tuple{Int,Int}}
    node_cells::Vector{Tuple{Int,Int}}
    node_cells_offset::Vector{Int}
    node_dofs::Matrix{Int}
end

function Metadata(dh::DofHandler{dim}) where dim
    cell_dofs = get_cell_dofs_matrix(dh)
    dof_cells, dof_cells_offset = get_dof_cells_matrix(dh, cell_dofs)
    #node_first_cells = get_node_first_cells(dh)
    node_cells, node_cells_offset = get_node_cells(dh)
    node_dofs = reshape(dh.node_dofs, dim, getnnodes(dh.grid))

    #meta = Metadata(cell_dofs, dof_cells, dof_cells_offset, node_first_cells, node_dofs)
    meta = Metadata(cell_dofs, dof_cells, dof_cells_offset, node_cells, node_cells_offset, node_dofs)
end

function get_cell_dofs_matrix(dh)
    cell_dofs = zeros(Int, ndofs_per_cell(dh), getncells(dh.grid))
    for i in 1:size(cell_dofs, 2)
        r = dh.cell_dofs_offset[i]:dh.cell_dofs_offset[i+1]-1
        for j in 1:length(r)
            cell_dofs[j,i] = dh.cell_dofs[r[j]]
        end
    end
    cell_dofs
end

function get_dof_cells_matrix(dh, cell_dofs)
    dof_cells_vecofvecs = [Vector{Tuple{Int,Int}}() for i in 1:ndofs(dh)]
    l = 0
    for cellidx in 1:size(cell_dofs, 2)
        for localidx in 1:size(cell_dofs, 1)
            dofidx = cell_dofs[localidx, cellidx]
            push!(dof_cells_vecofvecs[dofidx], (cellidx, localidx))
            l += 1
        end
    end

    dof_cells = Tuple{Int,Int}[]
    sizehint!(dof_cells, l)
    dof_cells_offset = Int[1]
    sizehint!(dof_cells_offset, ndofs(dh)+1)

    for (dofidx, indsincells) in enumerate(dof_cells_vecofvecs)
        append!(dof_cells, indsincells)
        push!(dof_cells_offset, length(dof_cells)+1)
    end
    
    dof_cells, dof_cells_offset
end

function get_node_first_cells(dh)
    node_first_cells = fill((0,0), getnnodes(dh.grid))
    visited = BitVector(getnnodes(dh.grid))
    visited .= false
    for cellidx in 1:getncells(dh.grid)
        for (local_node_idx, global_node_idx) in enumerate(dh.grid.cells[cellidx].nodes)
            if !visited[global_node_idx]
                visited[global_node_idx] = true
                node_first_cells[global_node_idx] = (cellidx, local_node_idx)
            end
        end
    end
    return node_first_cells
end

function get_node_cells(dh)
    node_cells_vecofvecs = [Vector{Tuple{Int,Int}}() for i in 1:ndofs(dh)]
    l = 0
    for (cellidx, cell) in enumerate(CellIterator(dh))
        for (localidx, nodeidx) in enumerate(cell.nodes)
            push!(node_cells_vecofvecs[nodeidx], (cellidx, localidx))
            l += 1
        end
    end

    node_cells = Tuple{Int,Int}[] 
    node_cells_offset = Int[1]
    sizehint!(node_cells, l)
    sizehint!(node_cells_offset, getnnodes(dh.grid) + 1)

    for (nodeidx, indsincells) in enumerate(node_cells_vecofvecs)
        append!(node_cells, indsincells)
        push!(node_cells_offset, length(node_cells)+1)
    end

    return node_cells, node_cells_offset
end

function get_node_dofs(dh::DofHandler{dim}, node_first_cells) where dim
    node_dofs = fill(0, dim, getnnodes(dh.grid))
    xh = zeros(typeof(dh.grid.nodes[1].x), length(dh.grid.cells[1].nodes))
    _cell_dofs = zeros(Int, ndofs_per_cell(dh))
    bcv = dh.bc_values[1]
    interpolation = dh.field_interpolations[1]
    for nodeidx in 1:size(node_dofs, 2)
        found = false
        node_coord = dh.grid.nodes[nodeidx].x
        cellidx = node_first_cells[nodeidx][1]
        getcoordinates!(xh, dh.grid, cellidx)
        for (faceidx, face_field_points) in enumerate(JuAFEM.faces(interpolation))
            bcv.current_face[] = faceidx
            for p in 1:length(face_field_points)
                field_point_coord = spatial_coordinate(bcv, p, xh)
                if node_coord ≈ field_point_coord
                    found = true
                    celldofs!(_cell_dofs, dh, cellidx)
                    offset = (face_field_points[p]-1)*dim
                    for d in 1:dim
                        index_in_cell = offset + d
                        node_dofs[d, nodeidx] = _cell_dofs[index_in_cell]
                    end
                    break
                end
            end
            if found
                break
            end        
        end
    end
    return node_dofs
end
