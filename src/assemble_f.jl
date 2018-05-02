function assemble_f!(f::AbstractVector, problem, dloads)
    metadata = problem.metadata
    dof_cells = metadata.dof_cells
    dof_cells_offset = metadata.dof_cells_offset
    for dofidx in 1:ndofs(problem.ch.dh)
        r = dof_cells_offset[dofidx] : dof_cells_offset[dofidx+1]-1
        for i in r
            cellidx, localidx = dof_cells[i]
            f[dofidx] += dloads[cellidx][localidx]
        end
    end
    return f
end
