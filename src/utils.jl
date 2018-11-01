function find_varind(black, white, ::Type{TI}=Int) where TI
    nel = length(black)
    nel == length(white) || throw("Black and white vectors should be of the same length")
    varind = zeros(TI, nel)
    k = 1
    for i in 1:nel
        if !black[i] && !white[i]
            varind[i] = k
            k += 1
        end
    end
    return varind
end

function find_black_and_white(dh)
    black = BitVector(getncells(dh.grid))
    white = BitVector(getncells(dh.grid))
    black .= false
    white .= false
    if haskey(dh.grid.cellsets, "black")
        for c in grid.cellsets["black"]
            black[c] = true
        end
    end
    if haskey(dh.grid.cellsets, "white")
        for c in grid.cellsets["white"]
            white[c] = true
        end
    end
    
    return black, white
end

YoungsModulus(p) = getE(p)
PoissonRatio(p) = getν(p)

JuAFEM.getncells(problem::StiffnessTopOptProblem) = JuAFEM.getncells(getdh(problem).grid)

function compliance(Ke, u, dofs)
    comp = zero(eltype(u))
    for i in 1:length(dofs)
        for j in 1:length(dofs)
            comp += u[dofs[i]]*Ke[i,j]*u[dofs[j]]
        end
    end
    comp
end

function meandiag(K::AbstractMatrix)
    z = zero(eltype(K))
    for i in 1:size(K, 1)
        z += abs(K[i, i])
    end
    return z / size(K, 1)
end

density(var, xmin) = var*(1-xmin) + xmin

macro debug(expr)
    return quote
        if DEBUG[]
            $(esc(expr))
        end
    end
end
