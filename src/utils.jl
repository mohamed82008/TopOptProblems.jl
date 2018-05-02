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

YoungsModulus(p::RectilinearPointLoad) = p.E
YoungsModulus(inp::InpStiffness) = inp.inp_content.E
