##
##  p o l y n o m i a l . j l  Polynomials
##


function pval(p::Vector, x::Number)
    local n = length(p)
    if n == 0
        return NaN
    elseif n == 1
        return p[1]
    else
        y = p[1]
        for i in 2:n
            y = y * x + p[i]
        end
    end
    return y
end


function horner(p::Vector, x::Number)
    local n = length(p)
    if n == 0
        return NaN, NaN
    elseif n == 1
        return p[1], 0
    else
        y = p[1]; dy = 0
        for i in 2:n
            dy = dy * x + y
            y  =  y * x + p[i]
        end
    end
    return y, dy
end
