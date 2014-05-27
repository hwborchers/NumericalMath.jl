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


function pzero{T<:Real}(p::Vector{T}, x0::T)
    local x = x0, tol = 2*eps(x)
    px, dpx = horner(p, x)
    df = -px/dpx
    niter = 1
    while abs(df) >= tol && niter <= 100
        x += df
        px, dpx = horner(p, x)
        df = -px/dpx
        niter += 1
    end
    if niter > 100
        warn("Number of iterations exceeded.")
        x = NaN
    end
    return x + df
end


function pfit()
    error("Polynomial fit not yet implemented.")
end
