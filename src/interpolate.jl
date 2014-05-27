##
##  i n t e r p o l a t e . j l  Interpolation
##


function interp1d{T<:Real}(xi::Array{T,1}, yi::Array{T,1}, x::T;
                           method::Symbol = :cubic)
    n = length(xi)
    if length(yi) != n
        error("Vectors 'xi', 'yi' must be of same length.")
    end
    if n == 1
       error("At least two data points must be provided.")
    end
    if any(xi[1:end-1] .> xi[2:end])
        error("Elements in 'xi' must be strictly increasing.")
    end

    local k::Int = 0
    if x < xi[1]; k = 0
    elseif x >= xi[end]; k = n
    else
        k = maximum(find(xi[1:(end-1)] .<= x))
    end
    
    if method == :constant
        if k == 0; k = 1; end
        return yi[k]

    elseif method == :nearest
        if k == 0; k = 1
        elseif k < n
            if x >= (xi[k] + xi[k+1]) / 2.0
                k += 1
            end
        end
        return yi[k]
        
    elseif method == :linear
        if k == 0; k = 1
        elseif k == n; k = n-1; end
        s = x - xi[k]
        delta = (yi[k+1] - yi[k]) / (xi[k+1] - xi[k])
        return yi[k] + (x - xi[k])*delta

    elseif method == :spline
        error("Method $(method) is not yet implemented.")

    elseif method == :cubic
        return pchip(xi, yi, x)

    else
        error("Unknown interpolation method: $(method)")
    end
end


function interp1d{T<:Real}(xi::Array{T,1}, yi::Array{T,1}, x::Array{T,1};
                           method::Symbol = :cubic)
    local m::Int = length(x)
    y = zeros(m)
    for i = 1:m
        y[i] = interp1d(xi, yi, x[i], method = method)
    end
    return y
end


#-- Moler: Piecewise Cubic Hermitean Interpolation Polynomials
#
function pchip{T<:Real}(xi::Array{T,1}, yi::Array{T,1}, x::T)
    n = length(xi)
    if length(yi) != n
        error("Vectors 'xi', 'yi' must be of same length.")
    end
    if any(xi[1:end-1] .> xi[2:end])
        error("Elements in 'xi' must be strictly increasing.")
    end

    # First derivatives
    h = diff(xi)
    delta = diff(yi) ./ h
    d = pchip_slopes(h, delta)

    # Piecewise polynomial coefficients
    a = (3.0*delta .- 2.0*d[1:(n-1)] .- d[2:n]) ./ h
    b = (d[1:(n-1)] .- 2.0*delta .+ d[2:n]) ./ h.^2

    # Find subinterval index k so that xi[k] <= x < xi[k+1]
    if x < xi[1]
        k = 1
    else
        k = maximum(find(xi[1:(n-1)] .<= x))
    end

    # Evaluate interpolant
    s = x - xi[k]
    v = yi[k] + s*(d[k] + s*(a[k] + s*b[k]))

    return v
end

function pchip_slopes(h, delta)

    # Slopes at interior points
    local n = length(h) + 1
    d = zeros(n)
    k = find(sign(delta[1:(n-2)]) .* sign(delta[2:(n-1)]) .> 0.0) .+ 1
    w1 = 2.0*h[k] .+ h[k.-1]
    w2 = h[k] .+ 2.0*h[k.-1]
    d[k] = (w1+w2) ./ (w1./delta[k.-1] + w2./delta[k])

    # Slopes at endpoints
    d[1] = pchip_end(h[1], h[2], delta[1], delta[2])
    d[n] = pchip_end(h[n-1], h[n-2], delta[n-1], delta[n-2])

    return d
end

function pchip_end(h1, h2, del1, del2)

    # Noncentered, shape-preserving, three-point formula.
    d = ((2.0*h1 + h2)*del1 - h1*del2) / (h1 + h2)
    if sign(d) != sign(del1)
        d = 0
    elseif sign(del1) != sign(del2)  &&  abs(d) > abs(3*del1)
        d = 3.0*del1
    end

    return d
end
