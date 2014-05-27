##
##  s p e c i a l . j l  Special Functions
##


# Lambert W function, inverse of x -> x exp(x)
# its derivative is 1 / (1+zw) / exp(zw) where zw = lambertW(z)
#
function lambertW(x::Real)
    if x <  -1.0/exp(1.0); return NaN, NaN;  end
    if x == -1.0/exp(1.0); return -Inf, Inf; end

    local w0::Real = 1.0, w1::Real, w2::Real
    w1 = w0 - (w0 * exp(w0) - x)/((w0 + 1.0) * exp(w0) - 
        (w0 + 2.0) * (w0 * exp(w0) - x)/(2.0 * w0 + 2.0))

    while abs(w1 - w0) > 1e-15
        w0 = w1
        w1 = w0 - (w0 * exp(w0) - x)/((w0 + 1.0) * exp(w0) - 
            (w0 + 2.0) * (w0 * exp(w0) - x)/(2.0 * w0 + 2.0))
    end

    w2 = 1.0 / (1 + w1) / exp(w1)

    return w1, w2
end
