##
##  s p e c i a l . j l  Special Functions
##


# Lambert W function, inverse of x -> x exp(x)
# its derivative is 1 / (1+zw) / exp(zw) where zw = lambert_W(z)
#
function lambertW(x::Real)
    if x < -1/exp(1); return NaN; end
    w0 = 1.0
    w1 = w0 - (w0 * exp(w0) - x)/((w0 + 1) * exp(w0) - 
        (w0 + 2) * (w0 * exp(w0) - x)/(2 * w0 + 2))
    while abs(w1 - w0) > 1e-15
        w0 = w1
        w1 = w0 - (w0 * exp(w0) - x)/((w0 + 1) * exp(w0) - 
            (w0 + 2) * (w0 * exp(w0) - x)/(2 * w0 + 2))
    end
    return w1
end
