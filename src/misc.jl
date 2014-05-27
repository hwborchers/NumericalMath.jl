##
##  m i s c . j l
##


#-- Algebraic-geometric mean
function agm(a::Number, b::Number)
    local a0 = a, b0 = b
    local c  = max(abs(a), abs(b))
    if typeof(a) <: Real && a < 0; a0 = complex(a, 0); end
    if typeof(b) <: Real && b < 0; b0 = complex(b, 0); end

    while abs(b0-a0) >= eps(c)
        a0, b0 = (a0 + b0)/2.0, sqrt(a0 * b0)
    end
    return (a0 + b0) / 2.0
end


#-- Arc length of a curve defined through f:R --> R
function arc_length(f::Function, a::Real, b::Real;
                        nmax::Int=25, tol::Real=2.5e-15)
    local k::Int
    h = b - a;
    fa = f(a); fb = f(b)
    A = zeros(nmax, nmax)
    A[1, 1] = sqrt(sum((fb - fa)^2))
    for k = 1:(nmax-1)
        h /= 2
        x = linspace(a, b, 2^k + 1)
        y = map(f, x)
        X = [x y]
        dX = diff(X)
        A[k+1, 1] = sum(sqrt(sum(dX.^2, 2)))
        for (l in 1:k)
            A[k+1, l+1] = (4^l * A[k+1, l] - A[k, l]) / (4^l - 1)
        end
        if abs(A[k+1, k+1] - A[k, k]) < tol && k > 3
            break
        end
    end
    err = abs(A[k+1, k+1] - A[k, k])
    return A[k+1, k+1], k+1, err
end
