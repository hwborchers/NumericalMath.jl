##
##  l a . j l  Linear algebra routines
##


function trisolve{T<:Number}(d0::Vector{T},
                             d1::Vector{T},d2::Vector{T},rhs::Vector{T})
    local n = length(d0)
    if n < 3
        error("Length of 'a' must be greater or equal to 3.")
    end
    if length(d1) != n-1 || length(d2) != n-1
        error("Vectors 'dup' and 'dlo' must be of a length of length(a)-1.")
    end
    if length(rhs) != n
        error("Vector 'rhs' must be of the same length as 'a'.")
    end

    local x = copy(rhs), a = copy(d0), b = [d1, 0], d = copy(d2)
    for i = 1:(n-1)
        if d[i] != 0
            t = a[i]/d[i]
            si = 1/sqrt(1+t*t)
            co = t*si
            a[i] = a[i]*co + d[i]*si

            h = b[i]
            b[i] = h*co + a[i+1]*si
            a[i+1] = -h*si + a[i+1]*co
            d[i] = b[i+1]*si
            b[i+1] = b[i+1]*co

            h = x[i]
            x[i] = h*co + x[i+1]*si
            x[i+1] = -h*si + x[i+1]*co
        end
    end

    if any(a .== 0.0)
        error("Triangular matrix is singular---system not solvable.")
    end

    x[n] = x[n]/a[n]
    x[n-1] = ( x[n-1] - b[n-1]*x[n] ) / a[n-1]
    for i = (n-2):-1:1
        x[i] = (x[i] - b[i]*x[i+1] - d[i]*x[i+2]) / a[i]
    end

    return x
end
