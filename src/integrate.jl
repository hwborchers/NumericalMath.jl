##
##  i n t e g r a t i o n . j l  Integration Routines
##


function trapz{Tx<:Number, Ty<:Number}(x::Vector{Tx}, y::Vector{Ty})
    # Trapezoidal integration rule
    local n = length(x)
    if (length(y) != n)
        error("Vectors 'x', 'y' must be of same length")
    end
    r = zero(zero(Tx) + zero(Ty))
    if n == 1; return r; end
    for i in 2:n
        r += (x[i] - x[i-1]) * (y[i] + y[i-1])
    end
    #= correction -h^2/12 * (f'(b) - f'(a))
    ha = x[2] - x[1]
    he = x[end] - x[end-1]
    ra = (y[2] - y[1]) / ha
    re = (y[end] - y[end-1]) / he
    r/2 - ha*he/12 * (re - ra)
    =#
    return r/2
end


function romberg(f::Function, a::Real, b::Real;
                   maxit = 50, tol = 1.0e-15)
    if a == b; return 0.0, 0.0; end
    sgn = 1.0
    if a > b
        a, b = b, a
        sgn = -1.0
    end

    local n::Int = 1, iter::Int = 0
    I = zeros(maxit+1, maxit+1)

    err = 1.0
    while err > tol && iter < maxit
        iter += 1
        n *= 2
        h = (b - a)/ n
        S = f(a)
        for i = 1:(n-1)
            xi = a + h*i
            S += 2.0*f(xi)
        end
        I[iter+1, 1] = (S + f(b)) * h / 2.0
        for k = 2:(iter+1)
            j = 2 + iter - k
            I[j,k] = (4^(k-1)*I[j+1,k-1] - I[j,k-1]) / (4^(k-1)-1)
        end
        err = abs((I[1,iter+1] - I[2,iter]) / I[1,iter+1])
     end
     if iter >= maxit
         warn("Maximum number of iterations has been reached.")
     end
     if err == 0.0; err = tol; end

     return I[1, iter+1], err
end

