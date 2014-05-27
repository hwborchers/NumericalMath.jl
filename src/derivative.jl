##
##  d e r i v a t i v e . j l  Numerical Derivatives
##


#-- Numerical gradient of multivariate function ------------------------------
#
function fd_gradient{T<:Real}(f::Function, x0::Array{T,1}; h::Real = 0.0)
    local heps::Real = h
    if heps == 0.0; heps = eps()^(1.0/3.0); end
    local n::Int = length(x0)
    local hh = zeros(n), gr = zeros(n)
    for i = 1:n
        hh[i] = heps
        gr[i] = (f(x0 + hh) - f(x0 - hh)) / (2.0*heps)
        hh[i] = 0.0
    end
    return gr
end


function fd_jacobian{T<:Real}(f::Function, x0::Array{T,1}; h::Real = 0.0)
    local heps::Real = h
    if heps == 0.0; heps = eps()^(1.0/3.0); end
	local n = length(x0), m = length(f(x0))
	local hh = zeros(n)
	jacob = zeros(m, n)
	for i = 1:n
		hh[i] = heps
		jacob[:, i] = (f(x0 + hh) - f(x0 - hh)) / (2.0*heps)
		hh[i] = 0.0
	end
	return jacob
end


function fd_hessian{T<:Real}(f::Function, x0::Array{T,1}; h::Real = 0.0)
    local heps::Real = h
    if heps == 0.0; heps = eps()^(1.0/4.0); end
    local n::Int = length(x0)
    if length(f(x0)) != 1
        error("Function 'f' must a scalar function of n variables.")
    end

    if n == 1
        H = (f(x0.-heps) - 2.0*f(x0) + f(x0.+heps)) / heps^2

    else
        H = zeros(n, n)
        hh = diagm(heps*ones(n))
        for i = 1:(n-1)
            hi = hh[:, i]
            H[i, i] = (f(x0-hi) - 2.0*f(x0) + f(x0+hi)) / heps^2
            for j = (i+1):n
                hj = hh[:, j]
                H[i, j] = (f(x0+hi+hj) - f(x0+hi-hj) - f(x0-hi+hj) + f(x0-hi-hj)) / (4.0*heps^2)
                H[j, i] = H[i, j]
            end
        end
        hi = hh[:, n]
        H[n, n] = (f(x0-hi) - 2.0*f(x0) + f(x0+hi)) / heps^2
    end

    return H
end


function fd_laplacian{T<:Real}(f::Function, x0::Array{T,1}; h::Real = 0.0)
    local heps::Real = h
    if heps == 0.0; heps = eps()^(1.0/4.0); end
    local n::Int = length(x0)
    local hh::Real = zeros(n)

    L = 0.0
    for i = 1:n
        hh[i] = heps
        L += (f(x0+hh) + f(x0-hh) - 2.0*f(x0)) / heps^2
        hh[i] = 0.0
    end

    return L
end


#-- Central difference combined with Richardson approximation ----------------
#
function numderiv(f::Function, x0::Real; n::Int = 16, h::Real = 0.1)

    local epsilon = eps(), err = Inf
    local hstep::Real = h

    j::Int = 1
    D = zeros(n, n)
    D[1, 1] = (f(x0+hstep) - f(x0-hstep))/(2*hstep)

    while j < n
        hstep = hstep / 2.0
        D[j+1, 1] = (f(x0+hstep) - f(x0-hstep)) / (2*hstep)
        for k = 1:j
            D[j+1,k+1] = D[j+1,k] + (D[j+1,k] - D[j,k]) / (4^k - 1)
        end

        err_new = 2 * abs(D[j+1,j+1] - D[j,j]) /
                  (abs(D[j+1,j+1]) + abs(D[j,j]) + epsilon)

        if err_new >= err; break; end
        err = err_new
        j += 1
    end

    return D[j,j]
end


#-- Complex-step derivative approach -----------------------------------------
#
function complex_step(f::Function, x0::Real; h::Real = 1e-20)
    # Check whether function accepts and returns complex values
    if !isa(f(x0), Real)
        error("Function 'f' must take on a real(!) value at 'x0'.")
    end

    fx0hi = try
        f(x0 + h*1im)
    catch err
        error("Function 'f' does not appear to accept complex arguments.")
    end

    # f(x0) must be real, fx0hi complex
    if !isa(fx0hi, Complex)
        error("'f(x0 + h*1im)' must be complex for 'complex-step' to work.")
    end

    return imag(fx0hi) / h
end
