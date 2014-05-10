##
##  p o l y n o m i a l . j l  Polynomials
##


function horner{T<:Number}(p::Vector{T}, x::T)
    local n = length(p)
    if n == 0
        return 0,0
    elseif n == 1
        return p[1],0
    else
        y = p[1]; dy = 0
        for i in 2:n
            dy = dy * x + y
            y  =  y * x + p[i]
        end
    end
    return y,dy
end


function horner{T<:Number}(p::Vector{T}, x::Vector{T})
    local n = length(p)
    if n == 0
        return [0],[0]
    elseif n == 1
        return p,[0]
    else
        y = p[1]; dy = 0
        for i in 2:n
            dy = dy .* x .+ y
            y  =  y .* x .+ p[i]
        end
    end
    return y,dy
end
