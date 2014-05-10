module NumericalMath

# package code goes here

  export horner,
         ridders, bisect, brent_dekker,
         trapz,
         lambert_W

  include("polynomial.jl")
  include("fzero.jl")
  include("integration.jl")
  include("special.jl")

end # module
