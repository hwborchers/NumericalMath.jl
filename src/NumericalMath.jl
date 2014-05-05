module NumericalMath

# package code goes here

  export pval, horner
  export ridders, bisect, brent_dekker
  export trapz
  export lambert_W

  include("polynomial.jl")
  include("fzero.jl")
  include("integration.jl")
  include("special.jl")


end # module
