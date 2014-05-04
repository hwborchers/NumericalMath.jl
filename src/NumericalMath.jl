module NumericalMath

# package code goes here

  export ridders, bisect, brent_dekker
  export trapz
  export lambert_W

  include("fzero.jl")
  include("integration.jl")
  include("special.jl")


end # module
