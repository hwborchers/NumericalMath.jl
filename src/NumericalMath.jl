module NumericalMath

# package code goes here

  export trisolve
  export pval, horner, pzero
  export interp1d, pchip
  export ridders, brent_dekker
  export fd_gradient, fd_jacobian, fd_hessian, fd_laplacian,
         numderiv, complex_step
  export trapz, romberg, line_integral
  export agm, arc_length
  export lambertW
  
  include("la.jl")
  include("polynom.jl")
  include("interpolate.jl")
  include("fzero.jl")
  include("derivative.jl")
  include("integrate.jl")
  include("misc.jl")
  include("specfun.jl")

end # module
