#   NumericalMath

##  Introduction

Functions from Numerical Mathematics for Julia.

##  Functions

### Linear algebra

  * trisolve(d0, d1, d2, rhs) -- solves triangular systems of equations

### Polynomials

  * pval(p, x)      -- evaluates polynomial p in x (Matlab-like order of coefficients)
  * horner(p, x)    -- Horner schema returning the value of polynomial p at x
                         and the value of its derivative p' at x
  * pzero(p, x0)    -- finds a root of polynomial p near x0 applying Newton-Raphson
  * pfit(xi, yi, n) -- polynomial fitting of data points (of order n)

### Root finding

  * ridders(f, a, b)      -- Ridders' method for finding roots of univariate functions
  * brent_dekker(f, a, b) -- root finding using the Brent-Dekker approach

### Interpolation and approximation

  * interp1d(xs, ys, x) -- Interpolation of data points with methods
                             :constant, :nearest, :linear, :spline, :cubic
  * pchip(xs, ys, x)    -- Piecewise cubic hermitean interpolating polynomial

### Differentiation

  * fd_gradient(f, x0; h)  -- numerical gradient of a multivariate function
                                applying the "central difference formula"
  * fd_jacobian(f, x0; h)  -- Jacobian matrix of multivariate function at x0
  * fd_hessian(f, x0; h)   -- numerical Hessian, based on finite differences
  * fd_laplacian(f, x0; h) -- numerical Laplacian, based on finite differences
  * numderiv(f, x0; n, h)  -- Richardson method applied to central difference
  * complex_step(f, x0; h) -- complex-step derivative

### Integration

  * trapz(x, y)      -- Trapezoidal rule for integrating discrete points
                          (with end point correction terms)
  * romberg(f, a, b) -- Romberg integration (ie., utilizes Richardson's method)
  * line_integral(f, points)    -- line integral of complex functions
                                     along the path defined by 'points'

### Miscellaneous

  * agm(a, b)           -- algebraic-geometric mean of numbers a, b
  * arc_length(f, a, b) -- arc length of the curve defined by f:[a,b] --> R

### Special functions

  * lambertW(x)    -- Lambert W function and its first derivative
  * legendre(n, x) -- Legendre polynomials of degree <= n at x

## Examples
