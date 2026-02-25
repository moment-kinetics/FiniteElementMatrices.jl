# FiniteElementMatrices

A package for computing unassembled matrices appearing in $`C^0`$ continuous Galerkin
finite element models.

# Matrices in one coordinate (1D)

This package computes matrices of the following form.
```math
A_{ij} = \int^{x_u}_{x_l} P_i(x) Q_j(x) \rho(x) d x,
```
where $`P_i`$ and $`Q_i`$ are from the set of Lagrange polynomials basis
functions (and their first derivatives) used for
$`C^0`$ continuous Galerkin
finite element models, $`x_u`$ and $`x_l`$
are the upper and lower limits of the element in the
coordinate $`x`$, and $`\rho(x)`$ is a function. We provide an interface
where the user can specify the function $`\rho(x)`$ as an argument, or
the user can pass an integer $`p`$, in which case $`\rho(x) = x^p`$.

The basis functions we use here are
```math
\Phi(x) = \Theta(x-x_l)\Theta(x_u - x)l_j(z(x)),
```
where $`\Theta(x)`$ is the Heaviside function and $`l_j(z)`$
is the $`j^{\rm th}`$ Lagrange poylnomial defined by the reference collocation points $`\{z_j\}`$, where
```math
 z = (x - c)/s,
```
and $`s`$ and $`c`$ are the coordinate scale and shift factors,
respectively. This means that $`P_i`$ and $`Q_i`$ can be either
of $`\Phi(x)`$ or $`d\Phi/dx`$, depending on which matrix
is required for the problem of interest.

To compute the matrix $`A_{ij}`$, we use the reference coordinate
$`z`$ to write
```math
A_{ij} = \int^{1}_{-1} S_i(z) T_j(z) \rho(s z + c) s d z,
```
where $`S_i(z)`$, $`T_i(z)`$ may be either of $`l_i(z)`$ or $`(1/s)d l_i(z)/d z`$.
We can also compute nonlinear matrices which have the form
```math
N_{ijk} = \int^{x_u}_{x_l} P_i(x) Q_j(x) R_k(x) \rho(x) d x = \int^{1}_{-1} S_i(z) T_j(z) U_k(z) \rho(s z + c) s d z.
```
where $`R_i(x)`$ may be either of $`\Phi(x)`$ or $`d\Phi/dx`$ and
$`U_i(z)`$ may be either of $`l_i(z)`$ or $`(1/s)d l_i(z)/d z`$.

# p-adaptive Gaussian quadrature

In order to guarantee that the matrices $`A_{ij}`$ with non-polynomial kernels $`\rho(x)`$ are calculated with sufficient accuracy, we use p-adaptive quadrature with Gauss-Legendre quadrature rules from `FastGaussQuadrature` at each iteration. We iteratively calculate the matrices with an increasingly large number of quadrature points and iterate until the result is constant up to absolute and relative tolerances. The inputs used to control this process are described in the section on usage below. For polynomial kernels we can deduce the number of quadrature points for exact results, avoiding the need for the iterative calculation in the pure polynomial case.

# Usage

To create these matrices, the user must supply the key grid
information for the element of interest, i.e., the set
$`\{z_j\}`$, and the factors $`s`$ and $`c`$.
This is done in the following function call, defining a
type `ElementCoordinates`. We use [FastGaussQuadrature](https://juliaapproximation.github.io/FastGaussQuadrature.jl/stable/)
to make a set of reference nodes, but any set of nodes with good
interpolation properties would be adequate.
```
using FastGaussQuadrature: gausslobatto
using FiniteElementMatrices
ngrid = 5 # change this to change the order of the higher-order FEM scheme
x, w = gausslobatto(ngrid)
scale = 1.0 # change these to match the size and position of elements in the problem
shift = 0.0
# the coordinate information
coordinate = ElementCoordinates(x,scale,shift)
```
Once an instance of `ElementCoordinates` is available,
we can compute the higher-order finite element matrices
as follows. We use the enums
`lagrange_x`, and `d_lagrange_dx` and the integer `p` to choose the form of the desired matrix. For a mass matrix, with
Jacobian $`1`$ for example,
we can use the following commands.
```
using FiniteElementMatrices
coordinate = ElementCoordinates(x,scale,shift)
M = finite_element_matrix(lagrange_x,lagrange_x,0,coordinate)
```
Note that the order in which the `@enum` Lagrange polynomials are given corresponds to the indices of `M`. E.g., for a matrix describing a first derivative
```
l_i = lagrange_x
l_j = d_lagrange_dx
P = finite_element_matrix(l_i,l_j,0,coordinate)
```
the number `P[i,j]` corresponds to the matrix element
```math
P_{ij} = \int^{x_u}_{x_l} \Phi_i(x) \frac{d \Phi_j}{d x} d x.
```
For a stiffness matrix describing a second derivative, we could
use the following.
```
K = -finite_element_matrix(d_lagrange_dx,d_lagrange_dx,0,coordinate)
```
The second derivative of a function $`f`$ can be found by solving the system
```math
\sum_{j}K_{ij}f_j = \sum_{k}M_{ik}\left[\frac{d^2 f}{d x^2}\right]_k
```
with $`f_j=f(x_j)`$. Here, we have ignored boundary terms
that would be introduced by a nonzero $`df/dx`$ at the boundaries
in $`x`$.

For a system with Jacobian $`x`$ (e.g., cylindrical polar coordinates), we could form the appropriate matrices to evaluate
the 1D, radial Laplacian via
```
coordinate = ElementCoordinates(x,scale,shift)
M = finite_element_matrix(lagrange_x,lagrange_x,1,coordinate)
K = -finite_element_matrix(d_lagrange_dx,d_lagrange_dx,1,coordinate)
```
i.e.,
```math
\sum_{j}K_{ij}f_j = \sum_{k}M_{ik}\left[\frac{1}{x}\frac{d}{d x}\left(x \frac{d f}{d x}\right)\right]_k.
```
Finally, we can form a nonlinear stiffness matrix by inserting an extra enum argument, as follows.
```
N = finite_element_matrix(lagrange_x,lagrange_x,lagrange_x,1,coordinate)
```

To specify an arbitrary kernel function $`\rho(x)`$, the following function call should be used
```
coordinate = ElementCoordinates(x,scale,shift)
M = finite_element_matrix(lagrange_x,lagrange_x,coordinate,
            kernel_function=(v -> rho(v)),additional_quadrature_points=n,
            quadrature_increment=m,
            max_iterations=q, atol=atol, rtol=rtol)
```
where we specify some function `rho(v)` to be the kernel, and we specify that we `n` additional quadrature
points should be used in addition to the number assumed from the size of the reference grid `x`. If `q` is greater than zero, for each iteration we proceed to calculate the matrix again with a further additional `m` quadrature points, until the matrix is calculated to the specified absolute tolerance `atol` and relative tolerance `rtol`, or we call `error()` with an appropriate message for the user. Note that `rho(v)` should be specified in terms of the physical coordinate including the scale and shift factors `v = scale * x + shift`, not in terms of the reference coordinate `x`.

# Matrices in two coordinates (2D)

It is possible to write down weak forms where the multi-dimensional integrals do not separate into a product of one-dimensional integrals, even for a tensor-product basis. This occurs where the kernel is not a separable function, e.g., for the matrix
```math
K_{ikjn} = \int^{x_{2u}}_{x_{2l}}\int^{x_{1u}}_{x_{1l}} P_i(x_1) R_k(x_2) Q_j(x_1) S_n(x_2)\rho(x_1,x_2) d x_1 d x_2,
```
where $`P_i(x_1)`$, $`R_i(x_2)`$, $`Q_i(x_1)`$, $`S_i(x_2)`$ are from the set of Lagrange polynomials basis
functions (and their first derivatives) in $`x_1`$ and $`x_2`$, respectively.

The syntax that we use to create matrices in two coordinates is as follows
```
l_i = lagrange_x
l_j = lagrange_x
l_k = d_lagrange_dx
l_n = d_lagrange_dx
K_2D = finite_element_matrix(l_i,l_j,coordinate_x1,
                            l_k,l_n,coordinate_x2;
                            kernel_function=((x1,x2)->rho(x1,x2)),
                            max_iterations=10, quadrature_increment=10,
                            atol=1.0e-13, rtol=1.0e-13)
```
where we have two coordinates, `x1` (with polynomials `l_i` and `l_j`) and `x2` (with polynomials `l_k` and `l_n`). Note that the `@enum` variables used to select the Lagrange polynomials are the same between the two coordinates.
The matrix created above is indexed in the order `K_2D[i,k,j,n]`, indexing over the first polynomial `l_i` of `x1`, then the first polynomial `l_k` of `x2`, then the second polynomial `l_j` of `x1` with the final index representing the polynomial `l_n` of `x2`. The same adaptive quadrature keyword arguments are used as in the 1D case.
# Examples

There are several examples of taking first and second derivatives
in the tests of this package `FiniteElementMatrices/test/runtests.jl`, as well as more complicated usages in the package
implementing [Fokker-Planck Coulomb collisions](https://github.com/moment-kinetics/FokkerPlanck).