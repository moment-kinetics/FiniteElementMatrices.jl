# FiniteElementMatrices

A package for computing unassembled matrices appearing in $`C^0`$ continuous Galerkin
finite element models.

# Methods

This package computes matrices of the following form.
```math
A_{ij} = \int^{x_u}_{x_l} P_i(x) Q_j(x) x^p d x,
```
where $`P_i`$ and $`Q_i`$ are from the set of Lagrange polynomials basis
functions (and their first derivatives) used for
$`C^0`$ continuous Galerkin
finite element models, $`x_u`$ and $`x_l`$
are the upper and lower limits of the element in the
coordinate $`x`$, and $`p`$ is an integer.

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
A_{ij} = \int^{1}_{-1} P_i(z) Q_j(z) (s z + c)^p d z.
```
We can also compute nonlinear matrices which have the form
```math
N_{ijk} = \int^{x_u}_{x_l} P_i(x) Q_j(x) R_k(x) x^p d x = \int^{1}_{-1} P_i(z) Q_j(z) Q_k(z) (s z + c)^p d z.
```

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

# Examples

There are several examples of taking first and second derivatives
in the tests of this package `FiniteElementMatrices/test/runtests.jl`, as well as more complicated usages in the package
implementing [Fokker-Planck Coulomb collisions](https://github.com/moment-kinetics/FokkerPlanck).