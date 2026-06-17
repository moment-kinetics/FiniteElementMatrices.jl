"""
"""
module FiniteElementMatrices

using LagrangePolynomials: lagrange_poly,
                           lagrange_poly_derivative,
                           LagrangePolyData,
                           LagrangePolyDataBase
using FastGaussQuadrature: gausslegendre

export lagrange_x,
       d_lagrange_dx,
       finite_element_matrix,
       ElementCoordinates,
       GLSpecifiedLimits

@enum LagrangeFunctionType begin
    lagrange_x
    d_lagrange_dx
end

abstract type Abstract1DQuadrature end
struct Default1DQuadrature <: Abstract1DQuadrature end
struct GLSpecifiedLimits <: Abstract1DQuadrature
    # minimum of range of integration
    # in physical grid v = s x + c
    # with x in [-1, 1]
    v_min::Float64
    # maximum of range of integration
    v_max::Float64
end

struct ElementCoordinates
    # precomputed data for calculating
    # the Lagrange polynomials, including
    # the reference x nodes
    lpoly_data::LagrangePolyData
    # scale factor s for transform
    # v = s x + c
    # from reference grid x on [-1,1] (or (-1,1], [-1,1), (-1,1))
    # to physical grid v
    scale::Float64
    # shift factor c
    shift::Float64
    """
    `x_nodes`: reference grid x on [-1,1]
            (or (-1,1], [-1,1), (-1,1), if endpoints
            are not included)
    `scale` : normalisation factor, see below.
    `shift` : normalisation factor, such that the physical coordinate
            `v = `scale * x_nodes + shift` in this element.
    """
    function ElementCoordinates(x_nodes::AbstractArray{Float64,1},
                                    scale::Float64,
                                    shift::Float64)
        # initialise and save Lagrange polynomial
        # data for this 1D element.
        lpoly_data = LagrangePolyData(x_nodes)
        return new(lpoly_data,
                    scale,
                    shift)
    end
end

function select_lagrange_function(fn_type::LagrangeFunctionType,scale::Float64)
    # including the normalisation factors due to v = scale * x + shift
    if fn_type == lagrange_x
        func = ((lpoly_data::LagrangePolyDataBase,z::Float64) -> lagrange_poly(lpoly_data,z))
    elseif fn_type == d_lagrange_dx
        # 1/scale due to scale in x in d l(z(x)) /d x
        func = ((lpoly_data::LagrangePolyDataBase,z::Float64) -> (1.0/scale)*lagrange_poly_derivative(lpoly_data,z))
    end
    return func
end

function check_power(power::Int64, varname::String)
    if power < 0
        error("Require power > 0, whereas $varname=$power was passed")
    end
    return nothing
end
# Method for p-adaptive quadrature in 1D and 2D.
# Note that the the rate at which the number of quadrature
# points is increased is the same for each coordinate.
function finite_element_matrix(args...;
    atol::Float64=1.0e-12, rtol::Float64=1.0e-12, max_iterations::Int64=1,
    quadrature_increment::Int64=5, verbose::Bool=false, kwargs...)
    matrix = _finite_element_matrix(args...;
                adaptive_quadrature_points=0,
                kwargs...)
    iteration = 0
    maxvalue = maximum(abs.(matrix))
    maxerr = maxvalue
    tolerance_exceeded = max_iterations > 0
    while iteration < max_iterations && tolerance_exceeded
        # increment the iteration
        iteration += 1
        matrixnew = _finite_element_matrix(args...;
                        adaptive_quadrature_points=quadrature_increment*iteration,
                        kwargs...)
        maxvalue = maximum(abs.(matrixnew))
        maxerr = maximum(abs.(matrixnew .- matrix))
        tolerance_exceeded = maxerr > atol + rtol*maxvalue
        # update the matrix
        matrix .= matrixnew
    end
    if verbose && !(tolerance_exceeded)
        if iteration > 0
            println("finite_element_matrix() converged after $iteration iterations \n with maxerr = $maxerr, adaptive_quadrature_points=$(quadrature_increment*iteration) \n")
        else
            println("finite_element_matrix() returned a matrix without p-adaptive quadrature refinement.")
        end
    end
    if tolerance_exceeded
        error("Unable to achieve specified error bounds: \n
              maxerr < atol + rtol*maxvalue \n
              after max_iterations = $max_iterations with atol=$atol and rtol=$rtol we obtain \n
              maxerr = $maxerr > $atol + $rtol*$maxvalue")
    end
    return matrix
end

# Method for computing polynomial integrals exactly
# for matrices with two Lagrange polymomials in 1D
function finite_element_matrix(
    fn1_type::LagrangeFunctionType,
    fn2_type::LagrangeFunctionType,
    power::Int64,
    coordinate::ElementCoordinates
    )
    check_power(power, "power")
    return _finite_element_matrix(fn1_type, fn2_type, coordinate;
                kernel_function=((v -> v^power)),
                additional_quadrature_points=power)
end

# The basic method for matrices with two Lagrange polymomials in 1D
function _finite_element_matrix(
    fn1_type::LagrangeFunctionType,
    fn2_type::LagrangeFunctionType,
    coordinate::ElementCoordinates;
    # function of the "physical" coord v = s z + c
    # rather than of the reference coordinate z on [-1,1]
    kernel_function::TFunction=((v -> 1.0)),
    additional_quadrature_points::Int64=0,
    # argument to permit single implementation of
    adaptive_quadrature_points::Int64=0,
    quadrature_option::TQuad=Default1DQuadrature()
    ) where {TFunction, TQuad <: Abstract1DQuadrature}
    lpoly_data = coordinate.lpoly_data
    ngrid = length(coordinate.lpoly_data.x_nodes)
    scale = coordinate.scale
    shift = coordinate.shift
    # the finite element array to be returned
    matrix = zeros(Float64,ngrid,ngrid)
    # the function objects for the required polynomials
    lagrange1 = select_lagrange_function(fn1_type,scale)
    lagrange2 = select_lagrange_function(fn2_type,scale)
    # nquad chosen for exact results for default inputs
    # with kernel = 1.0 and zero additional quadrature points
    nquad = ngrid + additional_quadrature_points + adaptive_quadrature_points
    zz, wz = quadrature1D(quadrature_option,nquad,coordinate)
    # compute integral
    # int P_i(z) Q_j(z) poly(z) s d z
    # with poly(z) = (s z + c)^power
    # s = scale, c = shift and
    # P_i(z), Q_i(z) in [l_i(z), (1/s) d l_i(z) / d z]
    for j in 1:ngrid
        jth_lpoly_data = lpoly_data.lpoly_data[j]
        for i in 1:ngrid
            ith_lpoly_data = lpoly_data.lpoly_data[i]
            for l in 1:nquad
                zzl = zz[l]
                funcz = kernel_function(scale*zzl+shift)::Float64
                matrix[i,j] += (scale*wz[l]*funcz*
                           lagrange1(ith_lpoly_data,zzl)*
                           lagrange2(jth_lpoly_data,zzl))
            end
        end
    end
    return matrix
end

# Method for computing polynomial integrals exactly
# for matrices with three Lagrange polymomials in 1D
function finite_element_matrix(
    fn1_type::LagrangeFunctionType,
    fn2_type::LagrangeFunctionType,
    fn3_type::LagrangeFunctionType,
    power::Int64,
    coordinate::ElementCoordinates
    )
    check_power(power, "power")
    return _finite_element_matrix(fn1_type, fn2_type, fn3_type,
                coordinate; kernel_function=((v -> v^power)),
                additional_quadrature_points=power)
end

# The basic method for matrices with three Lagrange polymomials in 2D
function _finite_element_matrix(
    fn1_type::LagrangeFunctionType,
    fn2_type::LagrangeFunctionType,
    fn3_type::LagrangeFunctionType,
    coordinate::ElementCoordinates;
    # function of the "physical" coord v = s z + c
    # rather than of the reference coordinate z on [-1,1]
    kernel_function::TFunction=((v -> 1.0)),
    additional_quadrature_points::Int64=0,
    adaptive_quadrature_points::Int64=0,
    quadrature_option::TQuad=Default1DQuadrature(),
    ) where {TFunction, TQuad <: Abstract1DQuadrature}
    lpoly_data = coordinate.lpoly_data
    ngrid = length(coordinate.lpoly_data.x_nodes)
    scale = coordinate.scale
    shift = coordinate.shift
    # the finite element array to be returned
    matrix = zeros(Float64,ngrid,ngrid,ngrid)
    # the function objects for the required polynomials
    lagrange1 = select_lagrange_function(fn1_type,scale)
    lagrange2 = select_lagrange_function(fn2_type,scale)
    lagrange3 = select_lagrange_function(fn3_type,scale)
    # nquad chosen for exact results for default inputs
    # with kernel = 1.0 and zero additional quadrature points
    nquad = 2*ngrid + additional_quadrature_points + adaptive_quadrature_points
    zz, wz = quadrature1D(quadrature_option,nquad,coordinate)
    # compute integral
    # int P_i(z) Q_j(z) S_k(z) poly(z) d z
    # with poly(z) = (s z + c)^power
    # s = scale, c = shift and
    # P_i(z), Q_i(z), S_i(z) in [l_i(z), (1/s) d l_i(z) / d z]
    for k in 1:ngrid
        kth_lpoly_data = lpoly_data.lpoly_data[k]
        for j in 1:ngrid
            jth_lpoly_data = lpoly_data.lpoly_data[j]
            for i in 1:ngrid
                ith_lpoly_data = lpoly_data.lpoly_data[i]
                for l in 1:nquad
                    zzl = zz[l]
                    funcz = kernel_function(scale*zzl+shift)::Float64
                    matrix[i,j,k] += (scale*wz[l]*funcz*
                            lagrange1(ith_lpoly_data,zzl)*
                            lagrange2(jth_lpoly_data,zzl)*
                            lagrange3(kth_lpoly_data,zzl))
                end
            end
        end
    end
    return matrix
end

# Method for computing polynomial integrals exactly
# for matrices with two 2D basis functions
# (two Lagrange polymomials per coordinate) in 2D
function finite_element_matrix(
    fn1_x1_type::LagrangeFunctionType,
    fn2_x1_type::LagrangeFunctionType,
    power_x1::Int64,
    coordinate_x1::ElementCoordinates,
    fn1_x2_type::LagrangeFunctionType,
    fn2_x2_type::LagrangeFunctionType,
    power_x2::Int64,
    coordinate_x2::ElementCoordinates
    )
    check_power(power_x1, "power_x1")
    check_power(power_x2, "power_x2")
    return _finite_element_matrix(fn1_x1_type, fn2_x1_type, coordinate_x1,
                fn1_x2_type, fn2_x2_type, coordinate_x2;
                kernel_function=((v1,v2) -> (v1^power_x1)*(v2^power_x2)),
                additional_quadrature_points_x1=power_x1,
                additional_quadrature_points_x2=power_x2)
end

# Basic method for computing matrices with two 2D basis functions
# (two Lagrange polymomials per coordinate) in 2D
function _finite_element_matrix(
    fn1_x1_type::LagrangeFunctionType,
    fn2_x1_type::LagrangeFunctionType,
    coordinate_x1::ElementCoordinates,
    fn1_x2_type::LagrangeFunctionType,
    fn2_x2_type::LagrangeFunctionType,
    coordinate_x2::ElementCoordinates;
    # function of the "physical" coord v = s z + c
    # rather than of the reference coordinate z on [-1,1]
    kernel_function::TFunction=((v1,v2) -> 1.0),
    additional_quadrature_points_x1::Int64=0,
    additional_quadrature_points_x2::Int64=0,
    adaptive_quadrature_points::Int64=0,
    quadrature_option_x1::TQuad1=Default1DQuadrature(),
    quadrature_option_x2::TQuad2=Default1DQuadrature(),
    ) where {TFunction,
        TQuad1 <: Abstract1DQuadrature,
        TQuad2 <: Abstract1DQuadrature}
    # coordinate x1 data
    lpoly_data_x1 = coordinate_x1.lpoly_data
    ngrid_x1 = length(coordinate_x1.lpoly_data.x_nodes)
    scale_x1 = coordinate_x1.scale
    shift_x1 = coordinate_x1.shift
    # coordinate x2 data
    lpoly_data_x2 = coordinate_x2.lpoly_data
    ngrid_x2 = length(coordinate_x2.lpoly_data.x_nodes)
    scale_x2 = coordinate_x2.scale
    shift_x2 = coordinate_x2.shift
    # the finite element array to be returned
    matrix = zeros(Float64,ngrid_x1,ngrid_x2,ngrid_x1,ngrid_x2)
    # the function objects for the required polynomials
    lagrange11 = select_lagrange_function(fn1_x1_type,scale_x1)
    lagrange12 = select_lagrange_function(fn1_x2_type,scale_x2)
    lagrange21 = select_lagrange_function(fn2_x1_type,scale_x1)
    lagrange22 = select_lagrange_function(fn2_x2_type,scale_x2)
    # the normalisation factors due to v = scale * x + shift
    # nquad chosen for exact results for default inputs
    # with kernel = 1.0 and zero additional quadrature points
    nquad_x1 = ngrid_x1 + additional_quadrature_points_x1 + adaptive_quadrature_points
    zz_x1, wz_x1 = quadrature1D(quadrature_option_x1,nquad_x1,coordinate_x1)
    nquad_x2 = ngrid_x2 + additional_quadrature_points_x2 + adaptive_quadrature_points
    zz_x2, wz_x2 = quadrature1D(quadrature_option_x2,nquad_x2,coordinate_x2)
    # compute integral
    # \int \int \left(P1_i(z_1) Q1_j(z_1) P2_i(z_1) Q2_j(z_2)
    # kernel(s_1 z_1 + c_1, s_2 z_2 + c_2) \right) s_1 s_2 d z_1 d z_2
    # s = scale, c = shift and
    # P_i(z), Q_i(z) in [l_i(z), (1/s) d l_i(z) / d z]
    for jx2 in 1:ngrid_x2
        jx2_lpoly_data = lpoly_data_x2.lpoly_data[jx2]
        for jx1 in 1:ngrid_x1
            jx1_lpoly_data = lpoly_data_x1.lpoly_data[jx1]
            for ix2 in 1:ngrid_x2
                ix2_lpoly_data = lpoly_data_x2.lpoly_data[ix2]
                for ix1 in 1:ngrid_x1
                    ix1_lpoly_data = lpoly_data_x1.lpoly_data[ix1]
                    for lx2 in 1:nquad_x2
                        zzl2 = zz_x2[lx2]
                        wgt2 = scale_x2*wz_x2[lx2]
                        lagrange_factor2 = (lagrange12(ix2_lpoly_data,zzl2)*
                                               lagrange22(jx2_lpoly_data,zzl2))
                        for lx1 in 1:nquad_x1
                            zzl1 = zz_x1[lx1]
                            wgt1 = scale_x1*wz_x1[lx1]
                            lagrange_factor1 = (lagrange11(ix1_lpoly_data,zzl1)*
                                               lagrange21(jx1_lpoly_data,zzl1))
                            kernel_x1_x2 = kernel_function(scale_x1*zzl1+shift_x1,
                                                        scale_x2*zzl2+shift_x2)::Float64
                            matrix[ix1,ix2,jx1,jx2] += (wgt1*wgt2*kernel_x1_x2*
                                                        lagrange_factor1*
                                                        lagrange_factor2)
                        end
                    end
                end
            end
        end
    end
    return matrix
end

function quadrature1D(::Default1DQuadrature, nquad::Int64, coordinate::ElementCoordinates)
    # default quadrature running from [-1,1]
    zz, wz = gausslegendre(nquad)
    return zz, wz
end
function quadrature1D(qopt::GLSpecifiedLimits, nquad::Int64, coordinate::ElementCoordinates)
    # limits specified for integration
    v_min = qopt.v_min
    v_max = qopt.v_max
    # limits based on the range supported by the Lagrange polynomials
    v_lower_limit = coordinate.shift - coordinate.scale
    v_upper_limit = coordinate.shift + coordinate.scale
    # some checks on the inputs
    if !(v_min < v_upper_limit)
        error("invalid integration range: v_min >= v_upper_limit")
    elseif !(v_max > v_lower_limit)
        error("invalid integration range: v_max =< v_lower_limit")
    elseif !(v_min < v_max)
        error("invalid integration range: v_min >= v_max")
    end
    # obtain valid limits of integration
    # in the physical range [v_lower_limit,v_upper_limit]
    v_min = max(v_min, v_lower_limit)
    v_max = min(v_max, v_upper_limit)
    # in the reference range [-1,1]
    x_min = (v_min - coordinate.shift)/coordinate.scale
    x_max = (v_max - coordinate.shift)/coordinate.scale
    # default quadrature running from [-1,1]
    zz, wz = gausslegendre(nquad)
    # quadrature running from [x_min, x_max]
    zx = zeros(Float64, nquad)
    wx = zeros(Float64, nquad)
    x_scale = 0.5*(x_max - x_min)
    x_shift = 0.5*(x_max + x_min)
    @. zx = x_scale*zz + x_shift
    @. wx = x_scale*wz
    return zx, wx
end


end