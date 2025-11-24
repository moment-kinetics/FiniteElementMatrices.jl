"""
"""
module FiniteElementMatrices

using LagrangePolynomials: lagrange_poly,
                           lagrange_poly_derivative,
                           LagrangePolyData
using FastGaussQuadrature: gausslegendre

export lagrange_x,
       d_lagrange_dx,
       finite_element_matrix,
       element_coordinates

@enum lagrange_function_type begin
    lagrange_x
    d_lagrange_dx
end

struct element_coordinates
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
    function element_coordinates(x_nodes::AbstractArray{Float64,1},
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

function select_lagrange_function(fn_type::lagrange_function_type)
    if fn_type == lagrange_x
        func = lagrange_poly
    elseif fn_type == d_lagrange_dx
        func = lagrange_poly_derivative
    end
    return func
end

function select_lagrange_prefactor(fn_type::lagrange_function_type,
                                    scale::Float64)
    if fn_type == lagrange_x
        prefactor = 1.0
    elseif fn_type == d_lagrange_dx
        prefactor = 1.0/scale
    end
    return prefactor
end

function get_polyz(power::Int64,
                    scale::Float64,
                    shift::Float64,
                    zz::Float64)
    polyz = 0.0
    for q in 0:power
        polyz += (binomial(power,q)*
                    ((scale*zz)^q)*
                    (shift^(power-q)))
    end
    return polyz
end

function finite_element_matrix(
    fn1_type::lagrange_function_type,
    fn2_type::lagrange_function_type,
    power::Int64,
    coordinate::element_coordinates
    )
    lpoly_data = coordinate.lpoly_data
    ngrid = length(coordinate.lpoly_data.x_nodes)
    scale = coordinate.scale
    shift = coordinate.shift
    # the finite element array to be returned
    matrix = zeros(Float64,ngrid,ngrid)
    # the function objects for the required polynomials
    lagrange1 = select_lagrange_function(fn1_type)
    lagrange2 = select_lagrange_function(fn2_type)
    # the normalisation factors due to v = scale * x + shift
    prefactor1 = select_lagrange_prefactor(fn1_type,scale)
    prefactor2 = select_lagrange_prefactor(fn2_type,scale)
    # nquad chosen for exact results for all power, ngrid
    nquad = ngrid + power
    zz, wz = gausslegendre(nquad)
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
                polyz = get_polyz(power,scale,shift,zzl)
                matrix[i,j] += (scale*wz[l]*polyz*
                           prefactor1*lagrange1(ith_lpoly_data,zzl)*
                           prefactor2*lagrange2(jth_lpoly_data,zzl))
            end
        end
    end
    return matrix
end

function finite_element_matrix(
    fn1_type::lagrange_function_type,
    fn2_type::lagrange_function_type,
    fn3_type::lagrange_function_type,
    power::Int64,
    coordinate::element_coordinates
    )
    lpoly_data = coordinate.lpoly_data
    ngrid = length(coordinate.lpoly_data.x_nodes)
    scale = coordinate.scale
    shift = coordinate.shift
    # the finite element array to be returned
    matrix = zeros(Float64,ngrid,ngrid,ngrid)
    # the function objects for the required polynomials
    lagrange1 = select_lagrange_function(fn1_type)
    lagrange2 = select_lagrange_function(fn2_type)
    lagrange3 = select_lagrange_function(fn3_type)
    # the normalisation factors due to v = scale * x + shift
    prefactor1 = select_lagrange_prefactor(fn1_type,scale)
    prefactor2 = select_lagrange_prefactor(fn2_type,scale)
    prefactor3 = select_lagrange_prefactor(fn3_type,scale)
    # nquad chosen for exact results for all power, ngrid
    nquad = 2*ngrid + power
    zz, wz = gausslegendre(nquad)
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
                    polyz = get_polyz(power,scale,shift,zzl)
                    matrix[i,j,k] += (scale*wz[l]*polyz*
                            prefactor1*lagrange1(ith_lpoly_data,zzl)*
                            prefactor2*lagrange2(jth_lpoly_data,zzl)*
                            prefactor3*lagrange3(kth_lpoly_data,zzl))
                end
            end
        end
    end
    return matrix
end

end