module FiniteElementMatricesTests

using Test: @testset, @test
using FiniteElementMatrices: lagrange_x,
                             d_lagrange_dx,
                             finite_element_matrix,
                             ElementCoordinates
using FastGaussQuadrature: gausslobatto, gaussradau, gausslegendre
using LinearAlgebra: lu, ldiv!, mul!

@enum node_type begin
    GLL
    GLR
    GLe
end

function reference_nodes(nodes::node_type,
                        ngrid::Int64)
    if nodes == GLL
        x, w = gausslobatto(ngrid)
    elseif nodes == GLR
        x, w = gaussradau(ngrid)
    elseif nodes == GLe
        x, w = gausslegendre(ngrid)
    end
    return x
end

function test_polynomial(x;n=2)
    return (x-1)^n
end

function test_polynomial_prime(x;n=2)
    return n*(x-1)^(n-1)
end

function test_first_derivative(;nodes::node_type=GLL,
                            ngrid::Int64=20,
                            func::Function=sin,
                            dfunc::Function=cos,
                            y0::Float64=-1.0,
                            y1::Float64=1.0,
                            pp::Int64=0,
                            atol::Float64=2.0e-13)
    x = reference_nodes(nodes,ngrid)
    scale = 0.5*(y1-y0)
    shift = 0.5*(y0+y1)
    coordinate = ElementCoordinates(x,scale,shift)
    M = finite_element_matrix(lagrange_x,lagrange_x,0,coordinate)
    P = finite_element_matrix(lagrange_x,d_lagrange_dx,pp,coordinate)
    S = finite_element_matrix(d_lagrange_dx,lagrange_x,pp,coordinate)
    # check S_ij = P_ji to cover tests where
    # derivative acts on test function
    @test isapprox(S,transpose(P), atol=2.0e-15)
    result = Array{Float64,1}(undef,ngrid) 
    dummy = Array{Float64,1}(undef,ngrid) 
    err = Array{Float64,1}(undef,ngrid)  
    # function to test
    f = Array{Float64,1}(undef,ngrid)
    df = Array{Float64,1}(undef,ngrid)
    for i in 1:ngrid
        # change to physical coordinate y
        # from reference coordinate x
        y = scale*x[i] + shift
        f[i] = func(y)
        df[i] = (y^pp)*dfunc(y)
    end
    luM = lu(M)
    mul!(dummy,P,f)
    ldiv!(result,luM,dummy)
    @. err = result - df
    maxerr = maximum(abs.(err))
    @test maxerr < atol
    #println(maxerr)
end

function test_first_derivative_nonlinear_operators(
                            ;nodes::node_type=GLL,
                            ngrid::Int64=20,
                            func::Function=sin,
                            dfunc::Function=cos,
                            y0::Float64=-1.0,
                            y1::Float64=1.0,
                            pp::Int64=0,
                            atol::Float64=2.0e-13)
    x = reference_nodes(nodes,ngrid)
    scale = 0.5*(y1-y0)
    shift = 0.5*(y0+y1)
    coordinate = ElementCoordinates(x,scale,shift)
    M = finite_element_matrix(lagrange_x,lagrange_x,0,coordinate)
    Y100 = finite_element_matrix(d_lagrange_dx,lagrange_x,lagrange_x,pp,coordinate)
    Y010 = finite_element_matrix(lagrange_x,d_lagrange_dx,lagrange_x,pp,coordinate)
    Y001 = finite_element_matrix(lagrange_x,lagrange_x,d_lagrange_dx,pp,coordinate)
    # check Y100_ijk = Y010_jik
    @test isapprox(permutedims(Y100, [2,1,3]),Y010,atol=2.0e-15)
    # check Y100_ijk = Y001_kji
    @test isapprox(permutedims(Y100, [3,2,1]),Y001,atol=2.0e-15)
    # check that a first derivative can be carried out with Y100
    result = Array{Float64,1}(undef,ngrid) 
    dummy = Array{Float64,1}(undef,ngrid)
    err = Array{Float64,1}(undef,ngrid)  
    # function to test
    f = Array{Float64,1}(undef,ngrid)
    df = Array{Float64,1}(undef,ngrid)
    for i in 1:ngrid
        y = scale*x[i] + shift
        f[i] = func(y)
        df[i] = (y^pp)*dfunc(y)
    end
    luM = lu(M)
    @. dummy = 0.0
    for k in 1:ngrid
        for j in 1:ngrid
            for i in 1:ngrid
                dummy[k] += Y100[i,j,k]*f[i]
            end
        end
    end
    ldiv!(result,luM,dummy)
    @. err = result - df
    maxerr = maximum(abs.(err))
    @test maxerr < atol
    #println(maxerr)
end

function test_second_derivative(;nodes::node_type=GLL,
                            ngrid::Int64=50,
                            y0::Float64=-1.0,
                            y1::Float64=1.0,
                            atol::Float64=5.0e-10)
    x = reference_nodes(nodes,ngrid)
    scale = 0.5*(y1-y0)
    shift = 0.5*(y0+y1)
    coordinate = ElementCoordinates(x,scale,shift)
    M = finite_element_matrix(lagrange_x,lagrange_x,0,coordinate)
    S = -finite_element_matrix(d_lagrange_dx,lagrange_x,0,coordinate)
    K = -finite_element_matrix(d_lagrange_dx,d_lagrange_dx,0,coordinate)
    result = Array{Float64,1}(undef,ngrid) 
    dummy = Array{Float64,1}(undef,ngrid) 
    err = Array{Float64,1}(undef,ngrid)  
    # function to test
    f = Array{Float64,1}(undef,ngrid)
    df = Array{Float64,1}(undef,ngrid)
    d2f = Array{Float64,1}(undef,ngrid)
    # to test the weak first derivative, 
    # choose a fn that vanishes at x = +-1
    # to avoid including boundary terms
    for i in 1:ngrid
        y = scale*x[i] + shift
        f[i] = sin(pi*(y-shift)/scale)
        df[i] = (pi/scale)*cos(pi*(y-shift)/scale)
    end
    
    # test the performance of a first derivative
    # taken using the "weak" methods via 
    # integration by parts
    luM = lu(M)
    mul!(dummy,S,f)
    ldiv!(result,luM,dummy)
    
    @. err = result - df
    maxerr = maximum(abs.(err))
    @test maxerr < atol
    
    # to test the weak second derivative, 
    # choose a fn that where the derivative vanishes at x = +-1
    # to avoid including boundary terms
    for i in 1:ngrid
        y = scale*x[i] + shift
        f[i] = cos(pi*(y-shift)/scale)
        d2f[i] = -((pi/scale)^2)*cos(pi*(y-shift)/scale)
    end
    # test the performance of a second derivative
    # taken using the "weak" methods via 
    # integration by parts
    mul!(dummy,K,f)
    ldiv!(result,luM,dummy)
    @. err = result - d2f
    maxerr = maximum(abs.(err))
    @test maxerr < atol
    
end

function test_nonlinear_operators(;nodes::node_type=GLL,
                                ngrid::Int64=25,
                                y0::Float64=-1.0,
                                y1::Float64=1.0,
                                atol::Float64=8.0e-13)
    x = reference_nodes(nodes,ngrid)
    scale = 0.5*(y1-y0)
    shift = 0.5*(y0+y1)
    coordinate = ElementCoordinates(x,scale,shift)
    M = finite_element_matrix(lagrange_x,lagrange_x,0,coordinate)
    Y000 = finite_element_matrix(lagrange_x,lagrange_x,lagrange_x,0,coordinate)
    Y100 = finite_element_matrix(d_lagrange_dx,lagrange_x,lagrange_x,0,coordinate)
    Y010 = finite_element_matrix(lagrange_x,d_lagrange_dx,lagrange_x,0,coordinate)
    Y001 = finite_element_matrix(lagrange_x,lagrange_x,d_lagrange_dx,0,coordinate)
    result = Array{Float64,1}(undef,ngrid) 
    dummy = Array{Float64,1}(undef,ngrid)
    err = Array{Float64,1}(undef,ngrid)  
    # function to test
    f = Array{Float64,1}(undef,ngrid)
    u = Array{Float64,1}(undef,ngrid)
    v = Array{Float64,1}(undef,ngrid)
    for i in 1:ngrid
        y = scale*x[i] + shift
        v[i] = cos(pi*(y-shift)/scale)
        u[i] = sin(pi*(y-shift)/scale)
    end
    luM = lu(M)
    # test Y000
    # check that g defined by M*g = Y000*u*v
    # is equal to f = u*v evaluated at collocation points
    for i in 1:ngrid
        f[i] = u[i]*v[i]
    end
    @. dummy = 0.0
    for k in 1:ngrid
        for j in 1:ngrid
            for i in 1:ngrid
                dummy[k] += Y000[i,j,k]*u[i]*v[j]
            end
        end
    end
    ldiv!(result,luM,dummy)
    @. err = result - f
    maxerr = maximum(abs.(err))
    #println("Y000: ",maxerr)
    @test maxerr < atol

    # test Y100
    # check that g defined by M*g = Y100*u*v
    # is equal to f = u'*v evaluated at collocation points
    for i in 1:ngrid
        f[i] = (pi/scale)*v[i]*v[i]
    end
    @. dummy = 0.0
    for k in 1:ngrid
        for j in 1:ngrid
            for i in 1:ngrid
                dummy[k] += Y100[i,j,k]*u[i]*v[j]
            end
        end
    end
    ldiv!(result,luM,dummy)
    @. err = result - f
    maxerr = maximum(abs.(err))
    #println("Y100: ",maxerr)
    @test maxerr < atol

    # test Y010
    # check that g defined by M*g = Y010*u*v
    # is equal to f = u*v' evaluated at collocation points
    for i in 1:ngrid
        f[i] = -(pi/scale)*u[i]*u[i]
    end
    @. dummy = 0.0
    for k in 1:ngrid
        for j in 1:ngrid
            for i in 1:ngrid
                dummy[k] += Y010[i,j,k]*u[i]*v[j]
            end
        end
    end
    ldiv!(result,luM,dummy)
    @. err = result - f
    maxerr = maximum(abs.(err))
    #println("Y010: ",maxerr)
    @test maxerr < atol

    # test Y001
    # check that g defined by M*g = -Y001*u*v
    # is equal to f = (u*v)' evaluated at collocation points
    # we can ignore boundary terms as u = 0 at x = +-1.
    for i in 1:ngrid
        f[i] = (pi/scale)*(v[i]*v[i] - u[i]*u[i])
    end
    @. dummy = 0.0
    for k in 1:ngrid
        for j in 1:ngrid
            for i in 1:ngrid
                dummy[k] -= Y001[i,j,k]*u[i]*v[j]
            end
        end
    end
    ldiv!(result,luM,dummy)
    @. err = result - f
    maxerr = maximum(abs.(err))
    #println("Y001: ",maxerr)
    @test maxerr < atol
end

function runtests()
    @testset "FiniteElementMatrices" begin
        println("test FiniteElementMatrices")
        println("test first derivative:")
        nodes_list = [GLL, GLR, GLe]
        for nodes in nodes_list
            n = 20
            println("    -- test: dfdx $(nodes) ngrid=$(n) trig")
            test_first_derivative(; nodes=nodes, ngrid=n)
            for p in 0:3
                for n in 2+p:10
                    function polynomial(x)
                    return test_polynomial(x,n=n-1-p)
                    end
                    function polynomial_prime(x)
                        return test_polynomial_prime(x,n=n-1-p)
                    end
                    println("    -- test: dfdx $(nodes) ngrid=$(n) poly p=$(p)")
                    test_first_derivative(;
                                        nodes=nodes,
                                        func=polynomial,
                                        dfunc=polynomial_prime,
                                        y0=-0.3, y1=1.0, pp=p,
                                        ngrid=n, atol=2.0e-11)
                    test_first_derivative_nonlinear_operators(
                                ;nodes=nodes,
                                ngrid=n,
                                func=polynomial,
                                dfunc=polynomial_prime,
                                y0=-0.3, y1=1.0, pp=p,
                                atol=2.0e-11)
                end
            end
        end
        println("test second derivative:")
        nodes_list = [GLL, GLR, GLe]
        for nodes in nodes_list
            n = 20
            println("    -- test: $(nodes) ngrid=$(n) trig")
            test_second_derivative(; nodes=nodes,
                                    y0 = -0.5, y1=1.0,
                                    ngrid=n)
        end
        
        println("test nonlinear operators:")
        nodes_list = [GLL, GLR, GLe]
        for nodes in nodes_list
            n = 25
            println("    -- test: $(nodes) ngrid=$(n) trig")
            test_nonlinear_operators(; nodes=nodes,
                                    y0 = -0.6, y1=1.0,
                                    ngrid=n)
        end
    end
end

end

using .FiniteElementMatricesTests
FiniteElementMatricesTests.runtests()