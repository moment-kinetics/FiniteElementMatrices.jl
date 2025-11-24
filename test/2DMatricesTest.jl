using FiniteElementMatrices

using FastGaussQuadrature: gausslobatto, gaussradau, gausslegendre
using Test

# x1
# endpoints of grid
y1 = 1.2
y0 = -1.3
ngrid_x1 = 5
x, w = gausslegendre(ngrid_x1)
scale_x1 = 0.5*(y1-y0)
shift_x1 = 0.5*(y0+y1)
coordinate_x1 = ElementCoordinates(x,scale_x1,shift_x1)
# x2
# endpoints of grid
y1 = 0.9
y0 = -0.78
ngrid_x2 = 4
x, w = gausslegendre(ngrid_x2)
scale_x2 = 0.5*(y1-y0)
shift_x2 = 0.5*(y0+y1)
coordinate_x2 = ElementCoordinates(x,scale_x2,shift_x2)

M_x1 = zeros(ngrid_x1,ngrid_x1)
M_x2 = zeros(ngrid_x2,ngrid_x2)
M_2D_err = zeros(ngrid_x1,ngrid_x2,ngrid_x1,ngrid_x2)
M_2D_exact = zeros(ngrid_x1,ngrid_x2,ngrid_x1,ngrid_x2)


@testset "FiniteElementMatrices2Ddev" begin
    println("test 2D matrices:")
    for power1 in 1:3
        for power2 in 1:3
            for l1A in (lagrange_x, d_lagrange_dx)
                for l1B in (lagrange_x, d_lagrange_dx)
                    for l2A in (lagrange_x, d_lagrange_dx)
                        for l2B in (lagrange_x, d_lagrange_dx)
                            #1D matrices
                            M_x1 .= finite_element_matrix(l1A,l1B,power1,coordinate_x1)
                            M_x2 .= finite_element_matrix(l2A,l2B,power2,coordinate_x2)
                            #2D matrices
                            M_2D .= finite_element_matrix(l1A,l1B,power1,coordinate_x1,
                                                        l2A,l2B,power2,coordinate_x2)
                            for i2 in 1:ngrid_x2
                                for i1 in 1:ngrid_x1
                                    for j2 in 1:ngrid_x2
                                        for j1 in 1:ngrid_x1
                                            M_2D_exact[j1,j2,i1,i2] = M_x1[j1,i1]*M_x2[j2,i2]
                                        end
                                    end
                                end
                            end
                            @. M_2D_err = abs(M_2D-M_2D_exact)
                            @test maximum(M_2D_err) < 1.0e-13
                            #println(maximum(M_2D_err))
                            #println(maximum(M_2D))
                        end
                    end
                end
            end
        end
    end
end # FiniteElementMatrices2Ddev