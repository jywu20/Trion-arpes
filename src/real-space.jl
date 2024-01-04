using StaticArrays
using LinearAlgebra

struct UniformGrid2D 
    Lx::Float64
    Ly::Float64
    dx::Float64
    dy::Float64
    points::Vector{SVector{2, Float64}}
end

"""
A two-dimensional grid.
Lx and Ly are the size of the grid in the x and y directions
(and therefore x_min ≈ -L/2, etc.)
and Nx and Ny are the number of samples in the two directions.

It's a good idea to set Nx and Ny to even numbers 
if the (0, 0)  point is to be included into the grid.
"""
function UniformGrid2D(Lx, Ly, Nx, Ny)
    dx = Lx / Nx 
    dy = Ly / Ny
    x_start = - Lx/2 + dx 
    x_end   = Lx / 2
    y_start = - Ly/2 + dy
    y_end   = Ly / 2

    points = reshape(
        Iterators.product(x_start:dx:x_end, y_start:dy:y_end) |> collect |> permutedims, 
        Nx * Ny
    )
    
    UniformGrid2D(Lx, Ly, dx, dy, points)
end

ϕ_r_1s_def(r::SVector, a_ex::Float64) = sqrt(2 / (π * a_ex^2)) * exp(- norm(r) / a_ex)

function overlap_numerical(ψ_1, ψ_2, grid::UniformGrid2D)
    dA = grid.dx * grid.dy
    dA * sum(grid.points) do pt
        ψ_1(pt)' * ψ_2(pt) 
    end
end

function coulomb_2D_numerical(ψ_1, ψ_2, grid::UniformGrid2D)
    dA = grid.dx * grid.dy   
    dA * sum(grid.points) do pt
        r = norm(pt)
        if r == 0.0
            return 0.0 
        end
        ψ_1(pt)' * ψ_2(pt) * 1 / r 
    end
end

function coulomb_twopoint_2D_numerical(ψ_1, ψ_2, grid::UniformGrid2D)
    dA = grid.dx * grid.dy
    dA^2 * sum(Iterators.product(grid.points, grid.points)) do (pt1, pt2) 
        r12 = norm(pt1 - pt2)
        if r12 == 0.0
            return 0.0  
        end
        ψ_1(pt1)' * ψ_2(pt2) * 1 / r12
    end
end