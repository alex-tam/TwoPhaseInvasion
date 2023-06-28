# Control function for 2D Fisher-Stefan level-set solutions
# Nizhum Rahman & Alex Tam, 19/06/2023

# Load packages
using Parameters
using Printf
using Dierckx
using LinearAlgebra
using DifferentialEquations
using Measures
using LaTeXStrings
using DelimitedFiles

# Include external files
include("domain.jl")
include("fkpp.jl")
include("velocity_extension.jl")
include("interface_density.jl")
include("interface_speed.jl")
include("level-set.jl")
include("reinitialisation.jl")

"Data structure for parameters"
@with_kw struct Params
    Du::Float64 = 1.0 # [-] Diffusion coefficient (u)
    λu::Float64 = 1.0 # [-] Reaction rate (u)
    κu::Float64 = 0.2 # [-] Inverse Stefan number (u)
    Dv::Float64 = 1.0 # [-] Diffusion coefficient (v)
    λv::Float64 = 1.0 # [-] Reaction rate (v)
    κv::Float64 = 0.0 # [-] Inverse Stefan number (v)
    αu::Float64 = 0.5 # [-] Maximum initial density (u)
    αv::Float64 = 0.0 # [-] Maximum initial density (v)
    β::Float64 = 1.0 # [-] Initial interface position
    uf::Float64 = 0.0 # [-] Background density at interface
    θb::Float64 = 0.01 # [-] Threshold for whether a grid point is close to interface (relative to Δx)
    θ::Float64 = 1.99 # [-] Parameter for minmod flux-limiter
    Lx::Float64 = 20.0 # [-] Spatial domain limit (x)
    Ly::Float64 = 20.0 # [-] Spatial domain limit (y)
    T::Float64 = 10.0 # [-] End time
    Nx::Int = 201 # [-] Number of grid points (x)
    Ny::Int = 201 # [-] Number of grid points (y)
    Nt::Int = 1001 # [-] Number of time steps
    V_Iterations::Int = 20 # [-] Number of iterations for velocity extrapolation PDE
    ϕ_Iterations::Int = 20 # [-] Number of iterations for reinitialisation PDE
    γ::Float64 = 0.0 # [-] Surface tension coefficient
end

"Obtain initial condition"
function ic(par::Params, x, y)
    U = Array{Float64}(undef, par.Nx, par.Ny) # Pre-allocate 2D array of U
    ϕ = Array{Float64}(undef, par.Nx, par.Ny) # Pre-allocate 2D array of ϕ
    # Loop over grid points
    for i in eachindex(x)
        for j in eachindex(y)
            # Ellipse
            if (1/7*(x[i]-par.Lx/2)^2 + 1/3*(y[j]-par.Ly/2)^2)< par.β
                U[i,j] = par.αu # Constant density (inside Ω1)
            elseif (1/7*(x[i]-par.Lx/2)^2 + 1/3*(y[j]-par.Ly/2)^2) == par.β
                U[i,j] = par.uf # Constant density (on interface)
            else
                U[i,j] = par.αv # Initial density (outside Ω1)
            end
            ϕ[i,j] = (1/7*(x[i]-par.Lx/2)^2 + 1/3*(y[j]-par.Ly/2)^2)- par.β # Initial signed-distance
            
        end
    end
    return U, ϕ
end

"Build vector from matrix, ordered by entries in D"
function build_vector(U::Array{Float64}, D)
    u = Vector{Float64}() # Pre-allocate empty vector
    for gp in D
        push!(u, U[gp.xInd, gp.yInd])
    end
    return u
end

"Build matrix from vector ordered by entries in D"
function build_u_matrix(u::Vector, y, par, D)
    U = zeros(par.Nx, par.Ny) # Pre-allocate (incorporate a Dirichlet condition on right boundary)
    for i in eachindex(D)
        U[D[i].xInd, D[i].yInd] = u[i]
    end
    # Apply zero derivative conditions
    U[par.Nx,:] .= U[par.Nx-1,:]
    U[1,:] .= U[2,:]
    U[:,par.Ny] .= U[:,par.Ny-1,:]
    U[:,1] .= U[:,2]
    return U
end

"Build matrix from vector ordered by entries in D"
function build_v_matrix(v::Vector, par, D)
    V = zeros(par.Nx, par.Ny) # Pre-allocate (incorporate a Dirichlet condition on computational boundary)
    for i in eachindex(D)
        V[D[i].xInd, D[i].yInd] = v[i]
    end
    return V
end

"Compute a solution"
function fisher_stefan_2d()
    # Parameters and domain
    par = Params() # Initialise data structure of model parameters
    nx::Int = (par.Nx-1)/2; ny::Int = (par.Ny-1)/2 # Indices for slice plots
    x = range(0, par.Lx, length = par.Nx); dx = x[2] - x[1] # Computational domain (x)
    y = range(0, par.Ly, length = par.Ny); dy = y[2] - y[1] # Computational domain (y)
    t = range(0, par.T, length = par.Nt); dt = t[2] - t[1] # Time domain
    writedlm("x.csv", x); writedlm("y.csv", y); writedlm("t.csv", t) # Write data to files
    # Initial condition
    U, ϕ = ic(par, x, y) # Obtain initial density and ϕ
    writedlm("U-0.csv", U); writedlm("Phi-0.csv", ϕ) # Write data to files
    plot_times = Vector{Int}() # Vector of time-steps at which data is obtained
    writedlm("plot_times.csv", plot_times)
    # Time stepping
    for i = 1:par.Nt-1
        # 1. Find Ω, dΩ, and irregular grid points
        D = find_domain(par, ϕ)
        dΩ = find_interface(par, D, ϕ)
        # 2. Solve FKPP equation on Ω
        uf = interface_density(dΩ, ϕ, par, dx, dy) # Density on interface for BC
        @time U = fkpp(D, dΩ, U, ϕ, uf, y, par, dx, dy, dt, i)
        # 3. Compute extension velocity field
        V = extend_velocity(D, dΩ, U, ϕ, par, dx, dy)
        # 4. Solve level-set equation
        ϕ = level_set(V, ϕ, par, dx, dy, dt)
        # 5. Re-initialise level-set function as a signed-distance function
        if mod(i, 1) == 0
            ϕ = reinitialisation(ϕ, par, dx, dy, par.ϕ_Iterations)
        end
        # Optional: Post-processing
        if mod(i, 100) == 0
            writedlm("ux-$i.csv", U[:,ny])
            writedlm("uy-$i.csv", U[nx,:])
            writedlm("U-$i.csv", U)
            writedlm("V-$i.csv", V)
            writedlm("Phi-$i.csv", ϕ)
            push!(plot_times, i)
            writedlm("plot_times.csv", plot_times)
        end
    end
end

@time fisher_stefan_2d()