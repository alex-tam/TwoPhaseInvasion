# Linear stability analysis for the two-phase Fisher-Stefan model
# Alex Tam, 13/01/2023

# Load required packages
using Parameters
using LinearAlgebra
using Printf
using Plots
using Measures
using LaTeXStrings
using DelimitedFiles

##### General #####
"Data structure for model and numerical parameters"
@with_kw struct Params
    L::Float64 = 20.0 # [-] Domain half-width
    Nz::Int = 1001 # [-] Number of grid points
    uf::Float64 = 0.0 # [-] # Front density (u)
    vf::Float64 = 0.0 # [-] # Front density (v)
    D::Float64 = 1.0 # [-] Scaled diffusivity
    λ::Float64 = 1.0 # [-] Scaled reaction
    γ::Float64 = 0.0 # [-] Surface tension coefficient
    κu::Float64 = 0.2 # [-]
    κv::Float64 = 0.1 # [-]
end

"Main subroutine for linear stability analysis"
function main()
    gr()
    default(fontfamily = "Computer Modern", titlefontsize = 18, guidefontsize = 18, tickfontsize = 14, legendfontsize = 12)
    # 1. Parameters and domain
    par = Params() # Initialise structure for parameters
    zl = range(-par.L, 0.0, length = par.Nz)
    zr = range(0.0, par.L, length = par.Nz)
    dz = zl[2] - zl[1]
    Q = range(0.0, 5.0, length = 51) # Pre-allocate wave numbers
    Ω = Vector{Float64}() # Pre-allocate growth rates
    # 2. Obtain travelling wave solution using shooting method
    c, u0, v0 = tws(par, dz, 1.0, 1e-6)
    @printf("Wave speed: c = %f.\n", c)
    plot(zl, u0, linewidth = 2, xlabel = L"$\xi$", label = L"$u_0(\xi)$", legend=:bottomleft)
    plot!(zr, v0, linewidth = 2, label = L"$v_0(\xi)$")
    D = par.D; l = par.λ; ku = par.κu; kv = par.κv; uf = par.uf; vf = par.vf; g = par.γ
    savefig("TWS-($D,$l,$ku,$kv,$uf,$vf).png")
    # 3. Obtain dispersion relation using shooting method
    for q in Q
        ω, u1, v1 = lsa(par, c, u0, v0, q, dz, 1.0, 1e-6)
        push!(Ω, ω)
        @printf("Growth rate: ω = %f.\n", ω)
        plot(zl, u1, linewidth = 2, xlabel = L"$\xi$", label = L"$u_1(\xi)$", legend=:bottomleft)
        plot!(zr, v1, linewidth = 2, label = L"$v_1(\xi)$")
        savefig("Fisher-Stefan_TwoPhase_Corrections.png")
    end
    # 4. Plot dispersion relation
    writedlm("q-($D,$l,$ku,$kv,$uf,$vf,$g).csv", Q)
    writedlm("omega-($D,$l,$ku,$kv,$uf,$vf,$g).csv", Ω)
    plot(Q, Ω, linewidth=2, xlabel = L"$q$", ylabel = L"$\omega$", legend=false)
    savefig("q-($D,$l,$ku,$kv,$uf,$vf,$g).png")
end

"Finite difference approximation for du/dz at z = 0"
function dudz(u, dz)
    return (3*u[end] - 4*u[end-1] + u[end-2])/(2*dz) # Second-order
end

"Finite difference approximation for du/dz at z = 0"
function dvdz(v, dz)
    return (-3*v[1] + 4*v[2] - v[3])/(2*dz) # Second-order
end

##### Travelling Wave Solutions #####

"Implement shooting method for travelling wave solutions"
function tws(par, dz, c, ϵ)
    # Implement Newton's method
    for j = 1:20
        c_old = c # Store previous guess
        # Solve boundary-value problems
        u = bvp_u0(par, dz, c)
        v = bvp_v0(par, dz, c)
        up = bvp_u0(par, dz, c+ϵ)
        vp = bvp_v0(par, dz, c+ϵ)
        # Implement Newton's method
        f = c + par.κu*dudz(u,dz) + par.κv*dvdz(v,dz) # Interface condition
        fp = c + ϵ + par.κu*dudz(up,dz) + par.κv*dvdz(vp,dz) # Perturbed interface condition
        d = (fp-f)/ϵ # Approximate df/dϵ
        c = c - f/d # Newton iteration
        if abs(c - c_old) < 1e-6
            return c, bvp_u0(par, dz, c), bvp_v0(par, dz, c)
        end
    end
    @printf("Shooting method for TWS did not converge.\n")
end

"Nonlinear boundary-value problem for u0"
function bvp_u0(par, dz, c)
    u = ones(par.Nz) # Initial guess
    # Newton's method
    for i = 1:20
        u_old = u # Store previous iteration u
        F = Fu0(par, dz, u, c) # Construct F
        J = Ju0(par, dz, u, c) # Construct Jacobian
        u = u - J\F # Newton iteration
        if norm(u - u_old) < 1e-6
            return u
        end
    end
    @printf("TWS BVP for u did not converge.\n")
end

"Nonlinear boundary value problem for v0"
function bvp_v0(par, dz, c)
    v = ones(par.Nz) # Initial guess
    # Newton's method
    for i = 1:20
        v_old = v # Store previous iteration u
        F = Fv0(par, dz, v, c) # Construct F
        J = Jv0(par, dz, v, c) # Construct Jacobian
        v = v - J\F # Newton iteration
        if norm(v - v_old) < 1e-6
            return v
        end
    end
    @printf("TWS BVP for v did not converge.\n")
end

"Vector function for u"
function Fu0(par, dz, u, c)
    F = Vector{Float64}(undef, par.Nz) # Pre-allocate F
    F[1] = u[1] - 1.0 # Dirichlet condition on left boundary
    for i = 2:par.Nz-1
        F[i] = (u[i+1]-2*u[i]+u[i-1])/(dz^2) + c*(u[i+1]-u[i-1])/(2*dz) + u[i]*(1-u[i])
    end
    F[par.Nz] = u[par.Nz] - par.uf # Dirichlet condition at interface
    return F
end

"Jacobian matrix for u"
function Ju0(par, dz, u, c)
    J = zeros(par.Nz, par.Nz) # Pre-allocate Jacobian
    J[1,1] = 1.0 # Jacobian entry for left Dirichlet condition
    J[par.Nz, par.Nz] = 1.0 # Jacobian entry for right Dirichlet condition
    for i = 2:par.Nz-1
        J[i,i-1] = 1/(dz^2) - c/(2*dz)
        J[i,i] = -2/(dz^2) + 1 - 2*u[i]
        J[i,i+1] = 1/(dz^2) + c/(2*dz)
    end
    return J
end

"Vector function for v"
function Fv0(par, dz, v, c)
    F = Vector{Float64}(undef, par.Nz) # Pre-allocate F
    F[1] = v[1] - par.vf # Dirichlet condition on interface
    for i = 2:par.Nz-1
        F[i] = par.D*(v[i+1]-2*v[i]+v[i-1])/(dz^2) + c*(v[i+1]-v[i-1])/(2*dz) + par.λ*v[i]*(1-v[i])
    end
    F[par.Nz] = v[par.Nz] - 1.0 # Dirichlet condition at right boundary
    return F
end

"Jacobian matrix for v"
function Jv0(par, dz, v, c)
    J = zeros(par.Nz, par.Nz) # Pre-allocate Jacobian
    J[1,1] = 1.0 # Jacobian entry for left Dirichlet condition
    J[par.Nz, par.Nz] = 1.0 # Jacobian entry for right Dirichlet condition
    for i = 2:par.Nz-1
        J[i,i-1] = par.D/(dz^2) - c/(2*dz)
        J[i,i] = -2*par.D/(dz^2) + par.λ*(1 - 2*v[i])
        J[i,i+1] = par.D/(dz^2) + c/(2*dz)
    end
    return J
end

##### Linear stability analysis #####

"Implement shooting method for linear stability analysis"
function lsa(par, c, u0, v0, q, dz, ω, ϵ)
    # Implement Newton's method
    for j = 1:20
        ω_old = ω # Store previous guess
        # Solve boundary-value problems
        u = bvp_u1(par, dz, c, q, ω, u0)
        v = bvp_v1(par, dz, c, q, ω, v0)
        up = bvp_u1(par, dz, c, q, ω+ϵ, u0)
        vp = bvp_v1(par, dz, c, q, ω+ϵ, v0)
        # Implement Newton's method
        f = ω + par.κu*dudz(u,dz) + par.κv*dvdz(v,dz) # Interface condition
        fp = ω + ϵ + par.κu*dudz(up,dz) + par.κv*dvdz(vp,dz) # Perturbed interface condition
        d = (fp-f)/ϵ # Approximate df/dϵ
        ω = ω - f/d # Newton iteration
        if abs(ω - ω_old) < 1e-6
            return ω, bvp_u1(par, dz, c, q, ω, u0), bvp_v1(par, dz, c, q, ω, v0)
        end
    end
    @printf("LSA Shooting method did not converge.\n")
end

"Nonlinear boundary-value problem for u1"
function bvp_u1(par, dz, c, q, ω, u0)
    u = ones(par.Nz) # Initial guess
    # Newton's method
    for i = 1:20
        u_old = u # Store previous iteration u
        F = Fu1(par, dz, u, c, q, ω, u0) # Construct F
        J = Ju1(par, dz, c, q, ω, u0) # Construct Jacobian
        u = u - J\F # Newton iteration
        if norm(u - u_old) < 1e-6
            return u
        end
    end
    @printf("TWS BVP for u did not converge.\n")
end

"Nonlinear boundary value problem for v1"
function bvp_v1(par, dz, c, q, ω, v0)
    v = ones(par.Nz) # Initial guess
    # Newton's method
    for i = 1:20
        v_old = v # Store previous iteration u
        F = Fv1(par, dz, v, c, q, ω, v0) # Construct F
        J = Jv1(par, dz, c, q, ω, v0) # Construct Jacobian
        v = v - J\F # Newton iteration
        if norm(v - v_old) < 1e-6
            return v
        end
    end
    @printf("TWS BVP for v did not converge.\n")
end

"Vector function for u1"
function Fu1(par, dz, u, c, q, ω, u0)
    F = Vector{Float64}(undef, par.Nz) # Pre-allocate F
    F[1] = u[1]  # Dirichlet condition on left boundary
    for i = 2:par.Nz-1
        F[i] = (u[i+1]-2*u[i]+u[i-1])/(dz^2) + c*(u[i+1]-u[i-1])/(2*dz) + (1-ω-q^2-2*u0[i])*u[i] + (ω + q^2)*(u0[i+1]-u0[i-1])/(2*dz)
    end
    F[par.Nz] = u[par.Nz] + par.γ*q^2 # Dirichlet condition at interface
    return F
end

"Jacobian matrix for u1"
function Ju1(par, dz, c, q, ω, u0)
    J = zeros(par.Nz, par.Nz) # Pre-allocate Jacobian
    J[1,1] = 1.0 # Jacobian entry for left Dirichlet condition
    J[par.Nz, par.Nz] = 1.0 # Jacobian entry for right Dirichlet condition
    for i = 2:par.Nz-1
        J[i,i-1] = 1/(dz^2) - c/(2*dz)
        J[i,i] = -2/(dz^2) + (1-ω-q^2-2*u0[i])
        J[i,i+1] = 1/(dz^2) + c/(2*dz)
    end
    return J
end

"Vector function for v1"
function Fv1(par, dz, v, c, q, ω, v0)
    F = Vector{Float64}(undef, par.Nz) # Pre-allocate F
    F[1] = v[1] + par.γ*q^2 # Dirichlet condition on interface
    for i = 2:par.Nz-1
        F[i] = par.D*(v[i+1]-2*v[i]+v[i-1])/(dz^2) + c*(v[i+1]-v[i-1])/(2*dz) + (par.λ-ω-par.D*q^2-2*par.λ*v0[i])*v[i] + (ω+par.D*q^2)*(v0[i+1]-v0[i-1])/(2*dz)
    end
    F[par.Nz] = v[par.Nz] # Dirichlet condition at right boundary
    return F
end

"Jacobian matrix for v1"
function Jv1(par, dz, c, q, ω, v0)
    J = zeros(par.Nz, par.Nz) # Pre-allocate Jacobian
    J[1,1] = 1.0 # Jacobian entry for left Dirichlet condition
    J[par.Nz, par.Nz] = 1.0 # Jacobian entry for right Dirichlet condition
    for i = 2:par.Nz-1
        J[i,i-1] = par.D/(dz^2) - c/(2*dz)
        J[i,i] = -2*par.D/(dz^2) + (par.λ-ω-par.D*q^2-2*par.λ*v0[i])
        J[i,i+1] = par.D/(dz^2) + c/(2*dz)
    end
    return J
end

##### Run code #####
@time main()