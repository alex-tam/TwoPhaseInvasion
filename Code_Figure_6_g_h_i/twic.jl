# Solve Two-Phase BVP to generate initial condition
# Nizhum Rahman & Alex Tam, 20/06/2023

"Main function for travelling wave and first-order correction"
function twic(par)
    # 1. Parameters and domain
    zl = range(-10.0, 0.0, length = par.Nz)
    zr = range(0.0, 10.0, length = par.Nz)
    dz = zl[2] - zl[1]
    # 2. Obtain travelling wave solution using shooting method
    c, u0, v0 = tws(par, dz, 0.0, 1e-6)
    # 3. Obtain dispersion relation using shooting method
    ω, u1, v1 = lsa(par, c, u0, v0, par.q, dz, 1.0, 1e-6)
    return zl, zr, u0, u1, v0, v1
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
        F[i] = par.Du*(u[i+1]-2*u[i]+u[i-1])/(dz^2) + c*(u[i+1]-u[i-1])/(2*dz) + par.λu*u[i]*(1-u[i])
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
        J[i,i-1] = par.Du/(dz^2) - c/(2*dz)
        J[i,i] = -2*par.Du/(dz^2) + par.λu*(1 - 2*u[i])
        J[i,i+1] = par.Du/(dz^2) + c/(2*dz)
    end
    return J
end

"Vector function for v"
function Fv0(par, dz, v, c)
    F = Vector{Float64}(undef, par.Nz) # Pre-allocate F
    F[1] = v[1] - par.vf # Dirichlet condition on interface
    for i = 2:par.Nz-1
        F[i] = par.Dv*(v[i+1]-2*v[i]+v[i-1])/(dz^2) + c*(v[i+1]-v[i-1])/(2*dz) + par.λv*v[i]*(1-v[i])
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
        J[i,i-1] = par.Dv/(dz^2) - c/(2*dz)
        J[i,i] = -2*par.Dv/(dz^2) + par.λv*(1 - 2*v[i])
        J[i,i+1] = par.Dv/(dz^2) + c/(2*dz)
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
        F[i] = par.Du*(u[i+1]-2*u[i]+u[i-1])/(dz^2) + c*(u[i+1]-u[i-1])/(2*dz) + (par.λu-ω-par.Du*q^2-2*par.λu*u0[i])*u[i] + (ω + par.Du*q^2)*(u0[i+1]-u0[i-1])/(2*dz)
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
        J[i,i-1] = par.Du/(dz^2) - c/(2*dz)
        J[i,i] = -2*par.Du/(dz^2) + (par.λu-ω-par.Du*q^2-2*par.λu*u0[i])
        J[i,i+1] = par.Du/(dz^2) + c/(2*dz)
    end
    return J
end

"Vector function for v1"
function Fv1(par, dz, v, c, q, ω, v0)
    F = Vector{Float64}(undef, par.Nz) # Pre-allocate F
    F[1] = v[1] + par.γ*q^2 # Dirichlet condition on interface
    for i = 2:par.Nz-1
        F[i] = par.Dv*(v[i+1]-2*v[i]+v[i-1])/(dz^2) + c*(v[i+1]-v[i-1])/(2*dz) + (par.λv-ω-par.Dv*q^2-2*par.λv*v0[i])*v[i] + (ω+par.Dv*q^2)*(v0[i+1]-v0[i-1])/(2*dz)
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
        J[i,i-1] = par.Dv/(dz^2) - c/(2*dz)
        J[i,i] = -2*par.Dv/(dz^2) + (par.λv-ω-par.Dv*q^2-2*par.λv*v0[i])
        J[i,i+1] = par.Dv/(dz^2) + c/(2*dz)
    end
    return J
end