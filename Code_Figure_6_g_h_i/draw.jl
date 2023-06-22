# Plot solutions to 2D Fisher-Stefan models
# Nizhum Rahman & Alex Tam, 20/06/2023

using Plots
using CairoMakie
using Measures
using LaTeXStrings
using DelimitedFiles
using Printf
using Polynomials
using DifferentialEquations

"Plot solutions as 2D heat maps"
function draw_heat(x, y, U, ϕ, i, Lx, Ly)
    # Extract u and v populations
    u = [ ϕ[i,j] <= 0 ?  U[i,j] : U[i,j] for i = 1:201, j = 1:201]
    v = [ ϕ[i,j] > 0 ?  U[i,j] : NaN for i=1: 201, j =1: 201]
    # Generate and save density figure
    fig = CairoMakie.Figure(resolution = (800, 600), fontsize=36, font="Computer Modern")
    ax = CairoMakie.Axis(fig[1, 1], xlabel=L"x", ylabel=L"y", xticks=(0:5:Lx, [L"%$s" for s in 0:5:20]), yticks=(0:5:Ly, [L"%$s" for s in 0:5:20]), aspect=1)
    u_heatmap = CairoMakie.heatmap!(ax, x, y, u, colormap=:Blues_3, colorrange=(0, 1.0))
    v_heatmap = CairoMakie.heatmap!(ax, x, y, v, colormap=:OrRd_3, colorrange=(0, 1.0))
    u_colorbar = CairoMakie.Colorbar(fig[1, 2], u_heatmap, label=L"u", ticks=(0:0.2:1, [L"%$s" for s in 0:0.2:1]), labelsize=36, ticklabelsize=36)
    v_colorbar = CairoMakie.Colorbar(fig[2, 1], v_heatmap, label=L"v", ticks=(0:0.2:1, [L"%$s" for s in 0:0.2:1]), labelsize=36, ticklabelsize=36, vertical=false, flipaxis=false)
    CairoMakie.contour!(x, y, ϕ, levels=[0.0], color=:black, linewidth = 2)
    CairoMakie.colsize!(fig.layout, 1, Aspect(1, 1))
    CairoMakie.resize_to_layout!(fig)
    CairoMakie.save("Density-$i.pdf", fig)
end

function draw_heat(x, y, U, V, ϕ, i, Lx, Ly)
    # Extract u and v populations
    u = [ ϕ[i,j] <= 0 ? U[i,j] : U[i,j] for i=1:201, j=1:201]
    v = [ ϕ[i,j] > 0 ? U[i,j] : NaN for i=1:201, j=1:201]
    # Generate and save density figure
    fig = CairoMakie.Figure(resolution = (800, 600), fontsize=36, font="Computer Modern")
    ax = CairoMakie.Axis(fig[1, 1], xlabel=L"x", ylabel=L"y", xticks=(0:5:Lx, [L"%$s" for s in 0:5:20]), yticks=(0:5:Ly, [L"%$s" for s in 0:5:20]), aspect=1)
    u_heatmap = CairoMakie.heatmap!(ax, x, y, u, colormap=:Blues_3, colorrange=(0, 1.0))
    v_heatmap = CairoMakie.heatmap!(ax, x, y, v, colormap=:OrRd_3, colorrange=(0, 1.0))
    u_colorbar = CairoMakie.Colorbar(fig[1, 2], u_heatmap, label=L"u", ticks=(0:0.2:1, [L"%$s" for s in 0:0.2:1]), labelsize=36, ticklabelsize=36)
    v_colorbar = CairoMakie.Colorbar(fig[2, 1], v_heatmap, label=L"v", ticks=(0:0.2:1, [L"%$s" for s in 0:0.2:1]), labelsize=36, ticklabelsize=36, vertical=false, flipaxis=false)
    CairoMakie.contour!(x, y, ϕ, levels=[0.0], color=:black, linewidth = 2)
    CairoMakie.colsize!(fig.layout, 1, Aspect(1, 1))
    CairoMakie.resize_to_layout!(fig)
    CairoMakie.save("Density-$i.pdf", fig)
    # Plot extension velocity
    Plots.default(fontfamily = "Computer Modern", titlefontsize = 18, guidefontsize = 26, tickfontsize = 18, legendfontsize = 18)
    Plots.heatmap(x,y,transpose(V), xlabel = L"$x$", ylabel = L"$y$", margin=3mm, aspect_ratio=:equal, xlims=(0,Lx), ylims =(0,Ly), tick_direction=:out, c=:plasma)
    Plots.savefig("ExtensionVelocity-$i.pdf")
end

"Density slice plots"
function draw_slices(x, y, nx, ny, Lx, Ly, plot_times)
    Plots.gr() # Load GR plotting backend
    Plots.plot() # Clear previous plots
    Plots.default(fontfamily = "Computer Modern", titlefontsize = 18, guidefontsize = 26, tickfontsize = 18, legendfontsize = 18)
    # Plot slice in x-direction
    push!(plot_times, 0)
    for i in plot_times
        U = readdlm("U-$i.csv") # Import data
        ϕ = readdlm("Phi-$i.csv")
        # Separate u and v and obtain slices
        u = [ ϕ[i,j] <= 0 ? U[i,j] : 0 for i=1:201, j=1:201]
        v = [ ϕ[i,j] > 0 ? U[i,j] : 0 for i=1:201, j=1:201]
        ux = u[:, ny]
        vx = v[:, ny]
        # Plot data
        Plots.plot!(x, ux, xlabel = L"$x$", ylabel = L"$u(x,10,t), v(x,10,t)$", linecolor=:blue, linewidth = 1, aspect_ratio = Lx, grid = false, margin=3mm, legend=false, xlims=(0, Lx), ylims=(0,1.2))
        Plots.plot!(x, vx, linecolor=:red, linewidth = 1)
    end
    Plots.savefig("Slice_x.pdf")
    # Plot slice in y-direction
    Plots.plot() # Clear previous plots
    for i in plot_times
        U = readdlm("U-$i.csv") # Import data
        ϕ = readdlm("Phi-$i.csv")
        # Separate u and v and obtain slices
        u = [ ϕ[i,j] <= 0 ? U[i,j] : 0 for i=1:201, j=1:201]
        v = [ ϕ[i,j] > 0 ? U[i,j] : 0 for i=1:201, j=1:201]
        uy = u[nx, :]
        vy = v[nx, :]
        # Plot data
        Plots.plot!(y, uy, xlabel = L"$y$", ylabel = L"$u(10,y,t), v(10,y,t)$", linecolor=:blue, linewidth = 1, aspect_ratio = Ly, grid = false, margin=3mm, legend=false, xlims=(0, Ly), ylims=(0,1.2))
        Plots.plot!(y, vy, linecolor=:red, linewidth = 1)
    end
    Plots.savefig("Slice_y.pdf")
end

"Compute growth rate"
function draw_growth(t, Amp, t_min::Float64)
    Plots.gr()
    Plots.plot() # Load GR plotting backend and clear previous plots
    Plots.default(fontfamily = "Computer Modern", titlefontsize = 14, guidefontsize = 20, tickfontsize = 14, legendfontsize = 14)
    ind::Int = 0 # Index corresponding to t = 0.1
    for i in eachindex(t)
        if t[i] <= t_min
            ind += 1
        end
    end
    nn=10
    poly = fit(t[ind+nn:end-nn], log.(Amp)[ind+nn:end-nn], 1) # Fit straight line to data
    Plots.plot(t, log.(Amp), xlabel = L"$t$", ylabel = L"$\log(A)$", label = "Numerical Data", linewidth = 2, margin = 5mm) # Plot data
    Plots.scatter!([t[ind+nn], t[end-nn]], [log(Amp[ind+nn]), log(Amp[end-nn])], markersize = 5, markershape = :xcross, markercolor = :red, label = false) # Scatter plot of t_min
    Plots.plot!(t, poly.coeffs[2].*t .+ poly.coeffs[1], label = "Linear Fit", linestyle=:dash, linewidth = 2,legend=:topleft) # Plot linear trendline
    ω = poly.coeffs[2] # Obtain slope
    @printf("The numerical growth rate is: %f.\n", ω)
    Plots.savefig("perturbation_amplitude.pdf")
end

"Plot interface position versus time"
function draw_interface(t, L, ε, T, x, y, dx, Lx, Ly, Nx, Ny, plot_times)
    Plots.gr(); Plots.plot() # Load GR plotting backend and clear previous plots
    Plots.default(titlefont = (18, "Computer Modern"), guidefont = (26, "Computer Modern"), tickfont = (18, "Computer Modern"))
    Plots.plot(t, L, xlabel = L"$t$", ylabel = L"$L(t)$", margin=3mm, xlims=(0, maximum(t)), ylims=(0, maximum(L)), legend = false)
    Plots.savefig("L.pdf")
    if ε == 0.1
        @printf("Numerical travelling wave speed is %f.\n", (L[end]-L[1])/T)
    end
    # Optional: Draw interface
    Plots.plot() # Clear previous plots
    plot_interface(x, y, dx, Lx, Ly, Nx, Ny, plot_times)
    Plots.savefig("interface.pdf")
end

"Draw interface as series of points"
function plot_interface(x, y, dx, Lx, Ly, Nx, Ny, plot_times)
    pt = vcat(0, plot_times)
    for k in pt
        # Import level-set function at relevant time
        ϕ = readdlm("Phi-$k.csv")
        # Pre-allocate empty vectors
        xi = Vector{Float64}()
        yi = Vector{Float64}()
        # Locate interface grid points
        for j = 1:Ny
            ϕv = ϕ[:,j] # Obtain 1D vector of ϕ
            for i = 1:Nx
                if (ϕv[i] < 0) && (ϕv[i+1] >= 0)
                    θ = ϕv[i]/(ϕv[i] - ϕv[i+1])
                    push!(xi, x[i] + θ*dx)
                    push!(yi, y[j])
                end
            end
        end
        # Plot interface
        Plots.scatter!([xi],[yi], xlabel = L"$x$", ylabel = L"$y$", margin=3mm, xlims=(0,Lx), ylims=(0,Ly), color="green", legend = false, aspect_ratio=:equal, markersize=1)
    end
    Plots.savefig("interface.pdf")
end

"Control function for plotting"
function draw()
    # Import data
    plot_times = convert(Vector{Int}, vec(readdlm("plot_times.csv")))
    x = vec(readdlm("x.csv")); dx = x[2] - x[1]
    y = vec(readdlm("y.csv"))
    nx::Int = (length(x)+1)/2 # Index for slice plot
    ny::Int = (length(y)+1)/2 # Index for slice plot
    Nx = length(x); Ny = length(y)
    Lx = maximum(x); Ly = maximum(y)
    U = readdlm("U-0.csv")
    ϕ = readdlm("Phi-0.csv")
    t = vec(readdlm("t.csv"))
    L = vec(readdlm("L.csv"))
    Amp = vec(readdlm("Amp.csv"))
    ε = 0.1
    T = 30.0
    # Plot
    draw_heat(x, y, U, ϕ, 0, Lx, Ly)
    for i in plot_times
        U = readdlm("U-$i.csv")
        V = readdlm("V-$i.csv")
        ϕ = readdlm("Phi-$i.csv")
        draw_heat(x, y, U, V, ϕ, i, Lx, Ly) # Heat maps
    end
    draw_slices(x, y, nx, ny, Lx, Ly, plot_times) # Slice plots
    draw_interface(t, L, ε, T, x, y, dx, Lx, Ly, Nx, Ny, plot_times) # Interface plots
    #draw_growth(t, Amp, ε) # Growth rate plot
end

@time draw()