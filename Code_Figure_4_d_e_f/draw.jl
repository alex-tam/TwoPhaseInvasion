# Plot solutions to 2D Fisher-Stefan models
# Nizhum Rahman & Alex Tam, 19/06/2023

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
    #v_heatmap = CairoMakie.heatmap!(ax, x, y, v, colormap=:OrRd_3, colorrange=(0, 1.0))
    v_heatmap = CairoMakie.heatmap!(ax, x, y, v, colormap=:Blues_3, colorrange=(0, 1.0))
    u_colorbar = CairoMakie.Colorbar(fig[1, 2], u_heatmap, label=L"u", ticks=(0:0.2:1, [L"%$s" for s in 0:0.2:1]), labelsize=36, ticklabelsize=36)
    #v_colorbar = CairoMakie.Colorbar(fig[2, 1], v_heatmap, label=L"v", ticks=(0:0.2:1, [L"%$s" for s in 0:0.2:1]), labelsize=36, ticklabelsize=36, vertical=false, flipaxis=false)
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
    #v_heatmap = CairoMakie.heatmap!(ax, x, y, v, colormap=:OrRd_3, colorrange=(0, 1.0))
    v_heatmap = CairoMakie.heatmap!(ax, x, y, v, colormap=:Blues_3, colorrange=(0, 1.0))
    u_colorbar = CairoMakie.Colorbar(fig[1, 2], u_heatmap, label=L"u", ticks=(0:0.2:1, [L"%$s" for s in 0:0.2:1]), labelsize=36, ticklabelsize=36)
    #v_colorbar = CairoMakie.Colorbar(fig[2, 1], v_heatmap, label=L"v", ticks=(0:0.2:1, [L"%$s" for s in 0:0.2:1]), labelsize=36, ticklabelsize=36, vertical=false, flipaxis=false)
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
        Plots.plot!(x, ux, xlabel = L"$x$", ylabel = L"$u(x,10,t), v(x,10,t)$", linecolor=:black, linewidth = 1, aspect_ratio = Lx, grid = false, margin=3mm, legend=false, xlims=(0, Lx), ylims=(0,1.2))
       # Plots.plot!(x, vx, linecolor=:red, linewidth = 1)
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
        Plots.plot!(y, uy, xlabel = L"$y$", ylabel = L"$u(10,y,t), v(10,y,t)$", linecolor=:black, linewidth = 1, aspect_ratio = Ly, grid = false, margin=3mm, legend=false, xlims=(0, Ly), ylims=(0,1.2))
        #Plots.plot!(y, vy, linecolor=:red, linewidth = 1)
    end
    Plots.savefig("Slice_y.pdf")
end

"Control function for plotting"
function draw()
    # Import data
    plot_times = convert(Vector{Int}, vec(readdlm("plot_times.csv")))
    x = vec(readdlm("x.csv"))
    y = vec(readdlm("y.csv"))
    nx::Int = (length(x)+1)/2 # Index for slice plot
    ny::Int = (length(y)+1)/2 # Index for slice plot
    Lx = maximum(x)
    Ly = maximum(y)
    U = readdlm("U-0.csv")
    ϕ = readdlm("Phi-0.csv")
    # Plot
    draw_heat(x, y, U, ϕ, 0, Lx, Ly)
    for i in plot_times
        U = readdlm("U-$i.csv")
        V = readdlm("V-$i.csv")
        ϕ = readdlm("Phi-$i.csv")
        draw_heat(x, y, U, V, ϕ, i, Lx, Ly) # Heat maps
    end
    draw_slices(x, y, nx, ny, Lx, Ly, plot_times) # Slice plots
end

@time draw()