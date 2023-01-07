# This script plots the distribution of the orbits of
# all asteroids from the Minor Planet Center (MPC)
# database and simulates the dynamic of a family of
# asteroid orbiting Jupiter
#
# The main goal is observation of Kirkwood gaps in
# the MPC data and of stable / unstable orbits in
# the results of the simulation

"Parse the string, taking more arguments if it's quoted"
function parse_string(i)::String
    # Start from the first argument after the flag
    j = i
    # If the arguments starts with an apostrophe,
    s = if startswith(ARGS[i], "'")
        # Find the index of the argument
        # which ends with an apostrophe
        while !endswith(ARGS[j], "'")
            j += 1
        end
        # Join the arguments in one string
        # and remove the apostrophes
        chop(join(ARGS[i:j], ' '), head=1, tail=1)
    else
        # Return the next argument
        ARGS[i]
    end
    return s
end

# Define default values for optional arguments
POSTFIX = ""
FORCE = false

# Parse the options
for i in eachindex(ARGS)
    # A postfix for the names of output files
    if ARGS[i] == "--postfix"
        try
            global POSTFIX = " ($(parse_string(i+1)))"
        catch
            println("Couldn't parse the value of the `--postfix` argument.")
            exit(1)
        end
    end
end

# Prepare color codes
RESET = "\e[0m"
GREEN = "\e[32m"
YELLOW = "\e[33m"

# Check for required arguments
if "--help" in ARGS
    println("""
        $(YELLOW)USAGE:$(RESET)
        { julia --project=. | ./julia.bash } scripts/script.jl [--postfix <POSTFIX>]

        $(YELLOW)OPTIONS:$(RESET)
            $(GREEN)--postfix <POSTFIX>$(RESET)    A postfix for the names of output files"""
    )
    exit(1)
end

"Floating point type used across the script"
F = Float64

"Padding in the output"
pad = " "^4

# Define the paths
CURRENT_DIR = @__DIR__
ROOT_DIR = basename(CURRENT_DIR) == "scripts" ? dirname(CURRENT_DIR) : CURRENT_DIR
PLOTS_DIR = joinpath(ROOT_DIR, "plots")
DATA_FILE = joinpath(ROOT_DIR, "data", "data.json")
DATA_CUT_FILE = joinpath(ROOT_DIR, "data", "data_cut.json")

# Make sure the needed directories exist
mkpath(PLOTS_DIR)

println('\n', pad, "> Loading the packages...")

using DifferentialEquations
using Distributions
using JSON
using LaTeXStrings
using LinearAlgebra
using Plots
using Random

# Use the GR backend for plots
pgfplotsx()

# Change some of the default parameters for plots
default(
    fontfamily="Computer Modern",
    dpi=300,
    size=(400, 400),
)

println(pad, "> Loading the MPC data...")

# Load the data
data = if isfile(DATA_CUT_FILE)
    JSON.parsefile(DATA_CUT_FILE)
else
    # Save the cut version
    data = JSON.parsefile(DATA_FILE)[1:50000]
    open(DATA_CUT_FILE, "w") do io
        JSON.print(io, data)
    end
    data
end
a = map(e -> e["a"], data)
e = map(e -> e["e"], data)
i = map(e -> e["i"], data)
M = map(e -> deg2rad(e["M"]), data)

# Compute the true anomaly
ν = @. M +
       (2 * e - e^3 / 4) * sin(M) +
       5 / 4 * e^2 * sin(2 * M) +
       13 / 12 * e^3 * sin(3 * M)

# Define the parameters of the visualization (MPC)
r_steps = 100
r_start = 2.0
r_end = 3.5
r_h = (r_end - r_start) / r_steps

"Get indices of the groups"
function groups_indices(vec)
    inner = findall(a -> a < 2.5, vec)
    intermediate = findall(a -> 2.5 <= a < 2.82, vec)
    outer = findall(a -> 2.82 <= a, vec)
    return (inner, intermediate, outer)
end

# Get indices of the groups
inner, intermediate, outer = groups_indices(a)

"Add a vertical line to the plot"
function vline!(x, label)
    linewidth = 0.5
    ymax = ylims()[2]
    color = :black
    plot!(
        [x, x],
        [ymax * 0.65, 0],
        label="",
        linestyle=:dash;
        color,
        linewidth
    )
    scatter!(
        [x],
        [ymax * 0.675],
        label="",
        marker=:dtriangle,
        markerstrokewidth=linewidth;
        color
    )
    annotate!([x], [ymax * 0.71], text(label, 9))
end

println(pad, "> Plotting the MPC histogram...")

"Plot the histogram of the distribution by semi-major axes"
function plot_distribution(data, var, r_start, r_h, r_end, add_lines)
    # Plot the histogram of the distribution
    histogram(
        data,
        label=[
            "Inner main-belt";;
            "Intermediate main-belt";;
            "Outer main-belt"
        ],
        xlabel=latexstring("$(var)" * " \\; \\mathrm{[AU]}"),
        ylabel=L"N",
        bins=range(r_start, r_end; step=r_h),
        fillcolor=[1 2 3],
        xlims=(r_start, r_end),
        ylims=(0, Inf),
        legend=:topright,
        linewidth=0.1,
        linecolor=:match,
        grid=false,
        tick_direction=:out,
    )
    # Add notable real world resonances
    if add_lines
        vline!(2.502, "3:1")
        vline!(2.825, "5:2")
        vline!(2.958, "7:3")
        vline!(3.279, "2:1")
    end
end

# Plot the distribution
plot_distribution([a[inner], a[intermediate], a[outer]], :a, r_start, r_h, r_end, true)

# Save the figure
savefig(joinpath(PLOTS_DIR, "mpc_histogram$(POSTFIX).pdf"))

println(pad, "> Plotting the scatter plot...")

"Add the label to the plot's legend"
function label(label, color_index)
    scatter!(
        [0],
        [maximum(a) + 5],
        color=palette(:default).colors[color_index],
        markerstrokewidth=-1,
        legend_font_valign=:bottom;
        label
    )
end

# Plot the projections
left = 0.0
right = 4.0
xrange = range(0, 2π; step=π / 4)
scatter(
    [ν[inner], ν[intermediate], ν[outer]],
    [a[inner], a[intermediate], a[outer]],
    proj=:polar,
    label="",
    color=[1 2 3],
    ylims=(left, right),
    xticks=(xrange, map(x -> "$(rad2deg(x))°", xrange)),
    yticks=range(left, right - 1, step=1),
    ydraw_arrow=true,
    markersize=0.5,
    markerstrokewidth=-1,
    legend=(1.0, 1.0),
)
annotate!([4.0], [right - 0.4], text(L"a", 9))
label("Inner main-belt", 1)
label("Intermediate main-belt", 2)
label("Outer main-belt", 3)

# Save the figure
savefig(joinpath(PLOTS_DIR, "mpc_scatter$(POSTFIX).pdf"))

# Mark the data for garbage collection
data = nothing

println(pad, "> Running the simulation...")

"Convert the `m^3/s^2` value to `AU^3/days^2"
convert(x) = x * (24 * 60 * 60)^2 / 149597870700.0^3

# Define the parameters

# Standard gravitational parameter of the Sun [AU^3/days^2]
μ₁ = convert(1.32712440018e20)

# Standard gravitational parameter of the Jupiter [AU^3/days^2]
μ₂ = convert(1.26686534e17)

# Distance between the stationary bodies [AU]
R = 5.2

# Compute the constant orbital angular velocity [1/days]
ω = √((μ₁ + μ₂) / R^3)

# Compute the center of mass and
# positions of the stationary bodies
α = μ₂ / (μ₁ + μ₂)
x₁ = -α * R
x₂ = (1 - α) * R

"Compute the distance between the minor
body and any of the stationary bodies"
ρ(x, y, xₛ) = √((x - xₛ)^2 + y^2)

"The differential equation of the circular restricted
three-body problem with a co-rotating frame"
function df!(df, state, _, _)
    x, y, u, v = state
    ρ₁ = ρ(x, y, x₁)
    ρ₂ = ρ(x, y, x₂)
    df[1] = u
    df[2] = v
    df[3] = -μ₁ * (x - x₁) / ρ₁^3 - μ₂ * (x - x₂) / ρ₂^3 + ω^2 * x + 2 * ω * v
    df[4] = -μ₁ * y / ρ₁^3 - μ₂ * y / ρ₂^3 + ω^2 * y - 2 * ω * u
end

# Initialize a pseudo-random number generator
rng = Xoshiro(1)

# Define the parameters of the simulation
n = 10000
t_years = 1000
h = 14
tspan = (0.0, t_years * 365.0 / h)

# Define the parameters of the visualization (Simulation)
r_steps = 100
r_start = 2.0
r_end = 3.5
r_h = (r_end - r_start) / r_steps

# Define the limits
lim = 2 * R

"Compute a vector of initial values using a pseudo-random
number generator with uniform distribution"
function u₀()::Vector{F}
    # Get random polar coordinates
    θ = rand(rng, Uniform(0, 2π))
    r = sqrt(rand(rng, Uniform(r_start^2, r_end^2)))
    # Convert them to the Cartesian coordinate system
    x = r * cos(θ)
    y = r * sin(θ)
    # Set the initial velocities
    u = (√(μ₁ / r) - ω * R) * sin(θ)
    v = (√(μ₁ / r) - ω * R) * cos(θ)
    return [x, y, u, v]
end

# Prepare storage for initial values
init = Matrix{F}(undef, n, 4)
for i in 1:n
    init[i, :] = u₀()
end

# Prepare storage for the counts (in fractions)
counts = Vector{F}(undef, r_steps)

"Compute the distance"
distance(x, y) = √(x^2 + y^2)

"Change the initial values in the problem"
function prob_func(prob, i, _)
    remake(prob, u0=init[i, :])
end

# Define the problem
ode_prob = ODEProblem{true}(df!, init[1, :], tspan)
# Create an ensemble
ensemble_prob = EnsembleProblem(ode_prob; prob_func)
# Run the simulation
sol = solve(
    ensemble_prob,
    Tsit5(),
    trajectories=n,
    adaptive=false,
    dt=h,
    reltol=1e-8,
    save_start=true,
    save_everystep=false,
)
# Get the results
x₀ = sol[1, begin, :]
y₀ = sol[2, begin, :]
r₀ = distance.(x₀, y₀)
x = sol[1, end, :]
y = sol[2, end, :]
r = distance.(x, y)
i_inside = findall(r -> r <= lim, r)
x_inside = x[i_inside]
y_inside = y[i_inside]

# For each bin
for k in 1:r_steps
    # Find the indices for the bin
    indices = findall(r -> r_start + (k - 1) * r_h <= r < r_start + k * r_h, r₀)
    # Count how many are still in the limits (as a fraction)
    counts[k] = count(r -> r <= lim, r[indices]) / length(indices)
end

# Imitate the input data for the histogram
inner, intermediate, outer = groups_indices(r)

println(pad, "> Plotting the results of the simulation...")

"Plot the state of the system"
function plot_state(x, y)
    # Plot the state
    scatter(
        x,
        y,
        xlabel=L"x",
        ylabel=L"y",
        legend=nothing,
        markersize=2,
        markerstrokewidth=-1,
        xlims=(-lim, lim),
        ylims=(-lim, lim),
    )
    # Mark the position of the Sun
    scatter!([x₁], [0])
    # Mark the position of the Jupiter
    scatter!([x₂], [0], markersize=3)
end

println(pad, pad, "...the initial state...")

# Plot the initial values
plot_state(init[:, 1], init[:, 2])

# Save the figure
savefig(joinpath(PLOTS_DIR, "model_initial$(POSTFIX).pdf"))

println(pad, pad, "...the final state...")

# Plot the final state
plot_state(x_inside, y_inside)

# Save the figure
savefig(joinpath(PLOTS_DIR, "model_final$(POSTFIX).pdf"))

println(pad, pad, "...the distribution histogram...")

# Plot the distribution histogram
plot_distribution([r[inner], r[intermediate], r[outer]], :r, r_start, r_h, r_end, false)

# Save the figure
savefig(joinpath(PLOTS_DIR, "model_distribution$(POSTFIX).pdf"))

println(pad, pad, "...the fractions histogram...")

"Return a rectangle shape"
rectangle(w, h, x, y) = Shape(x .+ [0, w, w, 0], y .+ [0, 0, h, h])

# Plot the fractions histogram
scatter(
    repeat([[-1.0]], 3),
    repeat([[-1.0]], 3),
    xlims=(r_start, r_end),
    ylims=(0, Inf),
    label=[
        "Inner main-belt";;
        "Intermediate main-belt";;
        "Outer main-belt"
    ],
    color = [1 2 3],
    xlabel=L"r \; \mathrm{[AU]}",
    ylabel="Still in limits (fraction)",
    linewidth=0.1,
    grid=false,
    tick_direction=:out,
    legend = :topright,
)
for k in 1:r_steps
    # Compute the `x` coordinate
    rec_x = r_start + (k - 1) * r_h
    # Define the color
    colors = palette(:default).colors
    color = if rec_x < 2.5
        colors[1]
    elseif 2.5 <= rec_x <= 2.82
        colors[2]
    else
        colors[3]
    end
    # Plot the rectangle
    plot!(
        rectangle(r_h, counts[k], rec_x, 0),
        label="",
        linecolor=:match;
        color,
    )
end

# Save the figure
savefig(joinpath(PLOTS_DIR, "model_fractions$(POSTFIX).pdf"))

# Mark the results for garbage collection
GC.gc(false)

println()
