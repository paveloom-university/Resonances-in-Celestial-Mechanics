# This script plots the distribution of the orbits of
# all asteroids from the Minor Planet Center (MPC)
# database and models the dynamic of a family of
# minor planets orbiting Jupiter

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

using JSON
using LaTeXStrings
using Plots

# Use the GR backend for plots
pgfplotsx()

# Change some of the default parameters for plots
default(
    fontfamily="Computer Modern",
    dpi=300,
    size=(400, 400),
    legend=:topright,
)

println(pad, "> Loading the data...")

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

# Get indices of the groups
inner = findall(a -> a < 2.5, a)
intermediate = findall(a -> 2.5 <= a <= 2.82, a)
outer = findall(a -> 2.82 < a, a)

"Add a vertical line to the plot"
function vline(x, label)
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

println(pad, "> Plotting the histogram...")

# Plot the distribution
left = 2.0
right = 3.5
histogram(
    [a[inner], a[intermediate], a[outer]],
    label=[
        "Inner main-belt";;
        "Intermediate main-belt";;
        "Outer main-belt"
    ],
    xlabel=L"a \; \mathrm{[AU]}",
    ylabel=L"N",
    bins=range(left, right, step=0.005),
    fillcolor=[1 2 3],
    xlims=(left, right),
    ylims=(0, Inf),
    linewidth=0.1,
    linecolor=:match,
    grid=false,
    tick_direction=:out,
)
vline(2.502, "3:1")
vline(2.825, "5:2")
vline(2.958, "7:3")
vline(3.279, "2:1")

# Save the figure
savefig(joinpath(PLOTS_DIR, "histogram$(POSTFIX).pdf"))

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
savefig(joinpath(PLOTS_DIR, "scatter$(POSTFIX).pdf"))

# Mark the data for garbage collection
data = nothing

println()
