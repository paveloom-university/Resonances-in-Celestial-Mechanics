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

"Floating point type used across the script"
F = Float64

"Integer type used across the script"
I = Int64

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

println(pad, "> Plotting the distribution...")

function vline(x::F, label::String) where {F<:Real}
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

# Plot the distribution
left = 2.0
right = 3.5
histogram(
    [
        filter(a -> left <= a < 2.5, a),
        filter(a -> 2.5 <= a <= 2.82, a),
        filter(a -> 2.82 < a <= right, a),
    ],
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
savefig(joinpath(PLOTS_DIR, "a$(POSTFIX).pdf"))

# Mark the data for garbage collection
data = nothing

println()
