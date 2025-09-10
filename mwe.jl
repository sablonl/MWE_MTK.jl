using Pkg;
Pkg.activate("ModelingToolkit");

using ModelingToolkit
using Random
using Plots

"""
MTK Compilation Benchmarking Suite for Oceanic Box Models

Streamlined version focused on benchmark_scaling with integrated plotting.
"""

include("src/benchmark_utils.jl")

# Run benchmarks
println("MTK Compilation Benchmarking Suite")
println("==================================")
println()

# Precompilation run
benchmark_scaling([2]; verbose = false);

# Main benchmarking run
results = benchmark_scaling(10:10:200; verbose = true)

# Extract data for plotting
equations_data = [r.equations.original for r in results]
compilation_times = [r.times.compilation for r in results]
problem_times = [r.times.problem for r in results]

# Create plot
p = plot(equations_data, compilation_times,
	label = "Compilation Time",
	linewidth = 2,
	marker = :circle,
	markersize = 6,
	xlabel = "Number of Equations",
	ylabel = "Time (seconds)",
	title = "MTK Performance: Equations vs Compilation Time",
	legend = :topleft,
	grid = true)

plot!(p, equations_data, problem_times,
	label = "ODEProblem Creation Time",
	linewidth = 2,
	marker = :square,
	markersize = 6)

# Fit compilation_times = a*equations_data^3 + b*equations_data^2 + c*equations_data + d
a, b, c, d = [equations_data .^ 3 equations_data .^ 2 equations_data ones(length(equations_data))] \ compilation_times
poly_fit_line(x) = @. a * x^3 + b * x^2 + c * x + d
plot!(p, poly_fit_line, extrema(equations_data)...,
	label = "Cubic Fit",
	linewidth = 2,
	linestyle = :dot)
savefig(p, "mtk_benchmark_plot.png")

poly_fit_line(10_000) / 3600
