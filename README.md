# MTK Compilation Benchmarking Suite

A comprehensive benchmarking suite for evaluating ModelingToolkit (MTK) compilation performance on large-scale oceanic box models. This suite is designed to help analyze the performance characteristics of MTK's `structural_simplify` function (equivalent to `mtkcompile`) on realistic advection-reaction systems.

## Overview

This benchmarking suite creates generic oceanic box models that mimic the structure of real geochemical models like GEOCLIM. The models include:

- **Reactive species**: 4 species (a, b, c, d) each with N boxes representing different ocean regions
- **Linear terms**: Advection between boxes using ring topology  
- **Non-linear terms**: Michaelis-Menten kinetics for species interactions
- **Auxiliary variables**: Flux and reaction rate variables (non-reactive species)

## Quick Start

```julia
# Load the environment and run basic test
include("mwe.jl")

# Run full demonstration
include("demo.jl")
```

## Key Functions

### `benchmark_mtkcompile(N; verbose=true)`

Benchmarks a single system of size N, timing:
1. **Model creation**: Building equations and variables
2. **System compilation**: MTK's `structural_simplify` step
3. **ODEProblem creation**: Final problem setup

**Example:**
```julia
result = benchmark_mtkcompile(10)
# Output includes timing breakdown, equation counts, and performance metrics
```

### `benchmark_scaling(N_values; verbose=false)`

Compares performance across multiple system sizes to analyze scaling behavior.

**Example:**
```julia
results = benchmark_scaling([5, 10, 20, 50])
# Produces scaling analysis table showing performance vs system size
```

### `create_oceanic_box_model(N)`

Creates a model without benchmarking for manual analysis.

**Example:**
```julia
model = create_oceanic_box_model(10)
# Returns OceanicBoxModel struct with all components
```

## Model Structure

### Reactive Species (ODEs)
- `a[1:N], b[1:N], c[1:N], d[1:N]` - Concentrations in each box
- Each follows: `D(species[i]) = Sources[i] / Volume`

### Non-Reactive Species (Algebraic)
- **Advection fluxes**: `F_adv_*[1:N]` - Transport between boxes
- **Reaction rates**: `R_ab[1:N], R_cd[1:N], R_decay_a[1:N]` - Reaction kinetics
- **Source terms**: `S_*[1:N]` - Net sources/sinks for each species

### Model Equations
1. **Linear advection**: Ring topology with nearest-neighbor exchange
2. **Non-linear reactions**: 
   - `a + b → products` (Michaelis-Menten kinetics)
   - `c + d → products` (Michaelis-Menten kinetics)  
   - `a → decay` (First-order decay)
3. **Coupling**: Products from one reaction feed other species

## Benchmark Output

```
================================================================================
MTK Compilation Benchmark for Oceanic Box Model
================================================================================
System size (N): 10
Total reactive species variables: 40
Total auxiliary variables: 110

Step 1: Creating model equations...
  ✓ Created 150 equations
  ✓ Time: 0.0123 seconds

Step 2: Compiling system (structural_simplify)...
  ✓ Reduced to 40 equations  
  ✓ Reduction: 110 equations eliminated
  ✓ Time: 1.2345 seconds

Step 3: Creating ODEProblem...
  ✓ ODEProblem created
  ✓ Time: 0.0456 seconds

================================================================================
BENCHMARK SUMMARY
================================================================================
System size (N):                  10
Total variables:                   150
Original equations:                150
Compiled equations:                40
Equations eliminated:              110

Timing breakdown:
  Model creation:                  0.0123 s
  System compilation:              1.2345 s
  ODEProblem creation:             0.0456 s
  Total time:                      1.2924 s

Performance metrics:
  Compilation time per equation:   8.23 ms
  Equations processed per second:  121.5
================================================================================
```

## Files

- `mwe.jl` - Main entry point with quick test
- `src/benchmark_utils.jl` - Core benchmarking functions
- `demo.jl` - Full demonstration script
- `ModelingToolkit/Project.toml` - Dependencies

## Dependencies

- ModelingToolkit.jl - Symbolic modeling framework
- OrdinaryDiffEq.jl - ODE solver suite
- BenchmarkTools.jl - Performance measurement
- Random.jl - Reproducible random numbers
- LinearAlgebra.jl - Matrix operations

## Inspiration

This benchmarking suite is inspired by the [GEOCLIM_mtk](https://github.com/sablonl/GEOCLIM_mtk) project but uses:
- Latest ModelingToolkit version (not JuliaSimCompiler)
- Generic model structure for benchmarking
- Focus on compilation performance analysis

## Usage Tips

1. **Start small**: Begin with N=5-10 to understand the output format
2. **Memory considerations**: Large N values (>100) may require significant RAM
3. **Reproducibility**: Fixed random seed ensures consistent results
4. **Scaling analysis**: Use `benchmark_scaling()` to identify performance bottlenecks

## Customization

The model structure can be easily modified in `src/benchmark_utils.jl`:
- Add more reactive species
- Change reaction kinetics
- Modify advection topology
- Adjust parameter values

This provides a flexible framework for testing MTK performance on different model types.
