"""
    benchmark_utils.jl

Minimal utilities for benchmarking ModelingToolkit compilation performance.
Only contains functions needed for benchmark_scaling.
"""

using ModelingToolkit: @variables, @parameters, D_nounits as D, ODESystem, ODEProblem
using Random

# Simplified model structure
mutable struct OceanicBoxModel
    N::Int
    reactive_species::Dict
    non_reactive_species::Dict
    equations::Vector{Equation}
    parameters::Dict
    initial_conditions::Dict
    system::Union{Nothing, ODESystem}
end

function create_oceanic_box_model(N::Int)
    model = OceanicBoxModel(N, Dict(), Dict(), [], Dict(), Dict(), nothing)
    create_variables!(model)
    create_parameters!(model)
    create_equations!(model)
    create_initial_conditions!(model)
    return model
end

function create_variables!(model::OceanicBoxModel)
    N = model.N
    t = ModelingToolkit.t_nounits
    
    @variables a(t)[1:N] b(t)[1:N] c(t)[1:N] d(t)[1:N]
    model.reactive_species[:a] = a
    model.reactive_species[:b] = b
    model.reactive_species[:c] = c
    model.reactive_species[:d] = d
    
    @variables begin
        F_adv_a(t)[1:N]; F_adv_b(t)[1:N]; F_adv_c(t)[1:N]; F_adv_d(t)[1:N]
        R_ab(t)[1:N]; R_cd(t)[1:N]; R_decay_a(t)[1:N]
        S_a(t)[1:N]; S_b(t)[1:N]; S_c(t)[1:N]; S_d(t)[1:N]
    end
    
    model.non_reactive_species[:F_adv_a] = F_adv_a
    model.non_reactive_species[:F_adv_b] = F_adv_b
    model.non_reactive_species[:F_adv_c] = F_adv_c
    model.non_reactive_species[:F_adv_d] = F_adv_d
    model.non_reactive_species[:R_ab] = R_ab
    model.non_reactive_species[:R_cd] = R_cd
    model.non_reactive_species[:R_decay_a] = R_decay_a
    model.non_reactive_species[:S_a] = S_a
    model.non_reactive_species[:S_b] = S_b
    model.non_reactive_species[:S_c] = S_c
    model.non_reactive_species[:S_d] = S_d
end

function create_parameters!(model::OceanicBoxModel)
    @parameters k_adv=0.1 k_ab=0.01 k_cd=0.005 k_decay=0.001 Km_ab=1.0 Km_cd=2.0 V_box=1.0
    model.parameters[:k_adv] = k_adv
    model.parameters[:k_ab] = k_ab  
    model.parameters[:k_cd] = k_cd
    model.parameters[:k_decay] = k_decay
    model.parameters[:Km_ab] = Km_ab
    model.parameters[:Km_cd] = Km_cd
    model.parameters[:V_box] = V_box
end

michaelis_menten_rate(S1, S2, k, Km) = k * S1 * S2 / (Km + S1 + S2)
@register_symbolic michaelis_menten_rate(S1, S2, k, Km)

function create_equations!(model::OceanicBoxModel)
    N = model.N
    equations = []
    
    a = model.reactive_species[:a]
    b = model.reactive_species[:b]
    c = model.reactive_species[:c]
    d = model.reactive_species[:d]
    
    F_adv_a = model.non_reactive_species[:F_adv_a]
    F_adv_b = model.non_reactive_species[:F_adv_b]
    F_adv_c = model.non_reactive_species[:F_adv_c]
    F_adv_d = model.non_reactive_species[:F_adv_d]
    R_ab = model.non_reactive_species[:R_ab]
    R_cd = model.non_reactive_species[:R_cd]
    R_decay_a = model.non_reactive_species[:R_decay_a]
    S_a = model.non_reactive_species[:S_a]
    S_b = model.non_reactive_species[:S_b]
    S_c = model.non_reactive_species[:S_c]
    S_d = model.non_reactive_species[:S_d]
    
    k_adv = model.parameters[:k_adv]
    k_ab = model.parameters[:k_ab]
    k_cd = model.parameters[:k_cd]
    k_decay = model.parameters[:k_decay]
    Km_ab = model.parameters[:Km_ab]
    Km_cd = model.parameters[:Km_cd]
    V_box = model.parameters[:V_box]
    
    for i in 1:N
        next_box = (i % N) + 1
        prev_box = ((i - 2 + N) % N) + 1
        
        push!(equations, F_adv_a[i] ~ k_adv * (a[prev_box] - 2*a[i] + a[next_box]))
        push!(equations, F_adv_b[i] ~ k_adv * (b[prev_box] - 2*b[i] + b[next_box]))
        push!(equations, F_adv_c[i] ~ k_adv * (c[prev_box] - 2*c[i] + c[next_box]))
        push!(equations, F_adv_d[i] ~ k_adv * (d[prev_box] - 2*d[i] + d[next_box]))
        
        push!(equations, R_ab[i] ~ michaelis_menten_rate(a[i], b[i], k_ab, Km_ab))
        push!(equations, R_cd[i] ~ michaelis_menten_rate(c[i], d[i], k_cd, Km_cd))
        push!(equations, R_decay_a[i] ~ k_decay * a[i])
        
        push!(equations, S_a[i] ~ F_adv_a[i] - R_ab[i] - R_decay_a[i])
        push!(equations, S_b[i] ~ F_adv_b[i] - R_ab[i])
        push!(equations, S_c[i] ~ F_adv_c[i] - R_cd[i] + 0.5 * R_ab[i])
        push!(equations, S_d[i] ~ F_adv_d[i] - R_cd[i] + 0.3 * R_decay_a[i])
        
        push!(equations, D(a[i]) ~ S_a[i] / V_box)
        push!(equations, D(b[i]) ~ S_b[i] / V_box)
        push!(equations, D(c[i]) ~ S_c[i] / V_box)
        push!(equations, D(d[i]) ~ S_d[i] / V_box)
    end
    
    model.equations = equations
end

function create_initial_conditions!(model::OceanicBoxModel; seed=1234)
    Random.seed!(seed)
    N = model.N
    u0_map = Dict()
    
    for (name, var) in model.reactive_species
        for i in 1:N
            u0_map[var[i]] = 1.0 + 0.1 * randn()
        end
    end
    
    model.initial_conditions = u0_map
end

function compile_system!(model::OceanicBoxModel)
    @named sys = System(model.equations, ModelingToolkit.t_nounits)
    compiled_sys = mtkcompile(sys)
    model.system = compiled_sys
    return compiled_sys
end

function benchmark_mtkcompile(N::Int)
    model = create_oceanic_box_model(N)
    original_equations = length(model.equations)
    
    time_compilation = @elapsed compile_system!(model)
    time_problem = @elapsed ODEProblem(model.system, model.initial_conditions, (0.0, 10.0))
    
    compiled_equations = length(equations(model.system))
    
    return (
        equations = (original=original_equations, compiled=compiled_equations),
        times = (compilation=time_compilation, problem=time_problem),
        N = N
    )
end

function benchmark_scaling(N_values; verbose=false)
    results = []
    
    if verbose
        println("Running scaling benchmark for N = $N_values")
        println()
    end
    
    for N in N_values
        if verbose; println("Benchmarking N = $N..."); end
        result = benchmark_mtkcompile(N)
        push!(results, result)
        if verbose; println("  Compilation time: $(round(result.times.compilation, digits=4)) s"); end
    end
    
    if verbose
        println("\n" * "="^80)
        println("SCALING ANALYSIS")
        println("="^80)
        println("| N    | Variables | Equations | Compilation Time (s) | Problem Time (s) |")
        println("|------|-----------|-----------|---------------------|------------------|")
        
        for result in results
            N = result.N
            vars = 15*N
            eqs = result.equations.original
            comp_time = result.times.compilation
            prob_time = result.times.problem
            println("| $(lpad(N,4)) | $(lpad(vars,9)) | $(lpad(eqs,9)) | $(lpad(round(comp_time, digits=4),19)) | $(lpad(round(prob_time, digits=4),16)) |")
        end
        println("="^80)
    end
    
    return results
end
