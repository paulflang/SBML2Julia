using CSV
using DataFrames
using Ipopt
using JuMP

t_steps = 101 # Setting number of ODE discretisation steps

# Data
println("Reading measurement data...")
data_path = "/media/sf_DPhil_Project/Project07_Parameter Fitting/df_software/DisFit/tests/fixtures/0015_objectivePrior/_measurements.tsv"
df = CSV.read(data_path; use_mmap=false)
insert!(df, 1, (1:length(df[:,1])), :id)
dfg = groupby(df, :simulationConditionId)
data = []
for condition in keys(dfg)
    push!(data,unstack(dfg[condition], :time, :observableId, :measurement))
end

k_to_time_idx = []
for c in 1:length(keys(dfg))
    push!(k_to_time_idx, [])
    j = 0
    prev = "a"
    for t in dfg[c][:, :time]
        if t != prev
            j = j+1
        end
        push!(k_to_time_idx[c], j)
        prev = t
    end
end

t_exp = Vector(DataFrame(groupby(dfg[1], :observableId)[1])[!, :time])
t_sim = range(0, stop=t_exp[end], length=t_steps)
t_sim_to_exp = []
for i in 1:length(t_exp)
    idx = argmin(abs.(t_exp[i] .- t_sim))
    append!(t_sim_to_exp, idx)
end

results = Dict()
results["objective_value"] = Dict()
results["parameters"] = Dict()
results["species"] = Dict()
results["observables"] = Dict()

j_to_cond_par = [1, 2]
cond_without_preequ = [1]
cond_with_preequ = [2]
preequ = [3]

i_start = 1
    m = Model(with_optimizer(Ipopt.Optimizer))


    # Define global parameters
    println("Defining global parameters...")
    @variable(m, 0 <= a0 <= 10, start=0+(10-0)*rand(Float64))
    @variable(m, 0 <= b0 <= 10, start=0+(10-0)*rand(Float64))
    @variable(m, 0 <= k1_free <= 10, start=0+(10-0)*rand(Float64))
    @variable(m, 0 <= k2 <= 10, start=0+(10-0)*rand(Float64))
    @variable(m, 0 <= noise_A1 <= 10, start=0+(10-0)*rand(Float64))
    @variable(m, 0 <= noise_A2 <= 10, start=0+(10-0)*rand(Float64))
    @variable(m, noise_B == 0.1, start=0.1)
    @variable(m, 0 <= offset_B <= 10, start=0+(10-0)*rand(Float64))
    @variable(m, 0 <= scaling_B <= 10, start=0+(10-0)*rand(Float64))

    # Define condition-specific parameters
    println("Defining condition-specific parameters...")
    @variable(m, k1[1:3])
    @constraint(m, k1[1] == k1_free)
    @constraint(m, k1[2] == k1_free)
    @constraint(m, k1[3] == 0.1)


    # Define overrides
    # Define observable overrides
    println("Defining observableParameter overrides...")

    # Define noise overrides
    println("Defining noiseParameter overrides...")

    # Model compartments
    println("Defining compartments...")
    @variable(m, compartment == 1.0, start=1.0)

    # Model species
    println("Defining species...")
    @variable(m, 0.0 <= A[j in 1:2, k in 1:(length(t_sim)+1)] <= 11.0)
    @variable(m, 0.0 <= B[j in 1:2, k in 1:(length(t_sim)+1)] <= 11.0)

    # Model initial assignments
    println("Defining initial assignments...")
    @constraint(m, [j in 1:1], A[cond_without_preequ[j],1] == a0)
    @constraint(m, [j in 1:1], B[cond_without_preequ[j],1] == b0)

    # Model ODEs
    println("Defining ODEs...")
    println("A")
    @NLconstraint(m, [j in 1:2, k in 1:length(t_sim)-1],
        A[j, k+1] == A[j, k] + ( -1.0*( compartment * k1[j_to_cond_par[j]] * A[j, k+1] ) +1.0*( compartment * k2 * B[j, k+1] )     ) * ( t_sim[k+1] - t_sim[k] ) )
    println("B")
    @NLconstraint(m, [j in 1:2, k in 1:length(t_sim)-1],
        B[j, k+1] == B[j, k] + ( +1.0*( compartment * k1[j_to_cond_par[j]] * A[j, k+1] ) -1.0*( compartment * k2 * B[j, k+1] )     ) * ( t_sim[k+1] - t_sim[k] ) )

    # Pre-equilibration constraints
    println("Defining pre-equilibration constraints...")
    @constraint(m, [j in 1:1], A[cond_without_preequ[j], length(t_sim)+1] == 0)
    println("A")
    @NLconstraint(m, [j in 1:1],
        -1.0*( compartment * k1[preequ[j]] * A[cond_with_preequ[j], length(t_sim)+1] ) +1.0*( compartment * k2 * B[cond_with_preequ[j], length(t_sim)+1] )  == 0 )
    @constraint(m, [j in 1:1], A[cond_with_preequ[j], length(t_sim)+1] == A[cond_with_preequ[j], 1])
    @constraint(m, [j in 1:1], B[cond_without_preequ[j], length(t_sim)+1] == 0)
    println("B")
    @NLconstraint(m, [j in 1:1],
        +1.0*( compartment * k1[preequ[j]] * A[cond_with_preequ[j], length(t_sim)+1] ) -1.0*( compartment * k2 * B[cond_with_preequ[j], length(t_sim)+1] )  == 0 )
    @constraint(m, [j in 1:1], B[cond_with_preequ[j], length(t_sim)+1] == B[cond_with_preequ[j], 1])

    # Define observables
    println("Defining observables...")
    @variable(m, -0.19999999999999996 <= obs_a[j in 1:2, k in 1:length(t_exp)] <= 1.6, start=1.)
    @NLconstraint(m, [j in 1:2, k in 1:length(t_exp)], obs_a[j, k] == A[j, t_sim_to_exp[k]])
    @variable(m, -0.4 <= obs_b[j in 1:2, k in 1:length(t_exp)] <= 0.8, start=1.)
    @NLconstraint(m, [j in 1:2, k in 1:length(t_exp)], obs_b[j, k] == offset_B+scaling_B*B[j, t_sim_to_exp[k]])

    # Define objective
    println("Defining objective...")
    @NLobjective(m, Min, sum(0.5 * log(2*pi*(( noise_A1 + noise_A2 * obs_a[j, k] ))^2) + 0.5*((obs_a[j, k_to_time_idx[j][k]]-data[j][k, :obs_a])/(( noise_A1 + noise_A2 * obs_a[j, k] )))^2 for j in 1:2 for k in 1:length(t_exp))
        + sum(0.5 * log(2*pi*(( noise_B ))^2) + 0.5*((obs_b[j, k_to_time_idx[j][k]]-data[j][k, :obs_b])/(( noise_B )))^2 for j in 1:2 for k in 1:length(t_exp))
        )

    println("Optimizing:")
    optimize!(m)

    println("Transfering results to Python...")
    parameter_names = [a0, b0, k1_free, k2, noise_A1, noise_A2, noise_B, offset_B, scaling_B]
    parameter_values = Dict()
    for p in parameter_names
        if occursin("[", string(p))
            parameter_values[split(string(p[1]), "[")[1]] = JuMP.value.(p)
        else
            parameter_values[string(p)] = JuMP.value.(p)
        end
    end

    species_names = [A, B, ]
    species_values = Dict()
    for s in species_names
        species_values[split(string(s[1]), "[")[1]] = JuMP.value.(s)
    end

    observable_names = [obs_a, obs_b, ]
    observable_values = Dict()
    for o in observable_names
        observable_values[split(string(o[1]), "[")[1]] = Array(JuMP.value.(o))
    end

    objective_val = objective_value(m)

    results["objective_value"][string(i_start)] = objective_val
    results["parameters"][string(i_start)] = parameter_values
    results["species"][string(i_start)] = species_values
    results["observables"][string(i_start)] = observable_values

results