using CSV
using DataFrames
using Ipopt
using JuMP

n_steps = 101 # Setting number of ODE discretisation steps

# Data
println("Reading measurement data...")
data_path = "/media/sf_DPhil_Project/Project07_Parameter Fitting/df_software/SBML2JuliaMP/tests/fixtures/0015_objectivePrior/_measurements.tsv"
df = CSV.read(data_path)
insert!(df, 1, (1:length(df[:,1])), :id)

df_by_o = groupby(df, :observableId)
df_by_c = groupby(df, :simulationConditionId)

# Setting variables
t_sim = range(0, stop=maximum(df[!, :time]), length=n_steps)
t_exp = Dict()
m_exp = Dict()
t_sim_to_exp = Dict()
deduplicated_time_idx = Dict()

cond2idx = Dict()
for i_cond in 1:length(keys(df_by_c))
    cond_name = values(keys(df_by_c)[i_cond])[1]
    cond2idx[cond_name] = i_cond
end

for i_obs in 1:length(keys(df_by_o))
    obs_name = values(keys(df_by_o)[i_obs])[1]
    t_exp[obs_name] = Dict()
    m_exp[obs_name] = Dict()
    t_sim_to_exp[obs_name] = Dict()
    deduplicated_time_idx[obs_name] = Dict()
    df_by_o_c = groupby(DataFrame(df_by_o[i_obs]), :simulationConditionId)
    for i_cond in 1:length(keys(df_by_o_c))
        cond_name = values(keys(df_by_o_c)[i_cond])[1]
        cond_idx = cond2idx[cond_name]
        t_exp[obs_name][cond_idx] = df_by_o_c[i_cond][!, :time]
        m_exp[obs_name][cond_idx] = df_by_o_c[i_cond][!, :measurement]

        t = df_by_o_c[i_cond][!, :time]
        tmp = []
        for i in 1:length(t)
            idx = argmin(abs.(t[i] .- t_sim))
            append!(tmp, idx)
        end
        t_sim_to_exp[obs_name][cond_idx] = tmp

        deduplicated_time_idx[obs_name][cond_idx] = []
        idx = 0
        prev = "a"
        for t in df_by_o_c[i_cond][!, :time]
            if t != prev
                idx = idx + 1
            end
            push!(deduplicated_time_idx[obs_name][cond_idx], idx)
            prev = t
        end
    end
end

obs2conds = Dict()
for obs_dict in t_exp
    obs_name = obs_dict[1]
    obs2conds[obs_name] = []
    for cond_dict in obs_dict[2]
        cond_idx = cond_dict[1]
        push!(obs2conds[obs_name], cond_idx)
    end
    sort!(obs2conds[obs_name])
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

    species_dict = Dict()
    obs_dict = Dict()

    m = Model(with_optimizer(Ipopt.Optimizer))

    # Define global parameters
    println("Defining global parameters...")
    @variable(m, 0.0 <= a0 <= 10.0, start=0.0+(10.0-(0.0))*rand(Float64))
    @variable(m, 0.0 <= b0 <= 10.0, start=0.0+(10.0-(0.0))*rand(Float64))
    @variable(m, 0.0 <= k1_free <= 10.0, start=0.0+(10.0-(0.0))*rand(Float64))
    @variable(m, 0.0 <= k2 <= 10.0, start=0.0+(10.0-(0.0))*rand(Float64))
    @variable(m, 0.005 <= noise_A1 <= 0.1, start=0.005+(0.1-(0.005))*rand(Float64))
    @variable(m, 0.01 <= noise_A2 <= 0.2, start=0.01+(0.2-(0.01))*rand(Float64))
    @variable(m, 0.01 <= noise_B <= 0.2, start=0.01+(0.2-(0.01))*rand(Float64))
    @variable(m, 0.0 <= offset_B <= 5.0, start=0.0+(5.0-(0.0))*rand(Float64))
    @variable(m, 0.0 <= scaling_B <= 10.0, start=0.0+(10.0-(0.0))*rand(Float64))

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
    species_dict["A"] = @variable(m, 0.0 <= A[j in 1:2, k in 1:(length(t_sim)+1)] <= 20.038222743)
    species_dict["B"] = @variable(m, 0.0 <= B[j in 1:2, k in 1:(length(t_sim)+1)] <= 20.038222743)

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
    @constraint(m, [j in cond_without_preequ], A[j, length(t_sim)+1] == 0) # Dummy preequilibration for these conditions
    println("A")
    @NLconstraint(m, [j in 1:1],
        -1.0*( compartment * k1[preequ[j]] * A[cond_with_preequ[j], length(t_sim)+1] ) +1.0*( compartment * k2 * B[cond_with_preequ[j], length(t_sim)+1] )  == 0 )
    @constraint(m, [j in 1:1], A[cond_with_preequ[j], length(t_sim)+1] == A[cond_with_preequ[j], 1])
    @constraint(m, [j in cond_without_preequ], B[j, length(t_sim)+1] == 0) # Dummy preequilibration for these conditions
    println("B")
    @NLconstraint(m, [j in 1:1],
        +1.0*( compartment * k1[preequ[j]] * A[cond_with_preequ[j], length(t_sim)+1] ) -1.0*( compartment * k2 * B[cond_with_preequ[j], length(t_sim)+1] )  == 0 )
    @constraint(m, [j in 1:1], B[cond_with_preequ[j], length(t_sim)+1] == B[cond_with_preequ[j], 1])

    # Define observables
    println("Defining observables...")
    obs_dict["obs_a"] = @variable(m, -0.4096061150000001 <= obs_a[j in obs2conds["obs_a"], k in 1:length(t_exp["obs_a"][j])] <= 1.8169054180000002, start=1.)
    @NLconstraint(m, [j in obs2conds["obs_a"], k in 1:length(t_exp["obs_a"][j])], obs_a[j, k] == A[j, t_sim_to_exp["obs_a"][j][k]])
    obs_dict["obs_b"] = @variable(m, -0.6546633489999999 <= obs_b[j in obs2conds["obs_b"], k in 1:length(t_exp["obs_b"][j])] <= 3.059816594, start=1.)
    @NLconstraint(m, [j in obs2conds["obs_b"], k in 1:length(t_exp["obs_b"][j])], obs_b[j, k] == offset_B+scaling_B*B[j, t_sim_to_exp["obs_b"][j][k]])

    # Defining objectivePriors
    println("Defining objectivePriors")
    @variable(m, prior_mean[l in 1:4], start=1.)
    @constraint(m, prior_mean[1] == 1)
    @constraint(m, prior_mean[2] == 0)
    @constraint(m, prior_mean[3] == -1.22314355131)
    @constraint(m, prior_mean[4] == -0.51082562376)

    @variable(m, 0 <= prior_std[l in 1:4], start=1.)
    @constraint(m, prior_std[1] == 1)
    @constraint(m, prior_std[2] == 1)
    @constraint(m, prior_std[3] == 1)
    @constraint(m, prior_std[4] == 2)

    @variable(m, par_est[l in 1:4], start=1.)
    @constraint(m, par_est[1] == a0)
    @constraint(m, par_est[2] == b0)
    @constraint(m, par_est[3] == k1_free)
    @constraint(m, par_est[4] == k2)

    normal_priors = [1]
    laplace_priors = [2]
    logNormal_priors = [3]
    logLaplace_priors = [4]
    @variable(m, delta_par[l in [2]])
    @constraint(m, [l in [2]], delta_par[l] >= par_est[l] - prior_mean[l])
    @constraint(m, [l in [2]], delta_par[l] >= prior_mean[l] - par_est[l])

    @variable(m, delta_log_par[l in [4]])
    @NLconstraint(m, [l in [4]], delta_log_par[l] >= log(par_est[l]) - prior_mean[l])
    @NLconstraint(m, [l in [4]], delta_log_par[l] >= prior_mean[l] - log(par_est[l]))

    @variable(m, nlp_normal)
    @NLconstraint(m, nlp_normal == sum( 0.5 * log(2*pi*prior_std[l]^2) + 0.5*((par_est[l]-prior_mean[l])/prior_std[l])^2 for l in normal_priors))
    @variable(m, nlp_laplace)
    @NLconstraint(m, nlp_laplace == sum( log(2*prior_std[l]) + delta_par[l]/prior_std[l] for l in laplace_priors))
    @variable(m, nlp_logNormal)
    @NLconstraint(m, nlp_logNormal == sum( log(par_est[l]*prior_std[l]*sqrt(2*pi)) + 0.5*((log(par_est[l])-prior_mean[l])/prior_std[l])^2 for l in logNormal_priors))
    @variable(m, nlp_logLaplace)
    @NLconstraint(m, nlp_logLaplace == sum( log(2*prior_std[l]*par_est[l]) + delta_log_par[l]/prior_std[l] for l in logLaplace_priors))

    @variable(m, neg_log_priors)
    @constraint(m, neg_log_priors == nlp_laplace + nlp_logLaplace + nlp_logNormal + nlp_normal)

    # Define objective
    println("Defining objective...")
    @NLobjective(m, Min, sum(0.5 * log(2*pi*(( noise_A1 + noise_A2 * obs_a[j, k] ))^2) + 0.5*((obs_a[j, deduplicated_time_idx["obs_a"][j][k]]-m_exp["obs_a"][j][k])/(( noise_A1 + noise_A2 * obs_a[j, k] )))^2 for j in obs2conds["obs_a"] for k in 1:length(t_exp["obs_a"][j]))
        + sum(0.5 * log(2*pi*(( noise_B ))^2) + 0.5*((obs_b[j, deduplicated_time_idx["obs_b"][j][k]]-m_exp["obs_b"][j][k])/(( noise_B )))^2 for j in obs2conds["obs_b"] for k in 1:length(t_exp["obs_b"][j]))
        + neg_log_priors)

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

    species_names = ["A", "B", ]
    species_values = Dict(string(spec)=>Dict(cond_idx=>[value(species_dict[spec][cond_idx, val_idx]) for val_idx=1:length(t_sim)]
        for cond_idx in 1:length(keys(df_by_c))) for spec in species_names)

    observable_values = Dict(obs=>Dict(cond_idx=>[value(obs_dict[obs][cond_idx, val_idx]) for val_idx=1:length(list)]
        for (cond_idx, list) in dict) for (obs, dict) in m_exp)

    objective_val = objective_value(m)

    results["objective_value"][string(i_start)] = objective_val
    results["parameters"][string(i_start)] = parameter_values
    results["species"][string(i_start)] = species_values
    results["observables"][string(i_start)] = observable_values

results