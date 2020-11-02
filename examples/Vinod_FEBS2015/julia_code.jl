using CSV
using DataFrames
using Ipopt
using JuMP

n_steps = 121 # Setting number of ODE discretisation steps

# Data
println("Reading measurement data...")
data_path = "/media/sf_DPhil_Project/Project07_Parameter Fitting/df_software/sbml2julia/examples/Vinod_FEBS2015/measurementData_Vinod_FEBS2015.tsv"
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
end

results = Dict()
results["objective_value"] = Dict()
results["parameters"] = Dict()
results["species"] = Dict()
results["observables"] = Dict()

j_to_cond_par = [1, 2, 3, 4]
cond_without_preequ = [1, 2, 3, 4]
cond_with_preequ = []
preequ = []

i_start = 1
    m = Model(with_optimizer(Ipopt.Optimizer))

    set_optimizer_attribute(m,"linear_solver","MA57")

    # Define global parameters
    println("Defining global parameters...")
    @variable(m, 0.025 <= kDpEnsa <= 0.1, start=0.025+(0.1-(0.025))*rand(Float64))
    @variable(m, 0.5 <= kPhGw <= 2.0, start=0.5+(2.0-(0.5))*rand(Float64))
    @variable(m, 0.125 <= kDpGw1 <= 0.5, start=0.125+(0.5-(0.125))*rand(Float64))
    @variable(m, 5.0 <= kDpGw2 <= 20.0, start=5.0+(20.0-(5.0))*rand(Float64))
    @variable(m, 0.005 <= kWee1 <= 0.02, start=0.005+(0.02-(0.005))*rand(Float64))
    @variable(m, 0.495 <= kWee2 <= 1.98, start=0.495+(1.98-(0.495))*rand(Float64))
    @variable(m, 0.5 <= kPhWee <= 2.0, start=0.5+(2.0-(0.5))*rand(Float64))
    @variable(m, 5.0 <= kDpWee <= 20.0, start=5.0+(20.0-(5.0))*rand(Float64))
    @variable(m, 0.05 <= kCdc25_1 <= 0.2, start=0.05+(0.2-(0.05))*rand(Float64))
    @variable(m, 0.45 <= kCdc25_2 <= 1.8, start=0.45+(1.8-(0.45))*rand(Float64))
    @variable(m, 0.5 <= kPhCdc25 <= 2.0, start=0.5+(2.0-(0.5))*rand(Float64))
    @variable(m, 5.0 <= kDpCdc25 <= 20.0, start=5.0+(20.0-(5.0))*rand(Float64))
    @variable(m, 0.0034 <= kDipEB55 <= 0.0136, start=0.0034+(0.0136-(0.0034))*rand(Float64))
    @variable(m, 28.5 <= kAspEB55 <= 114.0, start=28.5+(114.0-(28.5))*rand(Float64))
    @variable(m, 1.0 <= fCb <= 4.0, start=1.0+(4.0-(1.0))*rand(Float64))
    @variable(m, 0.05 <= jiWee <= 0.2, start=0.05+(0.2-(0.05))*rand(Float64))
    @variable(m, 0.5 <= fB55_wt <= 2.0, start=0.5+(2.0-(0.5))*rand(Float64))
    @variable(m, 0.45 <= fB55_iWee <= 1.8, start=0.45+(1.8-(0.45))*rand(Float64))
    @variable(m, 0.55 <= fB55_Cb_low <= 2.2, start=0.55+(2.2-(0.55))*rand(Float64))
    @variable(m, 0.5 <= fB55_pGw_weak <= 2.0, start=0.5+(2.0-(0.5))*rand(Float64))
    @variable(m, 0.05 <= kPhEnsa_wt <= 0.2, start=0.05+(0.2-(0.05))*rand(Float64))
    @variable(m, 0.05 <= kPhEnsa_iWee <= 0.2, start=0.05+(0.2-(0.05))*rand(Float64))
    @variable(m, 0.05 <= kPhEnsa_Cb_low <= 0.2, start=0.05+(0.2-(0.05))*rand(Float64))
    @variable(m, 0.045 <= kPhEnsa_pGw_weak <= 0.18, start=0.045+(0.18-(0.045))*rand(Float64))

    # Define condition-specific parameters
    println("Defining condition-specific parameters...")
    @variable(m, iWee_0[1:4])
    @constraint(m, iWee_0[1] == 0.0)
    @constraint(m, iWee_0[2] == 0.1)
    @constraint(m, iWee_0[3] == 0.0)
    @constraint(m, iWee_0[4] == 0.0)

    @variable(m, Cb_0[1:4])
    @constraint(m, Cb_0[1] == 0.8)
    @constraint(m, Cb_0[2] == 0.8)
    @constraint(m, Cb_0[3] == 0.75)
    @constraint(m, Cb_0[4] == 0.8)

    @variable(m, fB55[1:4])
    @constraint(m, fB55[1] == fB55_wt)
    @constraint(m, fB55[2] == fB55_iWee)
    @constraint(m, fB55[3] == fB55_Cb_low)
    @constraint(m, fB55[4] == fB55_pGw_weak)

    @variable(m, kPhEnsa[1:4])
    @constraint(m, kPhEnsa[1] == kPhEnsa_wt)
    @constraint(m, kPhEnsa[2] == kPhEnsa_iWee)
    @constraint(m, kPhEnsa[3] == kPhEnsa_Cb_low)
    @constraint(m, kPhEnsa[4] == kPhEnsa_pGw_weak)


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
    @variable(m, 0 <= Cb[j in 1:4, k in 1:(length(t_sim)+1)] <= 1.1)
    @variable(m, 0 <= pCb[j in 1:4, k in 1:(length(t_sim)+1)] <= 1.1)
    @variable(m, 0 <= Wee[j in 1:4, k in 1:(length(t_sim)+1)] <= 1.1)
    @variable(m, 0 <= pWee[j in 1:4, k in 1:(length(t_sim)+1)] <= 1.1)
    @variable(m, 0 <= Cdc25[j in 1:4, k in 1:(length(t_sim)+1)] <= 1.1)
    @variable(m, 0 <= pCdc25[j in 1:4, k in 1:(length(t_sim)+1)] <= 1.1)
    @variable(m, 0 <= Gw[j in 1:4, k in 1:(length(t_sim)+1)] <= 1.1)
    @variable(m, 0 <= pGw[j in 1:4, k in 1:(length(t_sim)+1)] <= 1.1)
    @variable(m, 0 <= Ensa[j in 1:4, k in 1:(length(t_sim)+1)] <= 1.1)
    @variable(m, 0 <= pEnsa[j in 1:4, k in 1:(length(t_sim)+1)] <= 1.1)
    @variable(m, 0 <= pEB55[j in 1:4, k in 1:(length(t_sim)+1)] <= 1.1)
    @variable(m, 0 <= B55[j in 1:4, k in 1:(length(t_sim)+1)] <= 1.1)
    @variable(m, iWee[j in 1:4, k in 1:(length(t_sim)+1)])

    # Model initial assignments
    println("Defining initial assignments...")
    @constraint(m, [j in 1:4], Cb[cond_without_preequ[j],1] == Cb_0[cond_without_preequ[j]])
    @constraint(m, [j in 1:4], iWee[cond_without_preequ[j],1] == iWee_0[cond_without_preequ[j]])

    # Model ODEs
    println("Defining ODEs...")
    println("Cb")
    @NLconstraint(m, [j in 1:4, k in 1:length(t_sim)-1],
        Cb[j, k+1] == Cb[j, k] + ( -1.0*( compartment * kWee2 * Cb[j, k+1] * 0.5 * (Wee[j, k+1] - iWee[j] - jiWee + sqrt(^(-Wee[j, k+1] + iWee[j] + jiWee, 2) + 4 * jiWee * Wee[j, k+1])) ) -1.0*( compartment * kWee1 * Cb[j, k+1] ) +1.0*( compartment * kCdc25_1 * pCb[j, k+1] ) +1.0*( compartment * kCdc25_2 * pCb[j, k+1] * pCdc25[j, k+1] )     ) * ( t_sim[k+1] - t_sim[k] ) )
    println("pCb")
    @NLconstraint(m, [j in 1:4, k in 1:length(t_sim)-1],
        pCb[j, k+1] == pCb[j, k] + ( +1.0*( compartment * kWee2 * Cb[j, k+1] * 0.5 * (Wee[j, k+1] - iWee[j] - jiWee + sqrt(^(-Wee[j, k+1] + iWee[j] + jiWee, 2) + 4 * jiWee * Wee[j, k+1])) ) +1.0*( compartment * kWee1 * Cb[j, k+1] ) -1.0*( compartment * kCdc25_1 * pCb[j, k+1] ) -1.0*( compartment * kCdc25_2 * pCb[j, k+1] * pCdc25[j, k+1] )     ) * ( t_sim[k+1] - t_sim[k] ) )
    println("Wee")
    @NLconstraint(m, [j in 1:4, k in 1:length(t_sim)-1],
        Wee[j, k+1] == Wee[j, k] + ( -1.0*( compartment * kPhWee * Cb[j, k+1] * Wee[j, k+1] ) +1.0*( compartment * kDpWee * pWee[j, k+1] * B55[j, k+1] )     ) * ( t_sim[k+1] - t_sim[k] ) )
    println("pWee")
    @NLconstraint(m, [j in 1:4, k in 1:length(t_sim)-1],
        pWee[j, k+1] == pWee[j, k] + ( +1.0*( compartment * kPhWee * Cb[j, k+1] * Wee[j, k+1] ) -1.0*( compartment * kDpWee * pWee[j, k+1] * B55[j, k+1] )     ) * ( t_sim[k+1] - t_sim[k] ) )
    println("Cdc25")
    @NLconstraint(m, [j in 1:4, k in 1:length(t_sim)-1],
        Cdc25[j, k+1] == Cdc25[j, k] + ( -1.0*( compartment * kPhCdc25 * Cb[j, k+1] * Cdc25[j, k+1] ) +1.0*( compartment * kDpCdc25 * pCdc25[j, k+1] * B55[j, k+1] )     ) * ( t_sim[k+1] - t_sim[k] ) )
    println("pCdc25")
    @NLconstraint(m, [j in 1:4, k in 1:length(t_sim)-1],
        pCdc25[j, k+1] == pCdc25[j, k] + ( +1.0*( compartment * kPhCdc25 * Cb[j, k+1] * Cdc25[j, k+1] ) -1.0*( compartment * kDpCdc25 * pCdc25[j, k+1] * B55[j, k+1] )     ) * ( t_sim[k+1] - t_sim[k] ) )
    println("Gw")
    @NLconstraint(m, [j in 1:4, k in 1:length(t_sim)-1],
        Gw[j, k+1] == Gw[j, k] + ( -1.0*( compartment * kPhGw * Gw[j, k+1] * Cb[j, k+1] ) +1.0*( compartment * kDpGw1 * pGw[j, k+1] ) +1.0*( compartment * kDpGw2 * pGw[j, k+1] * B55[j, k+1] )     ) * ( t_sim[k+1] - t_sim[k] ) )
    println("pGw")
    @NLconstraint(m, [j in 1:4, k in 1:length(t_sim)-1],
        pGw[j, k+1] == pGw[j, k] + ( +1.0*( compartment * kPhGw * Gw[j, k+1] * Cb[j, k+1] ) -1.0*( compartment * kDpGw1 * pGw[j, k+1] ) -1.0*( compartment * kDpGw2 * pGw[j, k+1] * B55[j, k+1] )     ) * ( t_sim[k+1] - t_sim[k] ) )
    println("Ensa")
    @NLconstraint(m, [j in 1:4, k in 1:length(t_sim)-1],
        Ensa[j, k+1] == Ensa[j, k] + ( -1.0*( compartment * kPhEnsa[j_to_cond_par[j]] * Ensa[j, k+1] * pGw[j, k+1] ) +1.0*( compartment * kDpEnsa * pEB55[j, k+1] )     ) * ( t_sim[k+1] - t_sim[k] ) )
    println("pEnsa")
    @NLconstraint(m, [j in 1:4, k in 1:length(t_sim)-1],
        pEnsa[j, k+1] == pEnsa[j, k] + ( +1.0*( compartment * kPhEnsa[j_to_cond_par[j]] * Ensa[j, k+1] * pGw[j, k+1] ) -1.0*( compartment * kAspEB55 * B55[j, k+1] * pEnsa[j, k+1] ) +1.0*( compartment * kDipEB55 * pEB55[j, k+1] )     ) * ( t_sim[k+1] - t_sim[k] ) )
    println("pEB55")
    @NLconstraint(m, [j in 1:4, k in 1:length(t_sim)-1],
        pEB55[j, k+1] == pEB55[j, k] + ( +1.0*( compartment * kAspEB55 * B55[j, k+1] * pEnsa[j, k+1] ) -1.0*( compartment * kDipEB55 * pEB55[j, k+1] ) -1.0*( compartment * kDpEnsa * pEB55[j, k+1] )     ) * ( t_sim[k+1] - t_sim[k] ) )
    println("B55")
    @NLconstraint(m, [j in 1:4, k in 1:length(t_sim)-1],
        B55[j, k+1] == B55[j, k] + ( -1.0*( compartment * kAspEB55 * B55[j, k+1] * pEnsa[j, k+1] ) +1.0*( compartment * kDipEB55 * pEB55[j, k+1] ) +1.0*( compartment * kDpEnsa * pEB55[j, k+1] )     ) * ( t_sim[k+1] - t_sim[k] ) )
    @constraint(m, [j in 1:4, k in 1:length(t_sim)], iWee[j, k] == iWee_0[j])

    # Pre-equilibration constraints
    @constraint(m, [j in 1:4], Cb[j, length(t_sim)+1] == 0) # Dummy preequilibration
    @constraint(m, [j in 1:4], pCb[j, length(t_sim)+1] == 0) # Dummy preequilibration
    @constraint(m, [j in 1:4], Wee[j, length(t_sim)+1] == 0) # Dummy preequilibration
    @constraint(m, [j in 1:4], pWee[j, length(t_sim)+1] == 0) # Dummy preequilibration
    @constraint(m, [j in 1:4], Cdc25[j, length(t_sim)+1] == 0) # Dummy preequilibration
    @constraint(m, [j in 1:4], pCdc25[j, length(t_sim)+1] == 0) # Dummy preequilibration
    @constraint(m, [j in 1:4], Gw[j, length(t_sim)+1] == 0) # Dummy preequilibration
    @constraint(m, [j in 1:4], pGw[j, length(t_sim)+1] == 0) # Dummy preequilibration
    @constraint(m, [j in 1:4], Ensa[j, length(t_sim)+1] == 0) # Dummy preequilibration
    @constraint(m, [j in 1:4], pEnsa[j, length(t_sim)+1] == 0) # Dummy preequilibration
    @constraint(m, [j in 1:4], pEB55[j, length(t_sim)+1] == 0) # Dummy preequilibration
    @constraint(m, [j in 1:4], B55[j, length(t_sim)+1] == 0) # Dummy preequilibration
    @constraint(m, [j in 1:4], iWee[j, length(t_sim)+1] == 0) # Dummy preequilibration

    # Define observables
    println("Defining observables...")
    @variable(m, -1.2171211 <= obs_Cb[j in obs2conds["obs_Cb"], k in 1:length(t_exp["obs_Cb"][j])] <= 3.0085605500000003, start=1.)
    @NLconstraint(m, [j in obs2conds["obs_Cb"], k in 1:length(t_exp["obs_Cb"][j])], obs_Cb[j, k] == fCb*Cb[j, t_sim_to_exp["obs_Cb"][j][k]])
    @variable(m, -0.5130194399999999 <= obs_Gw[j in obs2conds["obs_Gw"], k in 1:length(t_exp["obs_Gw"][j])] <= 1.75650972, start=1.)
    @NLconstraint(m, [j in obs2conds["obs_Gw"], k in 1:length(t_exp["obs_Gw"][j])], obs_Gw[j, k] == 2*Gw[j, t_sim_to_exp["obs_Gw"][j][k]]/2)
    @variable(m, -0.585472911 <= obs_pEnsa[j in obs2conds["obs_pEnsa"], k in 1:length(t_exp["obs_pEnsa"][j])] <= 1.170945822, start=1.)
    @NLconstraint(m, [j in obs2conds["obs_pEnsa"], k in 1:length(t_exp["obs_pEnsa"][j])], obs_pEnsa[j, k] == 1+pEnsa[j, t_sim_to_exp["obs_pEnsa"][j][k]]-1)
    @variable(m, -0.994651293 <= obs_pWee[j in obs2conds["obs_pWee"], k in 1:length(t_exp["obs_pWee"][j])] <= 1.989302586, start=1.)
    @NLconstraint(m, [j in obs2conds["obs_pWee"], k in 1:length(t_exp["obs_pWee"][j])], obs_pWee[j, k] == pWee[j, t_sim_to_exp["obs_pWee"][j][k]])
    @variable(m, -0.989302586 <= obs_Cdc25[j in obs2conds["obs_Cdc25"], k in 1:length(t_exp["obs_Cdc25"][j])] <= 1.994651293, start=1.)
    @NLconstraint(m, [j in obs2conds["obs_Cdc25"], k in 1:length(t_exp["obs_Cdc25"][j])], obs_Cdc25[j, k] == Cdc25[j, t_sim_to_exp["obs_Cdc25"][j][k]])
    @variable(m, -0.670096254 <= obs_Ensa[j in obs2conds["obs_Ensa"], k in 1:length(t_exp["obs_Ensa"][j])] <= 1.8350481269999999, start=1.)
    @NLconstraint(m, [j in obs2conds["obs_Ensa"], k in 1:length(t_exp["obs_Ensa"][j])], obs_Ensa[j, k] == Ensa[j, t_sim_to_exp["obs_Ensa"][j][k]])
    @variable(m, -0.75650972 <= obs_pGw[j in obs2conds["obs_pGw"], k in 1:length(t_exp["obs_pGw"][j])] <= 1.51301944, start=1.)
    @NLconstraint(m, [j in obs2conds["obs_pGw"], k in 1:length(t_exp["obs_pGw"][j])], obs_pGw[j, k] == pGw[j, t_sim_to_exp["obs_pGw"][j][k]])
    @variable(m, -0.249575216 <= obs_pEB55[j in obs2conds["obs_pEB55"], k in 1:length(t_exp["obs_pEB55"][j])] <= 0.499150432, start=1.)
    @NLconstraint(m, [j in obs2conds["obs_pEB55"], k in 1:length(t_exp["obs_pEB55"][j])], obs_pEB55[j, k] == pEB55[j, t_sim_to_exp["obs_pEB55"][j][k]])
    @variable(m, -0.989302586 <= obs_Wee[j in obs2conds["obs_Wee"], k in 1:length(t_exp["obs_Wee"][j])] <= 1.994651293, start=1.)
    @NLconstraint(m, [j in obs2conds["obs_Wee"], k in 1:length(t_exp["obs_Wee"][j])], obs_Wee[j, k] == Wee[j, t_sim_to_exp["obs_Wee"][j][k]])
    @variable(m, -0.994651293 <= obs_pCdc25[j in obs2conds["obs_pCdc25"], k in 1:length(t_exp["obs_pCdc25"][j])] <= 1.989302586, start=1.)
    @NLconstraint(m, [j in obs2conds["obs_pCdc25"], k in 1:length(t_exp["obs_pCdc25"][j])], obs_pCdc25[j, k] == pCdc25[j, t_sim_to_exp["obs_pCdc25"][j][k]])
    @variable(m, -0.274235388 <= obs_B55[j in obs2conds["obs_B55"], k in 1:length(t_exp["obs_B55"][j])] <= 0.549617694, start=1.)
    @NLconstraint(m, [j in obs2conds["obs_B55"], k in 1:length(t_exp["obs_B55"][j])], obs_B55[j, k] == B55[j, t_sim_to_exp["obs_B55"][j][k]]*fB55[j])
    @variable(m, -0.695490749 <= obs_pCb[j in obs2conds["obs_pCb"], k in 1:length(t_exp["obs_pCb"][j])] <= 1.390981498, start=1.)
    @NLconstraint(m, [j in obs2conds["obs_pCb"], k in 1:length(t_exp["obs_pCb"][j])], obs_pCb[j, k] == pCb[j, t_sim_to_exp["obs_pCb"][j][k]])

    # Define objective
    println("Defining objective...")
    @NLobjective(m, Min, sum(0.5 * log(2*pi*(( 1 ))^2) + 0.5*((obs_Cb[j, deduplicated_time_idx["obs_Cb"][j][k]]-m_exp["obs_Cb"][j][k])/(( 1 )))^2 for j in obs2conds["obs_Cb"] for k in 1:length(t_exp["obs_Cb"][j]))
        + sum(0.5 * log(2*pi*(( 1 ))^2) + 0.5*((obs_Gw[j, deduplicated_time_idx["obs_Gw"][j][k]]-m_exp["obs_Gw"][j][k])/(( 1 )))^2 for j in obs2conds["obs_Gw"] for k in 1:length(t_exp["obs_Gw"][j]))
        + sum(0.5 * log(2*pi*(( 1 ))^2) + 0.5*((obs_pEnsa[j, deduplicated_time_idx["obs_pEnsa"][j][k]]-m_exp["obs_pEnsa"][j][k])/(( 1 )))^2 for j in obs2conds["obs_pEnsa"] for k in 1:length(t_exp["obs_pEnsa"][j]))
        + sum(0.5 * log(2*pi*(( 1 ))^2) + 0.5*((obs_pWee[j, deduplicated_time_idx["obs_pWee"][j][k]]-m_exp["obs_pWee"][j][k])/(( 1 )))^2 for j in obs2conds["obs_pWee"] for k in 1:length(t_exp["obs_pWee"][j]))
        + sum(0.5 * log(2*pi*(( 1 ))^2) + 0.5*((obs_Cdc25[j, deduplicated_time_idx["obs_Cdc25"][j][k]]-m_exp["obs_Cdc25"][j][k])/(( 1 )))^2 for j in obs2conds["obs_Cdc25"] for k in 1:length(t_exp["obs_Cdc25"][j]))
        + sum(0.5 * log(2*pi*(( 1 ))^2) + 0.5*((obs_Ensa[j, deduplicated_time_idx["obs_Ensa"][j][k]]-m_exp["obs_Ensa"][j][k])/(( 1 )))^2 for j in obs2conds["obs_Ensa"] for k in 1:length(t_exp["obs_Ensa"][j]))
        + sum(0.5 * log(2*pi*(( 1 ))^2) + 0.5*((obs_pGw[j, deduplicated_time_idx["obs_pGw"][j][k]]-m_exp["obs_pGw"][j][k])/(( 1 )))^2 for j in obs2conds["obs_pGw"] for k in 1:length(t_exp["obs_pGw"][j]))
        + sum(0.5 * log(2*pi*(( 1 ))^2) + 0.5*((obs_pEB55[j, deduplicated_time_idx["obs_pEB55"][j][k]]-m_exp["obs_pEB55"][j][k])/(( 1 )))^2 for j in obs2conds["obs_pEB55"] for k in 1:length(t_exp["obs_pEB55"][j]))
        + sum(0.5 * log(2*pi*(( 1 ))^2) + 0.5*((obs_Wee[j, deduplicated_time_idx["obs_Wee"][j][k]]-m_exp["obs_Wee"][j][k])/(( 1 )))^2 for j in obs2conds["obs_Wee"] for k in 1:length(t_exp["obs_Wee"][j]))
        + sum(0.5 * log(2*pi*(( 1 ))^2) + 0.5*((obs_pCdc25[j, deduplicated_time_idx["obs_pCdc25"][j][k]]-m_exp["obs_pCdc25"][j][k])/(( 1 )))^2 for j in obs2conds["obs_pCdc25"] for k in 1:length(t_exp["obs_pCdc25"][j]))
        + sum(0.5 * log(2*pi*(( 1 ))^2) + 0.5*((obs_B55[j, deduplicated_time_idx["obs_B55"][j][k]]-m_exp["obs_B55"][j][k])/(( 1 )))^2 for j in obs2conds["obs_B55"] for k in 1:length(t_exp["obs_B55"][j]))
        + sum(0.5 * log(2*pi*(( 1 ))^2) + 0.5*((obs_pCb[j, deduplicated_time_idx["obs_pCb"][j][k]]-m_exp["obs_pCb"][j][k])/(( 1 )))^2 for j in obs2conds["obs_pCb"] for k in 1:length(t_exp["obs_pCb"][j]))
        )

    println("Optimizing:")
    optimize!(m)

    println("Transfering results to Python...")
    parameter_names = [kDpEnsa, kPhGw, kDpGw1, kDpGw2, kWee1, kWee2, kPhWee, kDpWee, kCdc25_1, kCdc25_2, kPhCdc25, kDpCdc25, kDipEB55, kAspEB55, fCb, jiWee, fB55_wt, fB55_iWee, fB55_Cb_low, fB55_pGw_weak, kPhEnsa_wt, kPhEnsa_iWee, kPhEnsa_Cb_low, kPhEnsa_pGw_weak]
    parameter_values = Dict()
    for p in parameter_names
        if occursin("[", string(p))
            parameter_values[split(string(p[1]), "[")[1]] = JuMP.value.(p)
        else
            parameter_values[string(p)] = JuMP.value.(p)
        end
    end

    species_names = ["Cb", "pCb", "Wee", "pWee", "Cdc25", "pCdc25", "Gw", "pGw", "Ensa", "pEnsa", "pEB55", "B55", "iWee", ]
    species_values = Dict(string(spec)=>Dict(cond_idx=>[value(eval(Symbol(spec))[cond_idx, val_idx]) for val_idx=1:length(t_sim)]
        for cond_idx in 1:length(keys(df_by_c))) for spec in species_names)

    observable_values = Dict(obs=>Dict(cond_idx=>[value(eval(Symbol(obs))[cond_idx, val_idx]) for val_idx=1:length(list)]
        for (cond_idx, list) in dict) for (obs, dict) in m_exp)

    objective_val = objective_value(m)

    results["objective_value"][string(i_start)] = objective_val
    results["parameters"][string(i_start)] = parameter_values
    results["species"][string(i_start)] = species_values
    results["observables"][string(i_start)] = observable_values

results