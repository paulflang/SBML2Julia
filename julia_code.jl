using CSV
using DataFrames
using Ipopt
using JuMP

t_ratio = 2 # Setting number of ODE discretisation steps

# Data
data_path = "tests/fixtures/G2M_copasi/measurementData_G2M_copasi.tsv"
df = CSV.read(data_path)
dfg = groupby(df, :simulationConditionId)
data = []
for condition in keys(dfg)
    push!(data,unstack(dfg[condition], :time, :observableId, :measurement))
end

t_exp = Vector(DataFrame(groupby(dfg[1], :observableId)[1])[!, :time])
t_sim = range(0, stop=t_exp[end], length=t_exp[end]*t_ratio+1)

results = Dict()
results["objective_val"] = Dict()
results["x"] = Dict()
results["states"] = Dict()
results["observables"] = Dict()
for i_start in 1:1
    m = Model(with_optimizer(Ipopt.Optimizer,tol=1e-6))

    # Define condition-defined parameters
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

    # Define condition-local parameters
    @variable(m, fB55[1:4])
    @constraint(m, 0.025 <= fB55[1] <= 0.1)
    @constraint(m, 0.5 <= fB55[2] <= 2.0)
    @constraint(m, 0.125 <= fB55[3] <= 0.5)
    @constraint(m, 5.0 <= fB55[4] <= 20.0)

    @variable(m, kPhEnsa[1:4])
    @constraint(m, 0.025 <= kPhEnsa[1] <= 0.1)
    @constraint(m, 0.5 <= kPhEnsa[2] <= 2.0)
    @constraint(m, 0.125 <= kPhEnsa[3] <= 0.5)
    @constraint(m, 5.0 <= kPhEnsa[4] <= 20.0)

    # Define global parameters
    @variable(m, 0.025 <= kDpEnsa <= 0.1, start=0.025+(0.1-0.025)*rand(Float64))
    @variable(m, 0.5 <= kPhGw <= 2.0, start=0.5+(2.0-0.5)*rand(Float64))
    @variable(m, 0.125 <= kDpGw1 <= 0.5, start=0.125+(0.5-0.125)*rand(Float64))
    @variable(m, 5.0 <= kDpGw2 <= 20.0, start=5.0+(20.0-5.0)*rand(Float64))
    @variable(m, 0.005 <= kWee1 <= 0.02, start=0.005+(0.02-0.005)*rand(Float64))
    @variable(m, 0.495 <= kWee2 <= 1.98, start=0.495+(1.98-0.495)*rand(Float64))
    @variable(m, 0.5 <= kPhWee <= 2.0, start=0.5+(2.0-0.5)*rand(Float64))
    @variable(m, 5.0 <= kDpWee <= 20.0, start=5.0+(20.0-5.0)*rand(Float64))
    @variable(m, 0.05 <= kCdc25_1 <= 0.2, start=0.05+(0.2-0.05)*rand(Float64))
    @variable(m, 0.45 <= kCdc25_2 <= 1.8, start=0.45+(1.8-0.45)*rand(Float64))
    @variable(m, 0.5 <= kPhCdc25 <= 2.0, start=0.5+(2.0-0.5)*rand(Float64))
    @variable(m, 5.0 <= kDpCdc25 <= 20.0, start=5.0+(20.0-5.0)*rand(Float64))
    @variable(m, 0.0034 <= kDipEB55 <= 0.0136, start=0.0034+(0.0136-0.0034)*rand(Float64))
    @variable(m, 28.5 <= kAspEB55 <= 114.0, start=28.5+(114.0-28.5)*rand(Float64))
    @variable(m, 1.0 <= fCb <= 4.0, start=1.0+(4.0-1.0)*rand(Float64))
    @variable(m, 0.05 <= jiWee <= 0.2, start=0.05+(0.2-0.05)*rand(Float64))

    # Model states
    println("Defining states ...")
    @variable(m, 0 <= Cb[j in 1:4, k in 1:length(t_sim)] <= 1.1)
    @variable(m, 0 <= pCb[j in 1:4, k in 1:length(t_sim)] <= 1.1)
    @variable(m, 0 <= Wee[j in 1:4, k in 1:length(t_sim)] <= 1.1)
    @variable(m, 0 <= pWee[j in 1:4, k in 1:length(t_sim)] <= 1.1)
    @variable(m, 0 <= Cdc25[j in 1:4, k in 1:length(t_sim)] <= 1.1)
    @variable(m, 0 <= pCdc25[j in 1:4, k in 1:length(t_sim)] <= 1.1)
    @variable(m, 0 <= Gw[j in 1:4, k in 1:length(t_sim)] <= 1.1)
    @variable(m, 0 <= pGw[j in 1:4, k in 1:length(t_sim)] <= 1.1)
    @variable(m, 0 <= Ensa[j in 1:4, k in 1:length(t_sim)] <= 1.1)
    @variable(m, 0 <= pEnsa[j in 1:4, k in 1:length(t_sim)] <= 1.1)
    @variable(m, 0 <= pEB55[j in 1:4, k in 1:length(t_sim)] <= 1.1)
    @variable(m, 0 <= B55[j in 1:4, k in 1:length(t_sim)] <= 1.1)
    @variable(m, iWee[j in 1:4])

    # Model ODEs
    println("Defining ODEs ...")
    @NLconstraint(m, [j in 1:4, k in 1:length(t_sim)-1],
        Cb[j, k+1] == Cb[j, k] + ( -1.0*( kWee2 * Cb[j, k+1] * 0.5 *  Wee[j, k+1] - iWee[j] - jiWee + sqrt(*( Wee[j, k+1] + iWee[j] + jiWee, 2) + 4 * jiWee * Wee[j, k+1] ) ) -1.0*( kWee1 * Cb[j, k+1] ) +1.0*( kCdc25_1 * pCb[j, k+1] ) +1.0*( kCdc25_2 * pCb[j, k+1] * pCdc25[j, k+1] )     ) * ( t_sim[k+1] - t_sim[k] ) )
    @NLconstraint(m, [j in 1:4, k in 1:length(t_sim)-1],
        pCb[j, k+1] == pCb[j, k] + ( +1.0*( kWee2 * Cb[j, k+1] * 0.5 *  Wee[j, k+1] - iWee[j] - jiWee + sqrt(*( Wee[j, k+1] + iWee[j] + jiWee, 2) + 4 * jiWee * Wee[j, k+1] ) ) +1.0*( kWee1 * Cb[j, k+1] ) -1.0*( kCdc25_1 * pCb[j, k+1] ) -1.0*( kCdc25_2 * pCb[j, k+1] * pCdc25[j, k+1] )     ) * ( t_sim[k+1] - t_sim[k] ) )
    @NLconstraint(m, [j in 1:4, k in 1:length(t_sim)-1],
        Wee[j, k+1] == Wee[j, k] + ( -1.0*( kPhWee * Cb[j, k+1] * Wee[j, k+1] ) +1.0*( kDpWee * pWee[j, k+1] * B55[j, k+1] )     ) * ( t_sim[k+1] - t_sim[k] ) )
    @NLconstraint(m, [j in 1:4, k in 1:length(t_sim)-1],
        pWee[j, k+1] == pWee[j, k] + ( +1.0*( kPhWee * Cb[j, k+1] * Wee[j, k+1] ) -1.0*( kDpWee * pWee[j, k+1] * B55[j, k+1] )     ) * ( t_sim[k+1] - t_sim[k] ) )
    @NLconstraint(m, [j in 1:4, k in 1:length(t_sim)-1],
        Cdc25[j, k+1] == Cdc25[j, k] + ( -1.0*( kPhCdc25 * Cb[j, k+1] * Cdc25[j, k+1] ) +1.0*( kDpCdc25 * pCdc25[j, k+1] * B55[j, k+1] )     ) * ( t_sim[k+1] - t_sim[k] ) )
    @NLconstraint(m, [j in 1:4, k in 1:length(t_sim)-1],
        pCdc25[j, k+1] == pCdc25[j, k] + ( +1.0*( kPhCdc25 * Cb[j, k+1] * Cdc25[j, k+1] ) -1.0*( kDpCdc25 * pCdc25[j, k+1] * B55[j, k+1] )     ) * ( t_sim[k+1] - t_sim[k] ) )
    @NLconstraint(m, [j in 1:4, k in 1:length(t_sim)-1],
        Gw[j, k+1] == Gw[j, k] + ( -1.0*( kPhGw * Gw[j, k+1] * Cb[j, k+1] ) +1.0*( kDpGw1 * pGw[j, k+1] ) +1.0*( kDpGw2 * pGw[j, k+1] * B55[j, k+1] )     ) * ( t_sim[k+1] - t_sim[k] ) )
    @NLconstraint(m, [j in 1:4, k in 1:length(t_sim)-1],
        pGw[j, k+1] == pGw[j, k] + ( +1.0*( kPhGw * Gw[j, k+1] * Cb[j, k+1] ) -1.0*( kDpGw1 * pGw[j, k+1] ) -1.0*( kDpGw2 * pGw[j, k+1] * B55[j, k+1] )     ) * ( t_sim[k+1] - t_sim[k] ) )
    @NLconstraint(m, [j in 1:4, k in 1:length(t_sim)-1],
        Ensa[j, k+1] == Ensa[j, k] + ( -1.0*( kPhEnsa[j] * Ensa[j, k+1] * pGw[j, k+1] ) +1.0*( kDpEnsa * pEB55[j, k+1] )     ) * ( t_sim[k+1] - t_sim[k] ) )
    @NLconstraint(m, [j in 1:4, k in 1:length(t_sim)-1],
        pEnsa[j, k+1] == pEnsa[j, k] + ( +1.0*( kPhEnsa[j] * Ensa[j, k+1] * pGw[j, k+1] ) -1.0*( kAspEB55 * B55[j, k+1] * pEnsa[j, k+1] ) +1.0*( kDipEB55 * pEB55[j, k+1] )     ) * ( t_sim[k+1] - t_sim[k] ) )
    @NLconstraint(m, [j in 1:4, k in 1:length(t_sim)-1],
        pEB55[j, k+1] == pEB55[j, k] + ( +1.0*( kAspEB55 * B55[j, k+1] * pEnsa[j, k+1] ) -1.0*( kDipEB55 * pEB55[j, k+1] ) -1.0*( kDpEnsa * pEB55[j, k+1] )     ) * ( t_sim[k+1] - t_sim[k] ) )
    @NLconstraint(m, [j in 1:4, k in 1:length(t_sim)-1],
        B55[j, k+1] == B55[j, k] + ( -1.0*( kAspEB55 * B55[j, k+1] * pEnsa[j, k+1] ) +1.0*( kDipEB55 * pEB55[j, k+1] ) +1.0*( kDpEnsa * pEB55[j, k+1] )     ) * ( t_sim[k+1] - t_sim[k] ) )
    @constraint(m, [j in 1:4], iWee[j] == iWee_0[j])

    # Define observables
    println("Defining observables ...")
    @variable(m, -0.04513633 <= obs_Cb[j in 1:4, k in 1:length(t_sim)] <= 0.940856055)
    @NLconstraint(m, [j in 1:4, k in 1:length(t_sim)], obs_Cb[j,k] == fCb*Cb[j,k])
    @variable(m, 0.09218833600000001 <= obs_Gw[j in 1:4, k in 1:length(t_sim)] <= 1.151301944)
    @NLconstraint(m, [j in 1:4, k in 1:length(t_sim)], obs_Gw[j,k] == 2*Gw[j,k]/2)
    @variable(m, -0.1170945822 <= obs_pEnsa[j in 1:4, k in 1:length(t_sim)] <= 0.7025674932)
    @NLconstraint(m, [j in 1:4, k in 1:length(t_sim)], obs_pEnsa[j,k] == 1+pEnsa[j,k]-1)
    @variable(m, -0.19893025860000002 <= obs_pWee[j in 1:4, k in 1:length(t_sim)] <= 1.1935815516)
    @NLconstraint(m, [j in 1:4, k in 1:length(t_sim)], obs_pWee[j,k] == pWee[j,k])
    @variable(m, -0.1935815516 <= obs_Cdc25[j in 1:4, k in 1:length(t_sim)] <= 1.1989302586)
    @NLconstraint(m, [j in 1:4, k in 1:length(t_sim)], obs_Cdc25[j,k] == Cdc25[j,k])
    @variable(m, -0.0020577524000000014 <= obs_Ensa[j in 1:4, k in 1:length(t_sim)] <= 1.1670096254)
    @NLconstraint(m, [j in 1:4, k in 1:length(t_sim)], obs_Ensa[j,k] == Ensa[j,k])
    @variable(m, -0.15130194400000002 <= obs_pGw[j in 1:4, k in 1:length(t_sim)] <= 0.907811664)
    @NLconstraint(m, [j in 1:4, k in 1:length(t_sim)], obs_pGw[j,k] == pGw[j,k])
    @variable(m, -0.04991504320000001 <= obs_pEB55[j in 1:4, k in 1:length(t_sim)] <= 0.29949025920000005)
    @NLconstraint(m, [j in 1:4, k in 1:length(t_sim)], obs_pEB55[j,k] == pEB55[j,k])
    @variable(m, -0.1935815516 <= obs_Wee[j in 1:4, k in 1:length(t_sim)] <= 1.1989302586)
    @NLconstraint(m, [j in 1:4, k in 1:length(t_sim)], obs_Wee[j,k] == Wee[j,k])
    @variable(m, -0.19893025860000002 <= obs_pCdc25[j in 1:4, k in 1:length(t_sim)] <= 1.1935815516)
    @NLconstraint(m, [j in 1:4, k in 1:length(t_sim)], obs_pCdc25[j,k] == pCdc25[j,k])
    @variable(m, -0.049490259200000004 <= obs_B55[j in 1:4, k in 1:length(t_sim)] <= 0.2999150432)
    @NLconstraint(m, [j in 1:4, k in 1:length(t_sim)], obs_B55[j,k] == B55[j,k])
    @variable(m, -0.13909814980000001 <= obs_pCb[j in 1:4, k in 1:length(t_sim)] <= 0.8345888988000001)
    @NLconstraint(m, [j in 1:4, k in 1:length(t_sim)], obs_pCb[j,k] == pCb[j,k])

    # Define objective
    println("Defining objective ...")
    @NLobjective(m, Min,sum((obs_Cb[j, (k-1)*t_ratio+1]-data[j][k, :obs_Cb])^2 for j in 1:4 for k in 1:length(t_exp))
        + sum((obs_Gw[j, (k-1)*t_ratio+1]-data[j][k, :obs_Gw])^2 for j in 1:4 for k in 1:length(t_exp))
        + sum((obs_pEnsa[j, (k-1)*t_ratio+1]-data[j][k, :obs_pEnsa])^2 for j in 1:4 for k in 1:length(t_exp))
        + sum((obs_pWee[j, (k-1)*t_ratio+1]-data[j][k, :obs_pWee])^2 for j in 1:4 for k in 1:length(t_exp))
        + sum((obs_Cdc25[j, (k-1)*t_ratio+1]-data[j][k, :obs_Cdc25])^2 for j in 1:4 for k in 1:length(t_exp))
        + sum((obs_Ensa[j, (k-1)*t_ratio+1]-data[j][k, :obs_Ensa])^2 for j in 1:4 for k in 1:length(t_exp))
        + sum((obs_pGw[j, (k-1)*t_ratio+1]-data[j][k, :obs_pGw])^2 for j in 1:4 for k in 1:length(t_exp))
        + sum((obs_pEB55[j, (k-1)*t_ratio+1]-data[j][k, :obs_pEB55])^2 for j in 1:4 for k in 1:length(t_exp))
        + sum((obs_Wee[j, (k-1)*t_ratio+1]-data[j][k, :obs_Wee])^2 for j in 1:4 for k in 1:length(t_exp))
        + sum((obs_pCdc25[j, (k-1)*t_ratio+1]-data[j][k, :obs_pCdc25])^2 for j in 1:4 for k in 1:length(t_exp))
        + sum((obs_B55[j, (k-1)*t_ratio+1]-data[j][k, :obs_B55])^2 for j in 1:4 for k in 1:length(t_exp))
        + sum((obs_pCb[j, (k-1)*t_ratio+1]-data[j][k, :obs_pCb])^2 for j in 1:4 for k in 1:length(t_exp))
        )

    println("Optimizing...")
    optimize!(m)

    println("Retreiving solution...")
    params = [kDpEnsa, kPhGw, kDpGw1, kDpGw2, kWee1, kWee2, kPhWee, kDpWee, kCdc25_1, kCdc25_2, kPhCdc25, kDpCdc25, kDipEB55, kAspEB55, fCb, jiWee, fB55] #, fB55_wt[1], fB55_iWee, fB55_Cb_low, fB55_pGw_weak, kPhEnsa_wt, kPhEnsa_iWee, kPhEnsa_Cb_low, kPhEnsa_pGw_weak]
    paramvalues = Dict()
    for p in params
        print(p)
        if occursin("[", string(p))
            paramvalues[split(string(p[1]), "[")[1]] = JuMP.value.(p)
        else
            paramvalues[string(p)] = JuMP.value.(p)
        end
    end

    println("variables")
    variables = [Cb, pCb, Wee, pWee, Cdc25, pCdc25, Gw, pGw, Ensa, pEnsa, pEB55, B55, iWee]
    variablevalues = Dict()
    println("entering for loop")
    for v in variables
        print("\n\n\n1\n")
        print(v)
        print("\n\n\n2\n")
        print(JuMP.value.(v))
        print("\n\n\n3\n")
        print(string(v[1]))
        print("\n\n\n4\n")
        print(split(string(v[1]), "[")[1])
        print("\n\n\n5\n")
        variablevalues[split(string(v[1]), "[")[1]] = JuMP.value.(v)
    end

    println("observables")
    observables = [obs_Cb, obs_Gw, obs_pEnsa, obs_pWee, obs_Cdc25, obs_Ensa, obs_pGw, obs_pEB55, obs_Wee, obs_pCdc25, obs_B55, obs_pCb]
    observablevalues = Dict()
    for o in observables
        observablevalues[split(string(o[1]), "[")[1]] = Array(JuMP.value.(o))
    end

    println("objective")
    v = objective_value(m)

    println("results")
    results["objective_val"][string(i_start)] = v
    results["x"][string(i_start)] = paramvalues
    results["states"][string(i_start)] = variablevalues
    results["observables"][string(i_start)] = observablevalues

end

results
