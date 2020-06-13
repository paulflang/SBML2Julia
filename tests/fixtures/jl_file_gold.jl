using CSV
using DataFrames
using Ipopt
using JuMP

t_ratio = 2 # Setting number of ODE discretisation steps

# Data
data_path = "/media/sf_DPhil_Project/Project07_Parameter Fitting/df_software/DisFit/tests/fixtures/G2M_copasi/measurementData_G2M_copasi.tsv"
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
    m = Model(with_optimizer(Ipopt.Optimizer, tol=1e-6))

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
    @constraint(m, fB55[1] == 1.0)
    @constraint(m, fB55[2] == 0.9)
    @constraint(m, fB55[3] == 1.1)
    @constraint(m, fB55[4] == 1.0)

    @variable(m, kPhEnsa[1:4])
    @constraint(m, kPhEnsa[1] == 0.1)
    @constraint(m, kPhEnsa[2] == 0.1)
    @constraint(m, kPhEnsa[3] == 0.1)
    @constraint(m, kPhEnsa[4] == 0.09)

    # Define global parameters
    @variable(m, kDpEnsa == 0.05)
    @variable(m, kPhGw == 1.0)
    @variable(m, kDpGw1 == 0.25)
    @variable(m, kDpGw2 == 10.0)
    @variable(m, kWee1 == 0.01)
    @variable(m, kWee2 == 0.99)
    @variable(m, kPhWee == 1.0)
    @variable(m, kDpWee == 10.0)
    @variable(m, kCdc25_1 == 0.1)
    @variable(m, kCdc25_2 == 0.9)
    @variable(m, kPhCdc25 == 1.0)
    @variable(m, kDpCdc25 == 10.0)
    @variable(m, kDipEB55 == 0.0068)
    @variable(m, kAspEB55 == 57.0)
    @variable(m, fCb == 2.0)
    @variable(m, jiWee == 0.1)

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
        Cb[j, k+1] == Cb[j, k] + ( -1.0*( kWee2 * Cb[j, k+1] * 0.5 * (Wee[j, k+1] - iWee[j] - jiWee + sqrt(^(-Wee[j, k+1] + iWee[j] + jiWee, 2) + 4 * jiWee * Wee[j, k+1])) ) -1.0*( kWee1 * Cb[j, k+1] ) +1.0*( kCdc25_1 * pCb[j, k+1] ) +1.0*( kCdc25_2 * pCb[j, k+1] * pCdc25[j, k+1] )     ) * ( t_sim[k+1] - t_sim[k] ) )
    @NLconstraint(m, [j in 1:4, k in 1:length(t_sim)-1],
        pCb[j, k+1] == pCb[j, k] + ( +1.0*( kWee2 * Cb[j, k+1] * 0.5 * (Wee[j, k+1] - iWee[j] - jiWee + sqrt(^(-Wee[j, k+1] + iWee[j] + jiWee, 2) + 4 * jiWee * Wee[j, k+1])) ) +1.0*( kWee1 * Cb[j, k+1] ) -1.0*( kCdc25_1 * pCb[j, k+1] ) -1.0*( kCdc25_2 * pCb[j, k+1] * pCdc25[j, k+1] )     ) * ( t_sim[k+1] - t_sim[k] ) )
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
    @variable(m, -0.09027266 <= obs_Cb[j in 1:4, k in 1:length(t_sim)] <= 1.88171211)
    @NLconstraint(m, [j in 1:4, k in 1:length(t_sim)], obs_Cb[j, k] == fCb*Cb[j, k])
    @variable(m, 0.09218833600000001 <= obs_Gw[j in 1:4, k in 1:length(t_sim)] <= 1.151301944)
    @NLconstraint(m, [j in 1:4, k in 1:length(t_sim)], obs_Gw[j, k] == 2*Gw[j, k]/2)
    @variable(m, -0.1170945822 <= obs_pEnsa[j in 1:4, k in 1:length(t_sim)] <= 0.7025674932)
    @NLconstraint(m, [j in 1:4, k in 1:length(t_sim)], obs_pEnsa[j, k] == 1+pEnsa[j, k]-1)
    @variable(m, -0.19893025860000002 <= obs_pWee[j in 1:4, k in 1:length(t_sim)] <= 1.1935815516)
    @NLconstraint(m, [j in 1:4, k in 1:length(t_sim)], obs_pWee[j, k] == pWee[j, k])
    @variable(m, -0.1935815516 <= obs_Cdc25[j in 1:4, k in 1:length(t_sim)] <= 1.1989302586)
    @NLconstraint(m, [j in 1:4, k in 1:length(t_sim)], obs_Cdc25[j, k] == Cdc25[j, k])
    @variable(m, -0.0020577524000000014 <= obs_Ensa[j in 1:4, k in 1:length(t_sim)] <= 1.1670096254)
    @NLconstraint(m, [j in 1:4, k in 1:length(t_sim)], obs_Ensa[j, k] == Ensa[j, k])
    @variable(m, -0.15130194400000002 <= obs_pGw[j in 1:4, k in 1:length(t_sim)] <= 0.907811664)
    @NLconstraint(m, [j in 1:4, k in 1:length(t_sim)], obs_pGw[j, k] == pGw[j, k])
    @variable(m, -0.04991504320000001 <= obs_pEB55[j in 1:4, k in 1:length(t_sim)] <= 0.29949025920000005)
    @NLconstraint(m, [j in 1:4, k in 1:length(t_sim)], obs_pEB55[j, k] == pEB55[j, k])
    @variable(m, -0.1935815516 <= obs_Wee[j in 1:4, k in 1:length(t_sim)] <= 1.1989302586)
    @NLconstraint(m, [j in 1:4, k in 1:length(t_sim)], obs_Wee[j, k] == Wee[j, k])
    @variable(m, -0.19893025860000002 <= obs_pCdc25[j in 1:4, k in 1:length(t_sim)] <= 1.1935815516)
    @NLconstraint(m, [j in 1:4, k in 1:length(t_sim)], obs_pCdc25[j, k] == pCdc25[j, k])
    @variable(m, -0.054541232800000006 <= obs_B55[j in 1:4, k in 1:length(t_sim)] <= 0.3299235388)
    @NLconstraint(m, [j in 1:4, k in 1:length(t_sim)], obs_B55[j, k] == B55[j, k]*fB55[j])
    @variable(m, -0.13909814980000001 <= obs_pCb[j in 1:4, k in 1:length(t_sim)] <= 0.8345888988000001)
    @NLconstraint(m, [j in 1:4, k in 1:length(t_sim)], obs_pCb[j, k] == pCb[j, k])

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
    params = [kDpEnsa, kPhGw, kDpGw1, kDpGw2, kWee1, kWee2, kPhWee, kDpWee, kCdc25_1, kCdc25_2, kPhCdc25, kDpCdc25, kDipEB55, kAspEB55, fCb, jiWee, fB55, kPhEnsa]
    paramvalues = Dict()
    for p in params
        if occursin("[", string(p))
            paramvalues[split(string(p[1]), "[")[1]] = JuMP.value.(p)
        else
            paramvalues[string(p)] = JuMP.value.(p)
        end
    end

    variables = [Cb, pCb, Wee, pWee, Cdc25, pCdc25, Gw, pGw, Ensa, pEnsa, pEB55, B55, iWee, ]
    variablevalues = Dict()
    for v in variables
        variablevalues[split(string(v[1]), "[")[1]] = JuMP.value.(v)
    end

    observables = [obs_Cb, obs_Gw, obs_pEnsa, obs_pWee, obs_Cdc25, obs_Ensa, obs_pGw, obs_pEB55, obs_Wee, obs_pCdc25, obs_B55, obs_pCb, ]
    observablevalues = Dict()
    for o in observables
        observablevalues[split(string(o[1]), "[")[1]] = Array(JuMP.value.(o))
    end

    v = objective_value(m)

    results["objective_val"][string(i_start)] = v
    results["x"][string(i_start)] = paramvalues
    results["states"][string(i_start)] = variablevalues
    results["observables"][string(i_start)] = observablevalues

end

results