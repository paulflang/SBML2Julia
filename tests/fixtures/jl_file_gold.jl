using CSV
using DataFrames
using Ipopt
using JuMP

t_ratio = 2 # Setting number of ODE discretisation steps

# Data
data_path = "/media/sf_DPhil_Project/Project07_Parameter Fitting/df_software/DisFit/tests/fixtures/G2M_copasi/measurementData_G2M_copasi.tsv"
df = CSV.read(data_path)
df = unstack(df, :time, :observableId, :measurement)
t_exp = Vector(df[!, :time]) # set of simulation times.)
t_sim = range(0, stop=t_exp[end], length=t_exp[end]*t_ratio+1)

results = Dict()
results["objective_val"] = Dict()
results["x"] = Dict()
results["states"] = Dict()
results["observables"] = Dict()
for i_start in 1:1
    m = Model(with_optimizer(Ipopt.Optimizer))

    @variable(m, 0.05 <= kPhEnsa <= 0.2, start=0.05+(0.2-0.05)*rand(Float64))
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

    # Model states
    println("Defining states ...")
    @variable(m, 0 <= Cb[k in 1:length(t_sim)] <= 1)
    @variable(m, 0 <= pCb[k in 1:length(t_sim)] <= 1)
    @variable(m, 0 <= Wee[k in 1:length(t_sim)] <= 1)
    @variable(m, 0 <= pWee[k in 1:length(t_sim)] <= 1)
    @variable(m, 0 <= Cdc25[k in 1:length(t_sim)] <= 1)
    @variable(m, 0 <= pCdc25[k in 1:length(t_sim)] <= 1)
    @variable(m, 0 <= Gw[k in 1:length(t_sim)] <= 1)
    @variable(m, 0 <= pGw[k in 1:length(t_sim)] <= 1)
    @variable(m, 0 <= Ensa[k in 1:length(t_sim)] <= 1)
    @variable(m, 0 <= pEnsa[k in 1:length(t_sim)] <= 1)
    @variable(m, 0 <= pEB55[k in 1:length(t_sim)] <= 1)
    @variable(m, 0 <= B55[k in 1:length(t_sim)] <= 1)

    # Model ODEs
    println("Defining ODEs ...")
    @NLconstraint(m, [k in 1:length(t_sim)-1],
        Cb[k+1] == Cb[k] + ( -1.0*(kWee2 * Cb[k+1] * Wee[k+1]) -1.0*(kWee1 * Cb[k+1]) +1.0*(kCdc25_1 * pCb[k+1]) +1.0*(kCdc25_2 * pCb[k+1] * pCdc25[k+1])     ) * ( t_sim[k+1] - t_sim[k] ) )
    @NLconstraint(m, [k in 1:length(t_sim)-1],
        pCb[k+1] == pCb[k] + ( +1.0*(kWee2 * Cb[k+1] * Wee[k+1]) +1.0*(kWee1 * Cb[k+1]) -1.0*(kCdc25_1 * pCb[k+1]) -1.0*(kCdc25_2 * pCb[k+1] * pCdc25[k+1])     ) * ( t_sim[k+1] - t_sim[k] ) )
    @NLconstraint(m, [k in 1:length(t_sim)-1],
        Wee[k+1] == Wee[k] + ( -1.0*(kPhWee * Cb[k+1] * Wee[k+1]) +1.0*(kDpWee * pWee[k+1] * B55[k+1])     ) * ( t_sim[k+1] - t_sim[k] ) )
    @NLconstraint(m, [k in 1:length(t_sim)-1],
        pWee[k+1] == pWee[k] + ( +1.0*(kPhWee * Cb[k+1] * Wee[k+1]) -1.0*(kDpWee * pWee[k+1] * B55[k+1])     ) * ( t_sim[k+1] - t_sim[k] ) )
    @NLconstraint(m, [k in 1:length(t_sim)-1],
        Cdc25[k+1] == Cdc25[k] + ( -1.0*(kPhCdc25 * Cb[k+1] * Cdc25[k+1]) +1.0*(kDpCdc25 * pCdc25[k+1] * B55[k+1])     ) * ( t_sim[k+1] - t_sim[k] ) )
    @NLconstraint(m, [k in 1:length(t_sim)-1],
        pCdc25[k+1] == pCdc25[k] + ( +1.0*(kPhCdc25 * Cb[k+1] * Cdc25[k+1]) -1.0*(kDpCdc25 * pCdc25[k+1] * B55[k+1])     ) * ( t_sim[k+1] - t_sim[k] ) )
    @NLconstraint(m, [k in 1:length(t_sim)-1],
        Gw[k+1] == Gw[k] + ( -1.0*(kPhGw * Gw[k+1] * Cb[k+1]) +1.0*(kDpGw1 * pGw[k+1]) +1.0*(kDpGw2 * pGw[k+1] * B55[k+1])     ) * ( t_sim[k+1] - t_sim[k] ) )
    @NLconstraint(m, [k in 1:length(t_sim)-1],
        pGw[k+1] == pGw[k] + ( +1.0*(kPhGw * Gw[k+1] * Cb[k+1]) -1.0*(kDpGw1 * pGw[k+1]) -1.0*(kDpGw2 * pGw[k+1] * B55[k+1])     ) * ( t_sim[k+1] - t_sim[k] ) )
    @NLconstraint(m, [k in 1:length(t_sim)-1],
        Ensa[k+1] == Ensa[k] + ( -1.0*(kPhEnsa * Ensa[k+1] * pGw[k+1]) +1.0*(kDpEnsa * pEB55[k+1])     ) * ( t_sim[k+1] - t_sim[k] ) )
    @NLconstraint(m, [k in 1:length(t_sim)-1],
        pEnsa[k+1] == pEnsa[k] + ( +1.0*(kPhEnsa * Ensa[k+1] * pGw[k+1]) -1.0*(kAspEB55 * B55[k+1] * pEnsa[k+1]) +1.0*(kDipEB55 * pEB55[k+1])     ) * ( t_sim[k+1] - t_sim[k] ) )
    @NLconstraint(m, [k in 1:length(t_sim)-1],
        pEB55[k+1] == pEB55[k] + ( +1.0*(kAspEB55 * B55[k+1] * pEnsa[k+1]) -1.0*(kDipEB55 * pEB55[k+1]) -1.0*(kDpEnsa * pEB55[k+1])     ) * ( t_sim[k+1] - t_sim[k] ) )
    @NLconstraint(m, [k in 1:length(t_sim)-1],
        B55[k+1] == B55[k] + ( -1.0*(kAspEB55 * B55[k+1] * pEnsa[k+1]) +1.0*(kDipEB55 * pEB55[k+1]) +1.0*(kDpEnsa * pEB55[k+1])     ) * ( t_sim[k+1] - t_sim[k] ) )

    # Define observables
    println("Defining observables ...")
    @variable(m, -0.16183341480000007 <= obs_Cb[k in 1:length(t_sim)] <= 2.2624076128)
    @NLconstraint(m, [k in 1:length(t_sim)], obs_Cb[k] == fCb*Cb[k])
    @variable(m, 0.01396170799999999 <= obs_Gw[k in 1:length(t_sim)] <= 1.3805782530000001)
    @NLconstraint(m, [k in 1:length(t_sim)], obs_Gw[k] == 2*Gw[k]/2)
    @variable(m, -0.1368000294 <= obs_pEnsa[k in 1:length(t_sim)] <= 0.8208001764)
    @NLconstraint(m, [k in 1:length(t_sim)], obs_pEnsa[k] == 1+pEnsa[k]-1)
    @variable(m, -0.23442151780000003 <= obs_pWee[k in 1:length(t_sim)] <= 1.4065291068)
    @NLconstraint(m, [k in 1:length(t_sim)], obs_pWee[k] == pWee[k])
    @variable(m, -0.21328375020000004 <= obs_Cdc25[k in 1:length(t_sim)] <= 1.3121887522)
    @NLconstraint(m, [k in 1:length(t_sim)], obs_Cdc25[k] == Cdc25[k])
    @variable(m, -0.0530869484 <= obs_Ensa[k in 1:length(t_sim)] <= 1.1980051134)
    @NLconstraint(m, [k in 1:length(t_sim)], obs_Ensa[k] == Ensa[k])
    @variable(m, -0.18421561320000002 <= obs_pGw[k in 1:length(t_sim)] <= 1.1052936791999999)
    @NLconstraint(m, [k in 1:length(t_sim)], obs_pGw[k] == pGw[k])
    @variable(m, -0.055821097 <= obs_pEB55[k in 1:length(t_sim)] <= 0.334926582)
    @NLconstraint(m, [k in 1:length(t_sim)], obs_pEB55[k] == pEB55[k])
    @variable(m, -0.2037955172 <= obs_Wee[k in 1:length(t_sim)] <= 1.2537439451999999)
    @NLconstraint(m, [k in 1:length(t_sim)], obs_Wee[k] == Wee[k])
    @variable(m, -0.254477912 <= obs_pCdc25[k in 1:length(t_sim)] <= 1.526867472)
    @NLconstraint(m, [k in 1:length(t_sim)], obs_pCdc25[k] == pCdc25[k])
    @variable(m, -0.0489286938 <= obs_B55[k in 1:length(t_sim)] <= 0.2959637828)
    @NLconstraint(m, [k in 1:length(t_sim)], obs_B55[k] == B55[k])
    @variable(m, -0.1581079262 <= obs_pCb[k in 1:length(t_sim)] <= 0.9486475571999999)
    @NLconstraint(m, [k in 1:length(t_sim)], obs_pCb[k] == pCb[k])

    # Define objective
    println("Defining objective ...")
    @NLobjective(m, Min,sum((obs_Cb[(k-1)*t_ratio+1]-df[k, :obs_Cb])^2 for k in 1:length(t_exp))
        + sum((obs_Gw[(k-1)*t_ratio+1]-df[k, :obs_Gw])^2 for k in 1:length(t_exp))
        + sum((obs_pEnsa[(k-1)*t_ratio+1]-df[k, :obs_pEnsa])^2 for k in 1:length(t_exp))
        + sum((obs_pWee[(k-1)*t_ratio+1]-df[k, :obs_pWee])^2 for k in 1:length(t_exp))
        + sum((obs_Cdc25[(k-1)*t_ratio+1]-df[k, :obs_Cdc25])^2 for k in 1:length(t_exp))
        + sum((obs_Ensa[(k-1)*t_ratio+1]-df[k, :obs_Ensa])^2 for k in 1:length(t_exp))
        + sum((obs_pGw[(k-1)*t_ratio+1]-df[k, :obs_pGw])^2 for k in 1:length(t_exp))
        + sum((obs_pEB55[(k-1)*t_ratio+1]-df[k, :obs_pEB55])^2 for k in 1:length(t_exp))
        + sum((obs_Wee[(k-1)*t_ratio+1]-df[k, :obs_Wee])^2 for k in 1:length(t_exp))
        + sum((obs_pCdc25[(k-1)*t_ratio+1]-df[k, :obs_pCdc25])^2 for k in 1:length(t_exp))
        + sum((obs_B55[(k-1)*t_ratio+1]-df[k, :obs_B55])^2 for k in 1:length(t_exp))
        + sum((obs_pCb[(k-1)*t_ratio+1]-df[k, :obs_pCb])^2 for k in 1:length(t_exp))
        )

    println("Optimizing...")
    optimize!(m)

    println("Retreiving solution...")
    params = [kPhEnsa, kDpEnsa, kPhGw, kDpGw1, kDpGw2, kWee1, kWee2, kPhWee, kDpWee, kCdc25_1, kCdc25_2, kPhCdc25, kDpCdc25, kDipEB55, kAspEB55, fCb]
    paramvalues = Dict()
    for param in params
        paramvalues[param] = JuMP.value.(param)
    end

    variables = [Cb, pCb, Wee, pWee, Cdc25, pCdc25, Gw, pGw, Ensa, pEnsa, pEB55, B55, ]
    variablevalues = Dict()
    for v in variables
        variablevalues[string(v[1])[1:end-3]] = Vector(JuMP.value.(v))
    end

    observables = [obs_Cb, obs_Gw, obs_pEnsa, obs_pWee, obs_Cdc25, obs_Ensa, obs_pGw, obs_pEB55, obs_Wee, obs_pCdc25, obs_B55, obs_pCb, ]
    observablevalues = Dict()
    for o in observables
        observablevalues[string(o[1])[1:end-3]] = Vector(JuMP.value.(o))
    end

    v = objective_value(m)

    results["objective_val"][string(i_start)] = v
    results["x"][string(i_start)] = paramvalues
    results["states"][string(i_start)] = variablevalues
    results["observables"][string(i_start)] = observablevalues

end

results