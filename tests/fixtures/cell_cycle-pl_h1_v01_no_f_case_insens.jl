using JuMP
using Ipopt
using CSV
using KNITRO

fc = 100 # Setting parameter search span
t_ratio = 1 # Setting ODE discretisation step
pow = 2
weight = 10

# Data
data_file_path = joinpath(pwd(), "tests", "fixtures", "G2M_copasi.csv")
df = CSV.read(data_file_path)
t_exp = Vector(df[!, :t]) # set of simulation times.)
t_sim = range(0, stop=t_exp[end], length=t_exp[end]*t_ratio+1)

# create JuMP model object")
# m = Model(with_optimizer(KNITRO.Optimizer, ms_enable=1, ms_maxtime_real=10))
m = Model(with_optimizer(Ipopt.Optimizer))

# Model parameters
println("Defining parameters ...")
@variable(m, 0.005544005544/fc <= kPhRbD <= 0.005544005544*fc)
@variable(m, 0.00693000693/fc <= kPhRbE <= 0.00693000693*fc)
@variable(m, 0.001386001386/fc <= kDpRb <= 0.001386001386*fc)
@variable(m, 0.02079002079/fc <= kPhRb <= 0.02079002079*fc)
@variable(m, 1.386001386e-05/fc <= kSyCe1 <= 1.386001386e-05*fc)
@variable(m, 0.0004158004158/fc <= kSyCe2 <= 0.0004158004158*fc)
@variable(m, 0.0001386001386/fc <= kDeCe <= 0.0001386001386*fc)
@variable(m, 0.0003465003465/fc <= kDiE2fRb <= 0.0003465003465*fc)
@variable(m, 0.3465003465/fc <= kAsE2fRb <= 0.3465003465*fc)
@variable(m, 6.93000693e-06/fc <= kSyE2f1 <= 6.93000693e-06*fc)
@variable(m, 6.93000693e-05/fc <= kSyE2f2 <= 6.93000693e-05*fc)
@variable(m, 6.93000693e-05/fc <= kDeE2f1 <= 6.93000693e-05*fc)
@variable(m, 0.001386001386/fc <= kDeE2f2 <= 0.001386001386*fc)
@variable(m, 0.001386001386/fc <= kPhE2f <= 0.001386001386*fc)
@variable(m, 0.0/fc <= kDpE2f1 <= 0.0*fc)
@variable(m, 0.0693000693/fc <= kDpE2f2 <= 0.0693000693*fc)
@variable(m, 3.465003465e-05/fc <= kSyEmi1 <= 3.465003465e-05*fc)
@variable(m, 0.0003465003465/fc <= kSyEmi2 <= 0.0003465003465*fc)
@variable(m, 0.0003465003465/fc <= kDeEmi1 <= 0.0003465003465*fc)
@variable(m, 0.0003465003465/fc <= kDeEmi2 <= 0.0003465003465*fc)
@variable(m, 0.003465003465/fc <= kDeEmi3 <= 0.003465003465*fc)
@variable(m, 3.465003465e-05/fc <= kPhEmiA <= 3.465003465e-05*fc)
@variable(m, 0.00017325017325/fc <= kPhEmiB <= 0.00017325017325*fc)
@variable(m, 3.465003465e-06/fc <= kDpEmi <= 3.465003465e-06*fc)
@variable(m, 6.93000693e-06/fc <= kSyCa1 <= 6.93000693e-06*fc)
@variable(m, 5.1975051975e-05/fc <= kSyCa2 <= 5.1975051975e-05*fc)
@variable(m, 6.93000693e-06/fc <= kDeCa1 <= 6.93000693e-06*fc)
@variable(m, 0.002772002772/fc <= kDeCa2 <= 0.002772002772*fc)
@variable(m, 0.000693000693/fc <= kDeCa3 <= 0.000693000693*fc)
@variable(m, 0.003465003465/fc <= kDpCdh <= 0.003465003465*fc)
@variable(m, 0.01732501732/fc <= kPhCdhA <= 0.01732501732*fc)
@variable(m, 0.000693000693/fc <= kPhCdhE <= 0.000693000693*fc)
@variable(m, 0.1732501733/fc <= kPhCdhB <= 0.1732501733*fc)
@variable(m, 0.0003465003465/fc <= kDiACE <= 0.0003465003465*fc)
@variable(m, 0.3465003465/fc <= kAsACE <= 0.3465003465*fc)
@variable(m, 0.0003465003465/fc <= kDiACdh <= 0.0003465003465*fc)
@variable(m, 0.3465003465/fc <= kAsACdh <= 0.3465003465*fc)
@variable(m, 0.00693000693/fc <= kPhCeA <= 0.00693000693*fc)
@variable(m, 6.93000693e-05/fc <= kPhCeE <= 6.93000693e-05*fc)
@variable(m, 0.000693000693/fc <= kPhFoxE <= 0.000693000693*fc)
@variable(m, 0.000693000693/fc <= kPhFoxA <= 0.000693000693*fc)
@variable(m, 0.003465003465/fc <= kPhFoxB <= 0.003465003465*fc)
@variable(m, 3.465003465e-05/fc <= kDpFox <= 3.465003465e-05*fc)
@variable(m, 6.93000693e-06/fc <= kSyCb1 <= 6.93000693e-06*fc)
@variable(m, 5.1975051975e-05/fc <= kSyCb2 <= 5.1975051975e-05*fc)
@variable(m, 6.93000693e-06/fc <= kDeCb1 <= 6.93000693e-06*fc)
@variable(m, 0.00086625086625/fc <= kDeCb2 <= 0.00086625086625*fc)
@variable(m, 0.00086625086625/fc <= kDeCb3 <= 0.00086625086625*fc)
@variable(m, 0.000693000693/fc <= kPhEnsa <= 0.000693000693*fc)
@variable(m, 0.0003465003465/fc <= kDpEnsa <= 0.0003465003465*fc)
@variable(m, 0.00693000693/fc <= kPhGw <= 0.00693000693*fc)
@variable(m, 0.0017325017325/fc <= kDpGw1 <= 0.0017325017325*fc)
@variable(m, 0.0693000693/fc <= kDpGw2 <= 0.0693000693*fc)
@variable(m, 6.93000693e-05/fc <= kWee1 <= 6.93000693e-05*fc)
@variable(m, 0.00693000693/fc <= kWee2 <= 0.00693000693*fc)
@variable(m, 0.000693000693/fc <= kPhWeeA <= 0.000693000693*fc)
@variable(m, 0.00693000693/fc <= kPhWeeB <= 0.00693000693*fc)
@variable(m, 0.0693000693/fc <= kDpWee <= 0.0693000693*fc)
@variable(m, 0.000693000693/fc <= kCdc25_1 <= 0.000693000693*fc)
@variable(m, 0.00693000693/fc <= kCdc25_2 <= 0.00693000693*fc)
@variable(m, 0.0/fc <= kPhC25A <= 0.0*fc)
@variable(m, 0.0693000693/fc <= kDpCdc25 <= 0.0693000693*fc)
@variable(m, 4.7124047124e-05/fc <= kDipEB55 <= 4.7124047124e-05*fc)
@variable(m, 0.39501039501/fc <= kAspEB55 <= 0.39501039501*fc)
@variable(m, 6.93000693e-06/fc <= kSyCdc_1 <= 6.93000693e-06*fc)
@variable(m, 6.93000693e-05/fc <= kSyCdc_2 <= 6.93000693e-05*fc)
@variable(m, 2.772002772e-05/fc <= kDeCdc_1 <= 2.772002772e-05*fc)
@variable(m, 0.0002772002772/fc <= kDeCdc_2 <= 0.0002772002772*fc)
@variable(m, 0.00086625086625/fc <= kDipACdc <= 0.00086625086625*fc)
@variable(m, 0.00693000693/fc <= kAspACdc <= 0.00693000693*fc)
@variable(m, 6.93000693e-05/fc <= kPhApcA <= 6.93000693e-05*fc)
@variable(m, 6.93000693e-05/fc <= kPhApcB <= 6.93000693e-05*fc)
@variable(m, 0.03465003465/fc <= kDpApc <= 0.03465003465*fc)
@variable(m, 6.93000693e-06/fc <= kSyFox1 <= 6.93000693e-06*fc)
@variable(m, 5.1975051975e-05/fc <= kSyFox2 <= 5.1975051975e-05*fc)
@variable(m, 6.93000693e-06/fc <= kDeFox1 <= 6.93000693e-06*fc)
@variable(m, 0.003465003465/fc <= kDeFox2 <= 0.003465003465*fc)
@variable(m, 0.3465003465/fc <= kAsEPx <= 0.3465003465*fc)
@variable(m, 0.03465003465/fc <= kDiEPx <= 0.03465003465*fc)
@variable(m, 0.3465003465/fc <= kAsFPcb <= 0.3465003465*fc)
@variable(m, 0.0693000693/fc <= kDiFPcb <= 0.0693000693*fc)
@variable(m, 0.3465003465/fc <= kAsFPcdc <= 0.3465003465*fc)
@variable(m, 0.17325017325/fc <= kDiFPcdc <= 0.17325017325*fc)
@variable(m, 0.00693000693/fc <= kPhC25B <= 0.00693000693*fc)
@variable(m, 0.00055440055/fc <= kPhCdhEf4 <= 0.00055440055*fc)
@variable(m, 0.01386001385/fc <= kPhCdhAf4 <= 0.01386001385*fc)
@variable(m, 0.13860013864/fc <= kPhCdhBf4 <= 0.13860013864*fc)

# Model states
println("Defining states ...")
@variable(m, 0 <= Ce[k in 1:length(t_sim)] <= 1)
@variable(m, 0 <= E2f[k in 1:length(t_sim)] <= 1)
@variable(m, 0 <= pE2f[k in 1:length(t_sim)] <= 1)
@variable(m, 0 <= Ca[k in 1:length(t_sim)] <= 1)
@variable(m, 0 <= Cb[k in 1:length(t_sim)] <= 1)
@variable(m, 0 <= E2fRb[k in 1:length(t_sim)] <= 1)
@variable(m, 0 <= pE2fRb[k in 1:length(t_sim)] <= 1)
@variable(m, 0 <= Rb[k in 1:length(t_sim)] <= 1)
@variable(m, 0 <= pRb[k in 1:length(t_sim)] <= 1)
@variable(m, 0 <= B55[k in 1:length(t_sim)] <= 1)
@variable(m, 0 <= E2fPx[k in 1:length(t_sim)] <= 1)
@variable(m, 0 <= Px[k in 1:length(t_sim)] <= 1)
@variable(m, 0 <= Apc[k in 1:length(t_sim)] <= 1)
@variable(m, 0 <= Cdc20[k in 1:length(t_sim)] <= 1)
@variable(m, 0 <= pApcCdc[k in 1:length(t_sim)] <= 1)
@variable(m, 0 <= pApc[k in 1:length(t_sim)] <= 1)
@variable(m, 0 <= ApcCdh[k in 1:length(t_sim)] <= 1)
@variable(m, 0 <= Cdh[k in 1:length(t_sim)] <= 1)
@variable(m, 0 <= pCdh[k in 1:length(t_sim)] <= 1)
@variable(m, 0 <= ACE[k in 1:length(t_sim)] <= 1)
@variable(m, 0 <= Emi[k in 1:length(t_sim)] <= 1)
@variable(m, 0 <= pACE[k in 1:length(t_sim)] <= 1)
@variable(m, 0 <= pApcCdh[k in 1:length(t_sim)] <= 1)
@variable(m, 0 <= pEmi[k in 1:length(t_sim)] <= 1)
@variable(m, 0 <= Fox[k in 1:length(t_sim)] <= 1)
@variable(m, 0 <= pFox[k in 1:length(t_sim)] <= 1)
@variable(m, 0 <= pFoxPcb[k in 1:length(t_sim)] <= 1)
@variable(m, 0 <= Procb[k in 1:length(t_sim)] <= 1)
@variable(m, 0 <= Pcdc[k in 1:length(t_sim)] <= 1)
@variable(m, 0 <= pFoxPcdc[k in 1:length(t_sim)] <= 1)
@variable(m, 0 <= Wee[k in 1:length(t_sim)] <= 1)
@variable(m, 0 <= pCb[k in 1:length(t_sim)] <= 1)
@variable(m, 0 <= pWee[k in 1:length(t_sim)] <= 1)
@variable(m, 0 <= pCdc25[k in 1:length(t_sim)] <= 1)
@variable(m, 0 <= Cdc25[k in 1:length(t_sim)] <= 1)
@variable(m, 0 <= Gw[k in 1:length(t_sim)] <= 1)
@variable(m, 0 <= pGw[k in 1:length(t_sim)] <= 1)
@variable(m, 0 <= Ensa[k in 1:length(t_sim)] <= 1)
@variable(m, 0 <= pEnsa[k in 1:length(t_sim)] <= 1)
@variable(m, 0 <= pEB55[k in 1:length(t_sim)] <= 1)

# Model ODEs
println("Defining ODEs ...")
@NLconstraint(m, [k in 1:length(t_sim)-1],
    Ce[k+1] == Ce[k] + ( +1.0*(kSyCe1) -1.0*(kDeCe * Ce[k+1]) -1.0*(kPhCeA * Ce[k+1] * Ca[k+1]) +1.0*(kSyCe2 * E2fPx[k+1]) ) * ( t_sim[k+1] - t_sim[k] ) )
@NLconstraint(m, [k in 1:length(t_sim)-1],
    E2f[k+1] == E2f[k] + ( +1.0*(kSyE2f1) +1.0*(kSyE2f2 * E2fPx[k+1]) -1.0*(kDeE2f1 * E2f[k+1]) +1.0*(kPhRbD * E2fRb[k+1] * Cd[k+1]) +1.0*(kPhRbE * E2fRb[k+1] * Ce[k+1]) +1.0*(kPhRb * E2fRb[k+1] * Ca[k+1]) +1.0*(kPhRb * E2fRb[k+1] * Cb[k+1]) -1.0*(kPhE2f * E2f[k+1] * Ca[k+1]) -1.0*(kPhE2f * E2f[k+1] * Cb[k+1]) +1.0*(kDpE2f1 * pE2f[k+1]) +1.0*(kDpE2f2 * pE2f[k+1] * B55[k+1]) -1.0*(kAsE2fRb * E2f[k+1] * Rb[k+1]) +1.0*(kDiE2fRb * E2fRb[k+1]) ) * ( t_sim[k+1] - t_sim[k] ) )
@NLconstraint(m, [k in 1:length(t_sim)-1],
    pE2f[k+1] == pE2f[k] + ( +1.0*(kPhE2f * E2f[k+1] * Ca[k+1]) +1.0*(kPhE2f * E2f[k+1] * Cb[k+1]) -1.0*(kDpE2f1 * pE2f[k+1]) -1.0*(kDeE2f1 * pE2f[k+1]) -1.0*(kDeE2f2 * pE2f[k+1]) +1.0*(kPhRbD * pE2fRb[k+1] * Cd[k+1]) +1.0*(kPhRb * pE2fRb[k+1] * Ca[k+1]) +1.0*(kPhRb * pE2fRb[k+1] * Cb[k+1]) +1.0*(kPhRbE * pE2fRb[k+1] * Ce[k+1]) -1.0*(kDpE2f2 * pE2f[k+1] * B55[k+1]) -1.0*(kAsE2fRb * pE2f[k+1] * Rb[k+1]) +1.0*(kDiE2fRb * pE2fRb[k+1]) ) * ( t_sim[k+1] - t_sim[k] ) )
@NLconstraint(m, [k in 1:length(t_sim)-1],
    Ca[k+1] == Ca[k] + ( +1.0*(kSyCa2 * E2fPx[k+1]) +1.0*(kSyCa1) -1.0*(kDeCa1 * Ca[k+1]) -1.0*(kDeCa2 * Ca[k+1] * ApcCdh[k+1]) -1.0*(kDeCa2 * Ca[k+1] * pApcCdh[k+1]) -1.0*(kDeCa3 * Ca[k+1] * pApcCdc[k+1]) ) * ( t_sim[k+1] - t_sim[k] ) )
@NLconstraint(m, [k in 1:length(t_sim)-1],
    Cb[k+1] == Cb[k] + ( +1.0*(kSyCb1) +1.0*(kSyCb2 * pFoxPcb[k+1]) -1.0*(kDeCb1 * Cb[k+1]) -1.0*(kDeCb2 * Cb[k+1] * ApcCdh[k+1]) -1.0*(kDeCb2 * Cb[k+1] * pApcCdh[k+1]) -1.0*(kDeCb3 * Cb[k+1] * pApcCdc[k+1]) -1.0*(kWee2 * Cb[k+1] * Wee[k+1]) -1.0*(kWee1 * Cb[k+1] * pWee[k+1]) +1.0*(kCdc25_2 * pCb[k+1] * pCdc25[k+1]) +1.0*(kCdc25_1 * pCb[k+1] * Cdc25[k+1]) ) * ( t_sim[k+1] - t_sim[k] ) )
@NLconstraint(m, [k in 1:length(t_sim)-1],
    E2fRb[k+1] == E2fRb[k] + ( -1.0*(kPhRbD * E2fRb[k+1] * Cd[k+1]) -1.0*(kPhRbE * E2fRb[k+1] * Ce[k+1]) -1.0*(kPhRb * E2fRb[k+1] * Ca[k+1]) -1.0*(kPhRb * E2fRb[k+1] * Cb[k+1]) +1.0*(kDpE2f2 * pE2fRb[k+1] * B55[k+1]) -1.0*(kPhE2f * E2fRb[k+1] * Ca[k+1]) -1.0*(kPhE2f * E2fRb[k+1] * Cb[k+1]) +1.0*(kDpE2f1 * pE2fRb[k+1]) -1.0*(kDeE2f1 * E2fRb[k+1]) +1.0*(kAsE2fRb * E2f[k+1] * Rb[k+1]) -1.0*(kDiE2fRb * E2fRb[k+1]) ) * ( t_sim[k+1] - t_sim[k] ) )
@NLconstraint(m, [k in 1:length(t_sim)-1],
    pE2fRb[k+1] == pE2fRb[k] + ( -1.0*(kDpE2f2 * pE2fRb[k+1] * B55[k+1]) -1.0*(kPhRbD * pE2fRb[k+1] * Cd[k+1]) -1.0*(kPhRb * pE2fRb[k+1] * Ca[k+1]) -1.0*(kPhRb * pE2fRb[k+1] * Cb[k+1]) -1.0*(kPhRbE * pE2fRb[k+1] * Ce[k+1]) +1.0*(kPhE2f * E2fRb[k+1] * Ca[k+1]) +1.0*(kPhE2f * E2fRb[k+1] * Cb[k+1]) -1.0*(kDpE2f1 * pE2fRb[k+1]) -1.0*(kDeE2f1 * pE2fRb[k+1]) -1.0*(kDeE2f2 * pE2fRb[k+1]) +1.0*(kAsE2fRb * pE2f[k+1] * Rb[k+1]) -1.0*(kDiE2fRb * pE2fRb[k+1]) ) * ( t_sim[k+1] - t_sim[k] ) )
@NLconstraint(m, [k in 1:length(t_sim)-1],
    Rb[k+1] == Rb[k] + ( +1.0*(kDeE2f1 * E2fRb[k+1]) +1.0*(kDeE2f1 * pE2fRb[k+1]) +1.0*(kDeE2f2 * pE2fRb[k+1]) -1.0*(kPhRbD * Rb[k+1] * Cd[k+1]) -1.0*(kPhRbE * Rb[k+1] * Ce[k+1]) -1.0*(kPhRb * Rb[k+1] * Ca[k+1]) -1.0*(kPhRb * Rb[k+1] * Cb[k+1]) +1.0*(kDpRb * pRb[k+1]) -1.0*(kAsE2fRb * E2f[k+1] * Rb[k+1]) +1.0*(kDiE2fRb * E2fRb[k+1]) -1.0*(kAsE2fRb * pE2f[k+1] * Rb[k+1]) +1.0*(kDiE2fRb * pE2fRb[k+1]) ) * ( t_sim[k+1] - t_sim[k] ) )
@NLconstraint(m, [k in 1:length(t_sim)-1],
    pRb[k+1] == pRb[k] + ( +1.0*(kPhRbD * E2fRb[k+1] * Cd[k+1]) +1.0*(kPhRbE * E2fRb[k+1] * Ce[k+1]) +1.0*(kPhRb * E2fRb[k+1] * Ca[k+1]) +1.0*(kPhRb * E2fRb[k+1] * Cb[k+1]) +1.0*(kPhRbD * pE2fRb[k+1] * Cd[k+1]) +1.0*(kPhRb * pE2fRb[k+1] * Ca[k+1]) +1.0*(kPhRb * pE2fRb[k+1] * Cb[k+1]) +1.0*(kPhRbE * pE2fRb[k+1] * Ce[k+1]) +1.0*(kPhRbD * Rb[k+1] * Cd[k+1]) +1.0*(kPhRbE * Rb[k+1] * Ce[k+1]) +1.0*(kPhRb * Rb[k+1] * Ca[k+1]) +1.0*(kPhRb * Rb[k+1] * Cb[k+1]) -1.0*(kDpRb * pRb[k+1]) ) * ( t_sim[k+1] - t_sim[k] ) )
@NLconstraint(m, [k in 1:length(t_sim)-1],
    B55[k+1] == B55[k] + ( +1.0*(kDpEnsa * pEB55[k+1]) -1.0*(kAspEB55 * pEnsa[k+1] * B55[k+1]) +1.0*(kDipEB55 * pEB55[k+1]) ) * ( t_sim[k+1] - t_sim[k] ) )
@NLconstraint(m, [k in 1:length(t_sim)-1],
    E2fPx[k+1] == E2fPx[k] + ( -1.0*(kDiEPx * E2fPx[k+1]) +1.0*(kAsEPx * E2f[k+1] * Px[k+1]) ) * ( t_sim[k+1] - t_sim[k] ) )
@NLconstraint(m, [k in 1:length(t_sim)-1],
    Px[k+1] == Px[k] + ( +1.0*(kDiEPx * E2fPx[k+1]) -1.0*(kAsEPx * E2f[k+1] * Px[k+1]) ) * ( t_sim[k+1] - t_sim[k] ) )
@NLconstraint(m, [k in 1:length(t_sim)-1],
    Apc[k+1] == Apc[k] + ( +1.0*(kDpApc * pApcCdc[k+1] * B55[k+1]) -1.0*(kPhApcA * Apc[k+1] * Ca[k+1]) -1.0*(kPhApcB * Apc[k+1] * Cb[k+1]) +1.0*(kDpApc * pApc[k+1] * B55[k+1]) +1.0*(kPhCdhA * ApcCdh[k+1] * Ca[k+1]) +1.0*(kPhCdhEf4 * ACE[k+1] * Ce[k+1]) +1.0*(kPhCdhAf4 * ACE[k+1] * Ca[k+1]) +1.0*(kPhCdhBf4 * ACE[k+1] * Cb[k+1]) +1.0*(kPhCdhE * ApcCdh[k+1] * Ce[k+1]) +1.0*(kPhCdhB * ApcCdh[k+1] * Cb[k+1]) -1.0*(kAsACdh * Cdh[k+1] * Apc[k+1]) +1.0*(kDiACdh * ApcCdh[k+1]) ) * ( t_sim[k+1] - t_sim[k] ) )
@NLconstraint(m, [k in 1:length(t_sim)-1],
    Cdc20[k+1] == Cdc20[k] + ( +1.0*(kDpApc * pApcCdc[k+1] * B55[k+1]) +1.0*(kSyCdc_1) +1.0*(kSyCdc_2 * pFoxPcdc[k+1]) -1.0*(kDeCdc_1 * Cdc20[k+1]) -1.0*(kDeCdc_2 * Cdc20[k+1] * ApcCdh[k+1]) -1.0*(kDeCdc_2 * Cdc20[k+1] * pApcCdh[k+1]) -1.0*(kAspACdc * pApc[k+1] * Cdc20[k+1]) +1.0*(kDipACdc * pApcCdc[k+1]) ) * ( t_sim[k+1] - t_sim[k] ) )
@NLconstraint(m, [k in 1:length(t_sim)-1],
    pApcCdc[k+1] == pApcCdc[k] + ( -1.0*(kDpApc * pApcCdc[k+1] * B55[k+1]) -1.0*(kDeCdc_1 * pApcCdc[k+1]) -1.0*(kDeCdc_2 * pApcCdc[k+1] * ApcCdh[k+1]) -1.0*(kDeCdc_2 * pApcCdc[k+1] * pApcCdh[k+1]) +1.0*(kAspACdc * pApc[k+1] * Cdc20[k+1]) -1.0*(kDipACdc * pApcCdc[k+1]) ) * ( t_sim[k+1] - t_sim[k] ) )
@NLconstraint(m, [k in 1:length(t_sim)-1],
    pApc[k+1] == pApc[k] + ( +1.0*(kPhApcA * Apc[k+1] * Ca[k+1]) +1.0*(kPhApcB * Apc[k+1] * Cb[k+1]) -1.0*(kDpApc * pApc[k+1] * B55[k+1]) +1.0*(kPhCdhEf4 * pACE[k+1] * Ce[k+1]) +1.0*(kPhCdhAf4 * pACE[k+1] * Ca[k+1]) +1.0*(kPhCdhBf4 * pACE[k+1] * Cb[k+1]) +1.0*(kPhCdhE * pApcCdh[k+1] * Ce[k+1]) +1.0*(kDeCdc_1 * pApcCdc[k+1]) +1.0*(kDeCdc_2 * pApcCdc[k+1] * ApcCdh[k+1]) +1.0*(kDeCdc_2 * pApcCdc[k+1] * pApcCdh[k+1]) +1.0*(kPhCdhA * pApcCdh[k+1] * Ca[k+1]) +1.0*(kPhCdhB * pApcCdh[k+1] * Cb[k+1]) -1.0*(kAsACdh * pApc[k+1] * Cdh[k+1]) +1.0*(kDiACdh * pApcCdh[k+1]) -1.0*(kAspACdc * pApc[k+1] * Cdc20[k+1]) +1.0*(kDipACdc * pApcCdc[k+1]) ) * ( t_sim[k+1] - t_sim[k] ) )
@NLconstraint(m, [k in 1:length(t_sim)-1],
    ApcCdh[k+1] == ApcCdh[k] + ( -1.0*(kPhCdhA * ApcCdh[k+1] * Ca[k+1]) +1.0*(kDpApc * pApcCdh[k+1] * B55[k+1]) +1.0*(kPhEmiA * ACE[k+1] * Ca[k+1]) +1.0*(kPhEmiB * ACE[k+1] * Cb[k+1]) -1.0*(kPhCdhE * ApcCdh[k+1] * Ce[k+1]) -1.0*(kPhCdhB * ApcCdh[k+1] * Cb[k+1]) +1.0*(kDeEmi1 * ACE[k+1]) +1.0*(kDeEmi2 * ACE[k+1]) -1.0*(kPhApcA * ApcCdh[k+1] * Ca[k+1]) -1.0*(kPhApcB * ApcCdh[k+1] * Cb[k+1]) +1.0*(kAsACdh * Cdh[k+1] * Apc[k+1]) -1.0*(kDiACdh * ApcCdh[k+1]) -1.0*(kAsACE * ApcCdh[k+1] * Emi[k+1]) +1.0*(kDiACE * ACE[k+1]) ) * ( t_sim[k+1] - t_sim[k] ) )
@NLconstraint(m, [k in 1:length(t_sim)-1],
    Cdh[k+1] == Cdh[k] + ( -1.0*(kPhCdhE * Cdh[k+1] * Ce[k+1]) -1.0*(kPhCdhA * Cdh[k+1] * Ca[k+1]) -1.0*(kPhCdhB * Cdh[k+1] * Cb[k+1]) +1.0*(kDpCdh * pCdh[k+1]) -1.0*(kAsACdh * Cdh[k+1] * Apc[k+1]) +1.0*(kDiACdh * ApcCdh[k+1]) -1.0*(kAsACdh * pApc[k+1] * Cdh[k+1]) +1.0*(kDiACdh * pApcCdh[k+1]) ) * ( t_sim[k+1] - t_sim[k] ) )
@NLconstraint(m, [k in 1:length(t_sim)-1],
    pCdh[k+1] == pCdh[k] + ( +1.0*(kPhCdhA * ApcCdh[k+1] * Ca[k+1]) +1.0*(kPhCdhEf4 * ACE[k+1] * Ce[k+1]) +1.0*(kPhCdhAf4 * ACE[k+1] * Ca[k+1]) +1.0*(kPhCdhBf4 * ACE[k+1] * Cb[k+1]) +1.0*(kPhCdhEf4 * pACE[k+1] * Ce[k+1]) +1.0*(kPhCdhAf4 * pACE[k+1] * Ca[k+1]) +1.0*(kPhCdhBf4 * pACE[k+1] * Cb[k+1]) +1.0*(kPhCdhE * pApcCdh[k+1] * Ce[k+1]) +1.0*(kPhCdhE * Cdh[k+1] * Ce[k+1]) +1.0*(kPhCdhA * Cdh[k+1] * Ca[k+1]) +1.0*(kPhCdhB * Cdh[k+1] * Cb[k+1]) -1.0*(kDpCdh * pCdh[k+1]) +1.0*(kPhCdhA * pApcCdh[k+1] * Ca[k+1]) +1.0*(kPhCdhB * pApcCdh[k+1] * Cb[k+1]) +1.0*(kPhCdhE * ApcCdh[k+1] * Ce[k+1]) +1.0*(kPhCdhB * ApcCdh[k+1] * Cb[k+1]) ) * ( t_sim[k+1] - t_sim[k] ) )
@NLconstraint(m, [k in 1:length(t_sim)-1],
    ACE[k+1] == ACE[k] + ( -1.0*(kPhCdhEf4 * ACE[k+1] * Ce[k+1]) -1.0*(kPhCdhAf4 * ACE[k+1] * Ca[k+1]) -1.0*(kPhCdhBf4 * ACE[k+1] * Cb[k+1]) -1.0*(kPhApcA * ACE[k+1] * Ca[k+1]) -1.0*(kPhApcB * ACE[k+1] * Cb[k+1]) +1.0*(kDpApc * pACE[k+1] * B55[k+1]) -1.0*(kPhEmiA * ACE[k+1] * Ca[k+1]) -1.0*(kPhEmiB * ACE[k+1] * Cb[k+1]) -1.0*(kDeEmi1 * ACE[k+1]) -1.0*(kDeEmi2 * ACE[k+1]) +1.0*(kAsACE * ApcCdh[k+1] * Emi[k+1]) -1.0*(kDiACE * ACE[k+1]) ) * ( t_sim[k+1] - t_sim[k] ) )
@NLconstraint(m, [k in 1:length(t_sim)-1],
    Emi[k+1] == Emi[k] + ( +1.0*(kPhCdhEf4 * ACE[k+1] * Ce[k+1]) +1.0*(kPhCdhAf4 * ACE[k+1] * Ca[k+1]) +1.0*(kPhCdhBf4 * ACE[k+1] * Cb[k+1]) +1.0*(kPhCdhEf4 * pACE[k+1] * Ce[k+1]) +1.0*(kPhCdhAf4 * pACE[k+1] * Ca[k+1]) +1.0*(kPhCdhBf4 * pACE[k+1] * Cb[k+1]) -1.0*(kPhEmiA * Emi[k+1] * Ca[k+1]) -1.0*(kPhEmiB * Emi[k+1] * Cb[k+1]) +1.0*(kDpEmi * pEmi[k+1]) +1.0*(kSyEmi2 * E2fPx[k+1]) +1.0*(kSyEmi1) -1.0*(kDeEmi1 * Emi[k+1]) -1.0*(kAsACE * Emi[k+1] * pApcCdh[k+1]) +1.0*(kDiACE * pACE[k+1]) -1.0*(kAsACE * ApcCdh[k+1] * Emi[k+1]) +1.0*(kDiACE * ACE[k+1]) ) * ( t_sim[k+1] - t_sim[k] ) )
@NLconstraint(m, [k in 1:length(t_sim)-1],
    pACE[k+1] == pACE[k] + ( -1.0*(kPhCdhEf4 * pACE[k+1] * Ce[k+1]) -1.0*(kPhCdhAf4 * pACE[k+1] * Ca[k+1]) -1.0*(kPhCdhBf4 * pACE[k+1] * Cb[k+1]) -1.0*(kDeEmi1 * pACE[k+1]) -1.0*(kDeEmi2 * pACE[k+1]) -1.0*(kPhEmiA * pACE[k+1] * Ca[k+1]) -1.0*(kPhEmiB * pACE[k+1] * Cb[k+1]) +1.0*(kPhApcA * ACE[k+1] * Ca[k+1]) +1.0*(kPhApcB * ACE[k+1] * Cb[k+1]) -1.0*(kDpApc * pACE[k+1] * B55[k+1]) +1.0*(kAsACE * Emi[k+1] * pApcCdh[k+1]) -1.0*(kDiACE * pACE[k+1]) ) * ( t_sim[k+1] - t_sim[k] ) )
@NLconstraint(m, [k in 1:length(t_sim)-1],
    pApcCdh[k+1] == pApcCdh[k] + ( -1.0*(kPhCdhE * pApcCdh[k+1] * Ce[k+1]) -1.0*(kDpApc * pApcCdh[k+1] * B55[k+1]) +1.0*(kDeEmi1 * pACE[k+1]) +1.0*(kDeEmi2 * pACE[k+1]) +1.0*(kPhEmiA * pACE[k+1] * Ca[k+1]) +1.0*(kPhEmiB * pACE[k+1] * Cb[k+1]) -1.0*(kPhCdhA * pApcCdh[k+1] * Ca[k+1]) -1.0*(kPhCdhB * pApcCdh[k+1] * Cb[k+1]) +1.0*(kPhApcA * ApcCdh[k+1] * Ca[k+1]) +1.0*(kPhApcB * ApcCdh[k+1] * Cb[k+1]) +1.0*(kAsACdh * pApc[k+1] * Cdh[k+1]) -1.0*(kDiACdh * pApcCdh[k+1]) -1.0*(kAsACE * Emi[k+1] * pApcCdh[k+1]) +1.0*(kDiACE * pACE[k+1]) ) * ( t_sim[k+1] - t_sim[k] ) )
@NLconstraint(m, [k in 1:length(t_sim)-1],
    pEmi[k+1] == pEmi[k] + ( +1.0*(kPhEmiA * pACE[k+1] * Ca[k+1]) +1.0*(kPhEmiB * pACE[k+1] * Cb[k+1]) +1.0*(kPhEmiA * Emi[k+1] * Ca[k+1]) +1.0*(kPhEmiB * Emi[k+1] * Cb[k+1]) -1.0*(kDpEmi * pEmi[k+1]) +1.0*(kPhEmiA * ACE[k+1] * Ca[k+1]) +1.0*(kPhEmiB * ACE[k+1] * Cb[k+1]) -1.0*(kDeEmi1 * pEmi[k+1]) -1.0*(kDeEmi3 * pEmi[k+1]) ) * ( t_sim[k+1] - t_sim[k] ) )
@NLconstraint(m, [k in 1:length(t_sim)-1],
    Fox[k+1] == Fox[k] + ( +1.0*(kSyFox1) +1.0*(kSyFox2 * E2fPx[k+1]) -1.0*(kDeFox1 * Fox[k+1]) -1.0*(kDeFox2 * Fox[k+1] * ApcCdh[k+1]) -1.0*(kDeFox2 * Fox[k+1] * pApcCdh[k+1]) -1.0*(kPhFoxE * Fox[k+1] * Ce[k+1]) -1.0*(kPhFoxA * Fox[k+1] * Ca[k+1]) -1.0*(kPhFoxB * Fox[k+1] * Cb[k+1]) +1.0*(kDpFox * pFox[k+1]) ) * ( t_sim[k+1] - t_sim[k] ) )
@NLconstraint(m, [k in 1:length(t_sim)-1],
    pFox[k+1] == pFox[k] + ( -1.0*(kDeFox2 * pFox[k+1] * pApcCdh[k+1]) -1.0*(kDeFox1 * pFox[k+1]) -1.0*(kDeFox2 * pFox[k+1] * ApcCdh[k+1]) +1.0*(kPhFoxE * Fox[k+1] * Ce[k+1]) +1.0*(kPhFoxA * Fox[k+1] * Ca[k+1]) +1.0*(kPhFoxB * Fox[k+1] * Cb[k+1]) -1.0*(kDpFox * pFox[k+1]) ) * ( t_sim[k+1] - t_sim[k] ) )
@NLconstraint(m, [k in 1:length(t_sim)-1],
    pFoxPcb[k+1] == pFoxPcb[k] + ( +1.0*(kAsFPcb * pFox[k+1] * Procb[k+1]) -1.0*(kDiFPcb * pFoxPcb[k+1]) ) * ( t_sim[k+1] - t_sim[k] ) )
@NLconstraint(m, [k in 1:length(t_sim)-1],
    Procb[k+1] == Procb[k] + ( -1.0*(kAsFPcb * pFox[k+1] * Procb[k+1]) +1.0*(kDiFPcb * pFoxPcb[k+1]) ) * ( t_sim[k+1] - t_sim[k] ) )
@NLconstraint(m, [k in 1:length(t_sim)-1],
    Pcdc[k+1] == Pcdc[k] + ( -1.0*(kAsFPcdc * pFox[k+1] * Pcdc[k+1]) +1.0*(kDiFPcdc * pFoxPcdc[k+1]) ) * ( t_sim[k+1] - t_sim[k] ) )
@NLconstraint(m, [k in 1:length(t_sim)-1],
    pFoxPcdc[k+1] == pFoxPcdc[k] + ( +1.0*(kAsFPcdc * pFox[k+1] * Pcdc[k+1]) -1.0*(kDiFPcdc * pFoxPcdc[k+1]) ) * ( t_sim[k+1] - t_sim[k] ) )
@NLconstraint(m, [k in 1:length(t_sim)-1],
    Wee[k+1] == Wee[k] + ( +1.0*(kDpWee * pWee[k+1] * B55[k+1]) -1.0*(kPhWeeA * Wee[k+1] * Ca[k+1]) -1.0*(kPhWeeB * Wee[k+1] * Cb[k+1]) ) * ( t_sim[k+1] - t_sim[k] ) )
@NLconstraint(m, [k in 1:length(t_sim)-1],
    pCb[k+1] == pCb[k] + ( +1.0*(kWee2 * Cb[k+1] * Wee[k+1]) +1.0*(kWee1 * Cb[k+1] * pWee[k+1]) -1.0*(kCdc25_2 * pCb[k+1] * pCdc25[k+1]) -1.0*(kCdc25_1 * pCb[k+1] * Cdc25[k+1]) -1.0*(kDeCb1 * pCb[k+1]) -1.0*(kDeCb2 * pCb[k+1] * ApcCdh[k+1]) -1.0*(kDeCb2 * pCb[k+1] * pApcCdh[k+1]) -1.0*(kDeCb3 * pCb[k+1] * pApcCdc[k+1]) ) * ( t_sim[k+1] - t_sim[k] ) )
@NLconstraint(m, [k in 1:length(t_sim)-1],
    pWee[k+1] == pWee[k] + ( -1.0*(kDpWee * pWee[k+1] * B55[k+1]) +1.0*(kPhWeeA * Wee[k+1] * Ca[k+1]) +1.0*(kPhWeeB * Wee[k+1] * Cb[k+1]) ) * ( t_sim[k+1] - t_sim[k] ) )
@NLconstraint(m, [k in 1:length(t_sim)-1],
    pCdc25[k+1] == pCdc25[k] + ( +1.0*(kPhC25B * Cdc25[k+1] * Cb[k+1]) +1.0*(kPhC25A * Cdc25[k+1] * Ca[k+1]) -1.0*(kDpCdc25 * pCdc25[k+1] * B55[k+1]) ) * ( t_sim[k+1] - t_sim[k] ) )
@NLconstraint(m, [k in 1:length(t_sim)-1],
    Cdc25[k+1] == Cdc25[k] + ( -1.0*(kPhC25B * Cdc25[k+1] * Cb[k+1]) -1.0*(kPhC25A * Cdc25[k+1] * Ca[k+1]) +1.0*(kDpCdc25 * pCdc25[k+1] * B55[k+1]) ) * ( t_sim[k+1] - t_sim[k] ) )
@NLconstraint(m, [k in 1:length(t_sim)-1],
    Gw[k+1] == Gw[k] + ( -1.0*(kPhGw * Gw[k+1] * Cb[k+1]) +1.0*(kDpGw1 * pGw[k+1]) +1.0*(kDpGw2 * pGw[k+1] * B55[k+1]) ) * ( t_sim[k+1] - t_sim[k] ) )
@NLconstraint(m, [k in 1:length(t_sim)-1],
    pGw[k+1] == pGw[k] + ( +1.0*(kPhGw * Gw[k+1] * Cb[k+1]) -1.0*(kDpGw1 * pGw[k+1]) -1.0*(kDpGw2 * pGw[k+1] * B55[k+1]) ) * ( t_sim[k+1] - t_sim[k] ) )
@NLconstraint(m, [k in 1:length(t_sim)-1],
    Ensa[k+1] == Ensa[k] + ( -1.0*(kPhEnsa * Ensa[k+1] * pGw[k+1]) +1.0*(kDpEnsa * pEB55[k+1]) ) * ( t_sim[k+1] - t_sim[k] ) )
@NLconstraint(m, [k in 1:length(t_sim)-1],
    pEnsa[k+1] == pEnsa[k] + ( +1.0*(kPhEnsa * Ensa[k+1] * pGw[k+1]) -1.0*(kAspEB55 * pEnsa[k+1] * B55[k+1]) +1.0*(kDipEB55 * pEB55[k+1]) ) * ( t_sim[k+1] - t_sim[k] ) )
@NLconstraint(m, [k in 1:length(t_sim)-1],
    pEB55[k+1] == pEB55[k] + ( -1.0*(kDpEnsa * pEB55[k+1]) +1.0*(kAspEB55 * pEnsa[k+1] * B55[k+1]) -1.0*(kDipEB55 * pEB55[k+1]) ) * ( t_sim[k+1] - t_sim[k] ) )

# Define objective
println("Defining objective ...")
@variable(m, aux)
numerical_error_penalty = sum(weight*(Ce[k] - Ce[k-1])^pow for k in 2:length(t_sim))
    + sum(weight*(E2f[k] - E2f[k-1])^pow for k in 2:length(t_sim))
    + sum(weight*(pE2f[k] - pE2f[k-1])^pow for k in 2:length(t_sim))
    + sum(weight*(Ca[k] - Ca[k-1])^pow for k in 2:length(t_sim))
    + sum(weight*(Cb[k] - Cb[k-1])^pow for k in 2:length(t_sim))
    + sum(weight*(E2fRb[k] - E2fRb[k-1])^pow for k in 2:length(t_sim))
    + sum(weight*(pE2fRb[k] - pE2fRb[k-1])^pow for k in 2:length(t_sim))
    + sum(weight*(Rb[k] - Rb[k-1])^pow for k in 2:length(t_sim))
    + sum(weight*(pRb[k] - pRb[k-1])^pow for k in 2:length(t_sim))
    + sum(weight*(B55[k] - B55[k-1])^pow for k in 2:length(t_sim))
    + sum(weight*(E2fPx[k] - E2fPx[k-1])^pow for k in 2:length(t_sim))
    + sum(weight*(Px[k] - Px[k-1])^pow for k in 2:length(t_sim))
    + sum(weight*(Apc[k] - Apc[k-1])^pow for k in 2:length(t_sim))
    + sum(weight*(Cdc20[k] - Cdc20[k-1])^pow for k in 2:length(t_sim))
    + sum(weight*(pApcCdc[k] - pApcCdc[k-1])^pow for k in 2:length(t_sim))
    + sum(weight*(pApc[k] - pApc[k-1])^pow for k in 2:length(t_sim))
    + sum(weight*(ApcCdh[k] - ApcCdh[k-1])^pow for k in 2:length(t_sim))
    + sum(weight*(Cdh[k] - Cdh[k-1])^pow for k in 2:length(t_sim))
    + sum(weight*(pCdh[k] - pCdh[k-1])^pow for k in 2:length(t_sim))
    + sum(weight*(ACE[k] - ACE[k-1])^pow for k in 2:length(t_sim))
    + sum(weight*(Emi[k] - Emi[k-1])^pow for k in 2:length(t_sim))
    + sum(weight*(pACE[k] - pACE[k-1])^pow for k in 2:length(t_sim))
    + sum(weight*(pApcCdh[k] - pApcCdh[k-1])^pow for k in 2:length(t_sim))
    + sum(weight*(pEmi[k] - pEmi[k-1])^pow for k in 2:length(t_sim))
    + sum(weight*(Fox[k] - Fox[k-1])^pow for k in 2:length(t_sim))
    + sum(weight*(pFox[k] - pFox[k-1])^pow for k in 2:length(t_sim))
    + sum(weight*(pFoxPcb[k] - pFoxPcb[k-1])^pow for k in 2:length(t_sim))
    + sum(weight*(Procb[k] - Procb[k-1])^pow for k in 2:length(t_sim))
    + sum(weight*(Pcdc[k] - Pcdc[k-1])^pow for k in 2:length(t_sim))
    + sum(weight*(pFoxPcdc[k] - pFoxPcdc[k-1])^pow for k in 2:length(t_sim))
    + sum(weight*(Wee[k] - Wee[k-1])^pow for k in 2:length(t_sim))
    + sum(weight*(pCb[k] - pCb[k-1])^pow for k in 2:length(t_sim))
    + sum(weight*(pWee[k] - pWee[k-1])^pow for k in 2:length(t_sim))
    + sum(weight*(pCdc25[k] - pCdc25[k-1])^pow for k in 2:length(t_sim))
    + sum(weight*(Cdc25[k] - Cdc25[k-1])^pow for k in 2:length(t_sim))
    + sum(weight*(Gw[k] - Gw[k-1])^pow for k in 2:length(t_sim))
    + sum(weight*(pGw[k] - pGw[k-1])^pow for k in 2:length(t_sim))
    + sum(weight*(Ensa[k] - Ensa[k-1])^pow for k in 2:length(t_sim))
    + sum(weight*(pEnsa[k] - pEnsa[k-1])^pow for k in 2:length(t_sim))
    + sum(weight*(pEB55[k] - pEB55[k-1])^pow for k in 2:length(t_sim))
@constraint(m, aux == numerical_error_penalty)

@NLobjective(m, Min,sum((Ce[(k-1)*t_ratio+1]-df[k, :Ce])^2 for k in 1:length(t_exp))
    + sum((E2f[(k-1)*t_ratio+1]-df[k, :E2f])^2 for k in 1:length(t_exp))
    + sum((pE2f[(k-1)*t_ratio+1]-df[k, :pE2f])^2 for k in 1:length(t_exp))
    + sum((Ca[(k-1)*t_ratio+1]-df[k, :Ca])^2 for k in 1:length(t_exp))
    + sum((Cb[(k-1)*t_ratio+1]-df[k, :Cb])^2 for k in 1:length(t_exp))
    + sum((E2fRb[(k-1)*t_ratio+1]-df[k, :E2fRb])^2 for k in 1:length(t_exp))
    + sum((pE2fRb[(k-1)*t_ratio+1]-df[k, :pE2fRb])^2 for k in 1:length(t_exp))
    + sum((Rb[(k-1)*t_ratio+1]-df[k, :Rb])^2 for k in 1:length(t_exp))
    + sum((pRb[(k-1)*t_ratio+1]-df[k, :pRb])^2 for k in 1:length(t_exp))
    + sum((B55[(k-1)*t_ratio+1]-df[k, :B55])^2 for k in 1:length(t_exp))
    + sum((E2fPx[(k-1)*t_ratio+1]-df[k, :E2fPx])^2 for k in 1:length(t_exp))
    + sum((Px[(k-1)*t_ratio+1]-df[k, :Px])^2 for k in 1:length(t_exp))
    + sum((Apc[(k-1)*t_ratio+1]-df[k, :Apc])^2 for k in 1:length(t_exp))
    + sum((Cdc20[(k-1)*t_ratio+1]-df[k, :Cdc20])^2 for k in 1:length(t_exp))
    + sum((pApcCdc[(k-1)*t_ratio+1]-df[k, :pApcCdc])^2 for k in 1:length(t_exp))
    + sum((pApc[(k-1)*t_ratio+1]-df[k, :pApc])^2 for k in 1:length(t_exp))
    + sum((ApcCdh[(k-1)*t_ratio+1]-df[k, :ApcCdh])^2 for k in 1:length(t_exp))
    + sum((Cdh[(k-1)*t_ratio+1]-df[k, :Cdh])^2 for k in 1:length(t_exp))
    + sum((pCdh[(k-1)*t_ratio+1]-df[k, :pCdh])^2 for k in 1:length(t_exp))
    + sum((ACE[(k-1)*t_ratio+1]-df[k, :ACE])^2 for k in 1:length(t_exp))
    + sum((Emi[(k-1)*t_ratio+1]-df[k, :Emi])^2 for k in 1:length(t_exp))
    + sum((pACE[(k-1)*t_ratio+1]-df[k, :pACE])^2 for k in 1:length(t_exp))
    + sum((pApcCdh[(k-1)*t_ratio+1]-df[k, :pApcCdh])^2 for k in 1:length(t_exp))
    + sum((pEmi[(k-1)*t_ratio+1]-df[k, :pEmi])^2 for k in 1:length(t_exp))
    + sum((Fox[(k-1)*t_ratio+1]-df[k, :Fox])^2 for k in 1:length(t_exp))
    + sum((pFox[(k-1)*t_ratio+1]-df[k, :pFox])^2 for k in 1:length(t_exp))
    + sum((pFoxPcb[(k-1)*t_ratio+1]-df[k, :pFoxPcb])^2 for k in 1:length(t_exp))
    + sum((Procb[(k-1)*t_ratio+1]-df[k, :Procb])^2 for k in 1:length(t_exp))
    + sum((Pcdc[(k-1)*t_ratio+1]-df[k, :Pcdc])^2 for k in 1:length(t_exp))
    + sum((pFoxPcdc[(k-1)*t_ratio+1]-df[k, :pFoxPcdc])^2 for k in 1:length(t_exp))
    + sum((Wee[(k-1)*t_ratio+1]-df[k, :Wee])^2 for k in 1:length(t_exp))
    + sum((pCb[(k-1)*t_ratio+1]-df[k, :pCb])^2 for k in 1:length(t_exp))
    + sum((pWee[(k-1)*t_ratio+1]-df[k, :pWee])^2 for k in 1:length(t_exp))
    + sum((pCdc25[(k-1)*t_ratio+1]-df[k, :pCdc25])^2 for k in 1:length(t_exp))
    + sum((Cdc25[(k-1)*t_ratio+1]-df[k, :Cdc25])^2 for k in 1:length(t_exp))
    + sum((Gw[(k-1)*t_ratio+1]-df[k, :Gw])^2 for k in 1:length(t_exp))
    + sum((pGw[(k-1)*t_ratio+1]-df[k, :pGw])^2 for k in 1:length(t_exp))
    + sum((Ensa[(k-1)*t_ratio+1]-df[k, :Ensa])^2 for k in 1:length(t_exp))
    + sum((pEnsa[(k-1)*t_ratio+1]-df[k, :pEnsa])^2 for k in 1:length(t_exp))
    + sum((pEB55[(k-1)*t_ratio+1]-df[k, :pEB55])^2 for k in 1:length(t_exp))
    + aux
    )

println("Optimizing...")
optimize!(m)


println("# Obtain the solution")
println("Retreiving solution...")
species_to_plot = []
params = [kPhRbD, kPhRbE, kDpRb, kPhRb, kSyCe1, kSyCe2, kDeCe, kDiE2fRb, kAsE2fRb, kSyE2f1, kSyE2f2, kDeE2f1, kDeE2f2, kPhE2f, kDpE2f1, kDpE2f2, kSyEmi1, kSyEmi2, kDeEmi1, kDeEmi2, kDeEmi3, kPhEmiA, kPhEmiB, kDpEmi, kSyCa1, kSyCa2, kDeCa1, kDeCa2, kDeCa3, kDpCdh, kPhCdhA, kPhCdhE, kPhCdhB, kDiACE, kAsACE, kDiACdh, kAsACdh, kPhCeA, kPhCeE, kPhFoxE, kPhFoxA, kPhFoxB, kDpFox, kSyCb1, kSyCb2, kDeCb1, kDeCb2, kDeCb3, kPhEnsa, kDpEnsa, kPhGw, kDpGw1, kDpGw2, kWee1, kWee2, kPhWeeA, kPhWeeB, kDpWee, kCdc25_1, kCdc25_2, kPhC25A, kDpCdc25, kDipEB55, kAspEB55, kSyCdc_1, kSyCdc_2, kDeCdc_1, kDeCdc_2, kDipACdc, kAspACdc, kPhApcA, kPhApcB, kDpApc, kSyFox1, kSyFox2, kDeFox1, kDeFox2, kAsEPx, kDiEPx, kAsFPcb, kDiFPcb, kAsFPcdc, kDiFPcdc, kPhC25B, kPhCdhEf4, kPhCdhAf4, kPhCdhBf4]
paramvalues = Dict()
for param in params
    paramvalues[param] = JuMP.value.(param)
end

variables = [ :Ce
 :E2f
 :pE2f
 :Ca
 :Cb
 :E2fRb
 :pE2fRb
 :Rb
 :pRb
 :B55
 :E2fPx
 :Px
 :Apc
 :Cdc20
 :pApcCdc
 :pApc
 :ApcCdh
 :Cdh
 :pCdh
 :ACE
 :Emi
 :pACE
 :pApcCdh
 :pEmi
 :Fox
 :pFox
 :pFoxPcb
 :Procb
 :Pcdc
 :pFoxPcdc
 :Wee
 :pCb
 :pWee
 :pCdc25
 :Cdc25
 :Gw
 :pGw
 :Ensa
 :pEnsa
 :pEB55
]
variablevalues = Dict()
for v in variables
    variablevalues[string(v)] = Vector(JuMP.value.(eval(v)))
end

data_matrix = zeros(length(t_sim), length(species_to_plot))
i = 1
for s in species_to_plot
    data_matrix[:, i] = variablevalues[string(s)]
i = i+1
end


using Plots
using DataFrames
using StatPlots
p = plot(t_sim, data_matrix, xlabel="Time (hr)", ylabel="Abundance", label=species_to_plot
, legend=:topleft)

@df df plot!(p, t_exp, seriestype=:scatter, markershape = :x,
    markersize = 2,
    markerstrokewidth = 0.1,
    markerstrokealpha = 1)