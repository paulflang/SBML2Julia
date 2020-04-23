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
@variable(m, 0.999999999999159/fc <= tRb <= 0.999999999999159*fc)
@variable(m, 0.999999999999977/fc <= tPcb <= 0.999999999999977*fc)
@variable(m, 5.00000000000149/fc <= tApc <= 5.00000000000149*fc)
@variable(m, 0.14868387326335/fc <= t_pApc <= 0.14868387326335*fc)
@variable(m, 1.00000000000128/fc <= tCdh <= 1.00000000000128*fc)
@variable(m, 1.00000000000001/fc <= tEnsa <= 1.00000000000001*fc)
@variable(m, 0.249999999999996/fc <= tB55 <= 0.249999999999996*fc)
@variable(m, 1.00000000000001/fc <= tGw <= 1.00000000000001*fc)
@variable(m, 1.0/fc <= tCdc25 <= 1.0*fc)
@variable(m, 0.999999999999995/fc <= tWee <= 0.999999999999995*fc)
@variable(m, 1.0/fc <= tPx <= 1.0*fc)
@variable(m, 1.0/fc <= tPcdc <= 1.0*fc)
@variable(m, 0.00939542544113453/fc <= tCb <= 0.00939542544113453*fc)
@variable(m, 0.111893832351909/fc <= tEmi <= 0.111893832351909*fc)
@variable(m, 0.2998240933904/fc <= tCdc20 <= 0.2998240933904*fc)
@variable(m, 0.001/fc <= kSyE2f1 <= 0.001*fc)
@variable(m, 50.0/fc <= kAsFPcb <= 50.0*fc)
@variable(m, 0.01/fc <= kSyE2f2 <= 0.01*fc)
@variable(m, 0.01/fc <= kDeE2f1 <= 0.01*fc)
@variable(m, 5.0/fc <= kDiEPx <= 5.0*fc)
@variable(m, 50.0/fc <= kAsE2fRb <= 50.0*fc)
@variable(m, 0.05/fc <= kDiE2fRb <= 0.05*fc)
@variable(m, 0.8/fc <= kPhRbD <= 0.8*fc)
@variable(m, 1.0/fc <= Cd <= 1.0*fc)
@variable(m, 1.0/fc <= kPhRbE <= 1.0*fc)
@variable(m, 3.0/fc <= kPhRb <= 3.0*fc)
@variable(m, 0.2/fc <= kPhE2f <= 0.2*fc)
@variable(m, 0.0/fc <= kDpE2f1 <= 0.0*fc)
@variable(m, 10.0/fc <= kDpE2f2 <= 10.0*fc)
@variable(m, 0.2/fc <= kDeE2f2 <= 0.2*fc)
@variable(m, 0.002/fc <= kSyCe1 <= 0.002*fc)
@variable(m, 0.0075/fc <= kSyCa2 <= 0.0075*fc)
@variable(m, 0.02/fc <= kDeCe <= 0.02*fc)
@variable(m, 0.01/fc <= kPhCeE <= 0.01*fc)
@variable(m, 0.2/fc <= kDpRb <= 0.2*fc)
@variable(m, 5.0/fc <= kDpApc <= 5.0*fc)
@variable(m, 0.01/fc <= kPhApcA <= 0.01*fc)
@variable(m, 0.01/fc <= kPhApcB <= 0.01*fc)
@variable(m, 50.0/fc <= kAsACdh <= 50.0*fc)
@variable(m, 0.05/fc <= kDiACdh <= 0.05*fc)
@variable(m, 2.5/fc <= kPhCdhA <= 2.5*fc)
@variable(m, 0.1/fc <= kPhCdhE <= 0.1*fc)
@variable(m, 25.0/fc <= kPhCdhB <= 25.0*fc)
@variable(m, 0.05/fc <= kDeEmi1 <= 0.05*fc)
@variable(m, 0.05/fc <= kDeEmi2 <= 0.05*fc)
@variable(m, 0.005/fc <= kPhEmiA <= 0.005*fc)
@variable(m, 0.025/fc <= kPhEmiB <= 0.025*fc)
@variable(m, 50.0/fc <= kAsACE <= 50.0*fc)
@variable(m, 0.05/fc <= kDiACE <= 0.05*fc)
@variable(m, 0.0005/fc <= kDpEmi <= 0.0005*fc)
@variable(m, 0.5/fc <= kDeEmi3 <= 0.5*fc)
@variable(m, 0.001/fc <= kSyFox1 <= 0.001*fc)
@variable(m, 0.05/fc <= kSyEmi2 <= 0.05*fc)
@variable(m, 0.5/fc <= kDpCdh <= 0.5*fc)
@variable(m, 0.001/fc <= kSyCa1 <= 0.001*fc)
@variable(m, 0.001/fc <= kDeCa1 <= 0.001*fc)
@variable(m, 1.0/fc <= kPhCeA <= 1.0*fc)
@variable(m, 0.4/fc <= kDeCa2 <= 0.4*fc)
@variable(m, 0.1/fc <= kDeCa3 <= 0.1*fc)
@variable(m, 0.0075/fc <= kSyFox2 <= 0.0075*fc)
@variable(m, 0.001/fc <= kDeFox1 <= 0.001*fc)
@variable(m, 0.5/fc <= kDeFox2 <= 0.5*fc)
@variable(m, 0.1/fc <= kPhFoxE <= 0.1*fc)
@variable(m, 0.1/fc <= kPhFoxA <= 0.1*fc)
@variable(m, 0.5/fc <= kPhFoxB <= 0.5*fc)
@variable(m, 0.005/fc <= kDpFox <= 0.005*fc)
@variable(m, 0.001/fc <= kSyCb1 <= 0.001*fc)
@variable(m, 0.0075/fc <= kSyCb2 <= 0.0075*fc)
@variable(m, 50.0/fc <= kAsEPx <= 50.0*fc)
@variable(m, 10.0/fc <= kDiFPcb <= 10.0*fc)
@variable(m, 0.001/fc <= kDeCb1 <= 0.001*fc)
@variable(m, 0.125/fc <= kDeCb2 <= 0.125*fc)
@variable(m, 0.125/fc <= kDeCb3 <= 0.125*fc)
@variable(m, 1.0/fc <= kAspACdc <= 1.0*fc)
@variable(m, 0.125/fc <= kDipACdc <= 0.125*fc)
@variable(m, 0.004/fc <= kDeCdc_1 <= 0.004*fc)
@variable(m, 0.04/fc <= kDeCdc_2 <= 0.04*fc)
@variable(m, 0.001/fc <= kSyCdc_1 <= 0.001*fc)
@variable(m, 0.01/fc <= kSyCdc_2 <= 0.01*fc)
@variable(m, 50.0/fc <= kAsFPcdc <= 50.0*fc)
@variable(m, 25.0/fc <= kDiFPcdc <= 25.0*fc)
@variable(m, 1.0/fc <= kWee2 <= 1.0*fc)
@variable(m, 0.01/fc <= kWee1 <= 0.01*fc)
@variable(m, 1.0/fc <= kCdc25_2 <= 1.0*fc)
@variable(m, 0.1/fc <= kCdc25_1 <= 0.1*fc)
@variable(m, 1.0/fc <= kPhC25B <= 1.0*fc)
@variable(m, 0.0/fc <= kPhC25A <= 0.0*fc)
@variable(m, 10.0/fc <= kDpCdc25 <= 10.0*fc)
@variable(m, 10.0/fc <= kDpWee <= 10.0*fc)
@variable(m, 0.1/fc <= kPhWeeA <= 0.1*fc)
@variable(m, 1.0/fc <= kPhWeeB <= 1.0*fc)
@variable(m, 1.0/fc <= kPhGw <= 1.0*fc)
@variable(m, 0.25/fc <= kDpGw1 <= 0.25*fc)
@variable(m, 10.0/fc <= kDpGw2 <= 10.0*fc)
@variable(m, 0.1/fc <= kPhEnsa <= 0.1*fc)
@variable(m, 57.0/fc <= kAspEB55 <= 57.0*fc)
@variable(m, 0.0068/fc <= kDipEB55 <= 0.0068*fc)
@variable(m, 0.05/fc <= kDpEnsa <= 0.05*fc)
@variable(m, 0.005/fc <= kSyEmi1 <= 0.005*fc)
@variable(m, 0.06/fc <= kSyCe2 <= 0.06*fc)

# Model states
println("Defining states ...")
@variable(m, 0 <= Rb[k in 1:length(t_sim)] <= 1)
@variable(m, 0 <= pRb[k in 1:length(t_sim)] <= 1)
@variable(m, 0 <= pE2fRb[k in 1:length(t_sim)] <= 1)
@variable(m, 0 <= E2fRb[k in 1:length(t_sim)] <= 1)
@variable(m, 0 <= Pcb[k in 1:length(t_sim)] <= 1)
@variable(m, 0 <= pFoxPcb[k in 1:length(t_sim)] <= 1)
@variable(m, 0 <= ACE[k in 1:length(t_sim)] <= 1)
@variable(m, 0 <= Apc[k in 1:length(t_sim)] <= 1)
@variable(m, 0 <= ApcCdh[k in 1:length(t_sim)] <= 1)
@variable(m, 0 <= pACE[k in 1:length(t_sim)] <= 1)
@variable(m, 0 <= pApc[k in 1:length(t_sim)] <= 1)
@variable(m, 0 <= pApcCdc[k in 1:length(t_sim)] <= 1)
@variable(m, 0 <= pApcCdh[k in 1:length(t_sim)] <= 1)
@variable(m, 0 <= Cdh[k in 1:length(t_sim)] <= 1)
@variable(m, 0 <= pCdh[k in 1:length(t_sim)] <= 1)
@variable(m, 0 <= Ensa[k in 1:length(t_sim)] <= 1)
@variable(m, 0 <= pEB55[k in 1:length(t_sim)] <= 1)
@variable(m, 0 <= pEnsa[k in 1:length(t_sim)] <= 1)
@variable(m, 0 <= B55[k in 1:length(t_sim)] <= 1)
@variable(m, 0 <= Gw[k in 1:length(t_sim)] <= 1)
@variable(m, 0 <= pGw[k in 1:length(t_sim)] <= 1)
@variable(m, 0 <= Cdc25[k in 1:length(t_sim)] <= 1)
@variable(m, 0 <= pCdc25[k in 1:length(t_sim)] <= 1)
@variable(m, 0 <= pWee[k in 1:length(t_sim)] <= 1)
@variable(m, 0 <= Wee[k in 1:length(t_sim)] <= 1)
@variable(m, 0 <= E2fPx[k in 1:length(t_sim)] <= 1)
@variable(m, 0 <= Px[k in 1:length(t_sim)] <= 1)
@variable(m, 0 <= Pcdc[k in 1:length(t_sim)] <= 1)
@variable(m, 0 <= pFoxPcdc[k in 1:length(t_sim)] <= 1)
@variable(m, 0 <= Cb[k in 1:length(t_sim)] <= 1)
@variable(m, 0 <= pCb[k in 1:length(t_sim)] <= 1)
@variable(m, 0 <= Emi[k in 1:length(t_sim)] <= 1)
@variable(m, 0 <= pEmi[k in 1:length(t_sim)] <= 1)
@variable(m, 0 <= Cdc20[k in 1:length(t_sim)] <= 1)
@variable(m, 0 <= E2f[k in 1:length(t_sim)] <= 1)
@variable(m, 0 <= pFox[k in 1:length(t_sim)] <= 1)
@variable(m, 0 <= Ce[k in 1:length(t_sim)] <= 1)
@variable(m, 0 <= Ca[k in 1:length(t_sim)] <= 1)
@variable(m, 0 <= pE2f[k in 1:length(t_sim)] <= 1)
@variable(m, 0 <= Fox[k in 1:length(t_sim)] <= 1)

# Model ODEs
println("Defining ODEs ...")
@NLconstraint(m, [k in 1:length(t_sim)-1],
    Rb[k+1] == Rb[k] + ( -1.0*(E2f * Rb - kDiE2fRb[k+1] * E2fRb[k+1]) -1.0*(pE2f * Rb - kDiE2fRb[k+1] * pE2fRb[k+1]) +1.0*(E2fRb) +1.0*(pE2fRb) +1.0*(pE2fRb) -1.0*(Rb * Cd[k+1]) -1.0*(Rb * Ce[k+1]) -1.0*(Rb * Ca[k+1]) -1.0*(Rb * Cb[k+1]) +1.0*(pRb) ) * ( t_sim[k+1] - t_sim[k] ) )
@NLconstraint(m, [k in 1:length(t_sim)-1],
    pRb[k+1] == pRb[k] + ( +1.0*(E2fRb * Cd[k+1]) +1.0*(E2fRb * Ce[k+1]) +1.0*(E2fRb * Ca[k+1]) +1.0*(E2fRb * Cb[k+1]) +1.0*(pE2fRb * Cd[k+1]) +1.0*(pE2fRb * Ca[k+1]) +1.0*(pE2fRb * Cb[k+1]) +1.0*(pE2fRb * Ce[k+1]) +1.0*(Rb * Cd[k+1]) +1.0*(Rb * Ce[k+1]) +1.0*(Rb * Ca[k+1]) +1.0*(Rb * Cb[k+1]) -1.0*(pRb) ) * ( t_sim[k+1] - t_sim[k] ) )
@NLconstraint(m, [k in 1:length(t_sim)-1],
    pE2fRb[k+1] == pE2fRb[k] + ( -1.0*(pE2fRb * B55[k+1]) +1.0*(pE2f * Rb - kDiE2fRb[k+1] * pE2fRb[k+1]) -1.0*(pE2fRb * Cd[k+1]) -1.0*(pE2fRb * Ca[k+1]) -1.0*(pE2fRb * Cb[k+1]) -1.0*(pE2fRb * Ce[k+1]) +1.0*(E2fRb * Ca[k+1]) +1.0*(E2fRb * Cb[k+1]) -1.0*(pE2fRb) -1.0*(pE2fRb) -1.0*(pE2fRb) ) * ( t_sim[k+1] - t_sim[k] ) )
@NLconstraint(m, [k in 1:length(t_sim)-1],
    E2fRb[k+1] == E2fRb[k] + ( +1.0*(E2f * Rb - kDiE2fRb[k+1] * E2fRb[k+1]) -1.0*(E2fRb * Cd[k+1]) -1.0*(E2fRb * Ce[k+1]) -1.0*(E2fRb * Ca[k+1]) -1.0*(E2fRb * Cb[k+1]) +1.0*(pE2fRb * B55[k+1]) -1.0*(E2fRb * Ca[k+1]) -1.0*(E2fRb * Cb[k+1]) +1.0*(pE2fRb) -1.0*(E2fRb) ) * ( t_sim[k+1] - t_sim[k] ) )
@NLconstraint(m, [k in 1:length(t_sim)-1],
    Pcb[k+1] == Pcb[k] + ( -1.0*(pFox * Pcb[k+1]) +1.0*(pFoxPcb) ) * ( t_sim[k+1] - t_sim[k] ) )
@NLconstraint(m, [k in 1:length(t_sim)-1],
    pFoxPcb[k+1] == pFoxPcb[k] + ( +1.0*(pFox * Pcb[k+1]) -1.0*(pFoxPcb) +1.0*(pFoxPcb) -1.0*(pFoxPcb) ) * ( t_sim[k+1] - t_sim[k] ) )
@NLconstraint(m, [k in 1:length(t_sim)-1],
    ACE[k+1] == ACE[k] + ( -1.0*(ACE * Ce[k+1] * 0.8[k+1]) -1.0*(ACE * Ca[k+1] * 0.8[k+1]) -1.0*(ACE * Cb[k+1] * 0.8[k+1]) -1.0*(ACE * Ca[k+1]) -1.0*(ACE * Cb[k+1]) +1.0*(pACE * B55[k+1]) -1.0*(ACE * Ca[k+1]) -1.0*(ACE * Cb[k+1]) +1.0*(ApcCdh * Emi - kDiACE[k+1] * ACE[k+1]) -1.0*(ACE) -1.0*(ACE) ) * ( t_sim[k+1] - t_sim[k] ) )
@NLconstraint(m, [k in 1:length(t_sim)-1],
    Apc[k+1] == Apc[k] + ( +1.0*(pApcCdc * B55[k+1]) -1.0*(Apc * Ca[k+1]) -1.0*(Apc * Cb[k+1]) +1.0*(pApc * B55[k+1]) -1.0*(Cdh * Apc - kDiACdh[k+1] * ApcCdh[k+1]) +1.0*(ApcCdh * Ca[k+1]) +1.0*(ACE * Ce[k+1] * 0.8[k+1]) +1.0*(ACE * Ca[k+1] * 0.8[k+1]) +1.0*(ACE * Cb[k+1] * 0.8[k+1]) +1.0*(ApcCdh * Ce[k+1]) +1.0*(ApcCdh * Cb[k+1]) ) * ( t_sim[k+1] - t_sim[k] ) )
@NLconstraint(m, [k in 1:length(t_sim)-1],
    ApcCdh[k+1] == ApcCdh[k] + ( +1.0*(Cdh * Apc - kDiACdh[k+1] * ApcCdh[k+1]) -1.0*(ApcCdh * Ca[k+1]) +1.0*(pApcCdh * B55[k+1]) +1.0*(ACE * Ca[k+1]) +1.0*(ACE * Cb[k+1]) -1.0*(ApcCdh * Emi - kDiACE[k+1] * ACE[k+1]) -1.0*(Ca * ApcCdh[k+1]) +1.0*(Ca * ApcCdh[k+1]) -1.0*(Fox * ApcCdh[k+1]) +1.0*(Fox * ApcCdh[k+1]) -1.0*(pFox * ApcCdh[k+1]) +1.0*(pFox * ApcCdh[k+1]) -1.0*(Cb * ApcCdh[k+1]) +1.0*(Cb * ApcCdh[k+1]) -1.0*(pApcCdc * ApcCdh[k+1]) +1.0*(pApcCdc * ApcCdh[k+1]) -1.0*(pCb * ApcCdh[k+1]) +1.0*(pCb * ApcCdh[k+1]) -1.0*(Cdc20 * ApcCdh[k+1]) +1.0*(Cdc20 * ApcCdh[k+1]) -1.0*(ApcCdh * Ce[k+1]) -1.0*(ApcCdh * Cb[k+1]) +1.0*(ACE) +1.0*(ACE) -1.0*(ApcCdh * Ca[k+1]) -1.0*(ApcCdh * Cb[k+1]) ) * ( t_sim[k+1] - t_sim[k] ) )
@NLconstraint(m, [k in 1:length(t_sim)-1],
    pACE[k+1] == pACE[k] + ( -1.0*(pACE * Ce[k+1] * 0.8[k+1]) -1.0*(pACE * Ca[k+1] * 0.8[k+1]) -1.0*(pACE * Cb[k+1] * 0.8[k+1]) -1.0*(pACE) -1.0*(pACE) -1.0*(pACE * Ca[k+1]) -1.0*(pACE * Cb[k+1]) +1.0*(Emi * pApcCdh - kDiACE[k+1] * pACE[k+1]) +1.0*(ACE * Ca[k+1]) +1.0*(ACE * Cb[k+1]) -1.0*(pACE * B55[k+1]) ) * ( t_sim[k+1] - t_sim[k] ) )
@NLconstraint(m, [k in 1:length(t_sim)-1],
    pApc[k+1] == pApc[k] + ( +1.0*(Apc * Ca[k+1]) +1.0*(Apc * Cb[k+1]) -1.0*(pApc * B55[k+1]) +1.0*(pACE * Ce[k+1] * 0.8[k+1]) +1.0*(pACE * Ca[k+1] * 0.8[k+1]) +1.0*(pACE * Cb[k+1] * 0.8[k+1]) +1.0*(pApcCdh * Ce[k+1]) -1.0*(pApc * Cdh - kDiACdh[k+1] * pApcCdh[k+1]) -1.0*(pApc * Cdc20 - kDipACdc[k+1] * pApcCdc[k+1]) +1.0*(pApcCdc) +1.0*(pApcCdc * ApcCdh[k+1]) +1.0*(pApcCdc * pApcCdh[k+1]) +1.0*(pApcCdh * Ca[k+1]) +1.0*(pApcCdh * Cb[k+1]) ) * ( t_sim[k+1] - t_sim[k] ) )
@NLconstraint(m, [k in 1:length(t_sim)-1],
    pApcCdc[k+1] == pApcCdc[k] + ( -1.0*(pApcCdc * B55[k+1]) -1.0*(Ca * pApcCdc[k+1]) +1.0*(Ca * pApcCdc[k+1]) -1.0*(Cb * pApcCdc[k+1]) +1.0*(Cb * pApcCdc[k+1]) +1.0*(pApc * Cdc20 - kDipACdc[k+1] * pApcCdc[k+1]) -1.0*(pApcCdc) -1.0*(pApcCdc * ApcCdh[k+1]) -1.0*(pApcCdc * pApcCdh[k+1]) -1.0*(pCb * pApcCdc[k+1]) +1.0*(pCb * pApcCdc[k+1]) ) * ( t_sim[k+1] - t_sim[k] ) )
@NLconstraint(m, [k in 1:length(t_sim)-1],
    pApcCdh[k+1] == pApcCdh[k] + ( -1.0*(pApcCdh * Ce[k+1]) -1.0*(pApcCdh * B55[k+1]) +1.0*(pApc * Cdh - kDiACdh[k+1] * pApcCdh[k+1]) +1.0*(pACE) +1.0*(pACE) +1.0*(pACE * Ca[k+1]) +1.0*(pACE * Cb[k+1]) -1.0*(Emi * pApcCdh - kDiACE[k+1] * pACE[k+1]) -1.0*(Ca * pApcCdh[k+1]) +1.0*(Ca * pApcCdh[k+1]) -1.0*(pFox * pApcCdh[k+1]) +1.0*(pFox * pApcCdh[k+1]) -1.0*(Fox * pApcCdh[k+1]) +1.0*(Fox * pApcCdh[k+1]) -1.0*(Cb * pApcCdh[k+1]) +1.0*(Cb * pApcCdh[k+1]) -1.0*(pApcCdc * pApcCdh[k+1]) +1.0*(pApcCdc * pApcCdh[k+1]) -1.0*(pCb * pApcCdh[k+1]) +1.0*(pCb * pApcCdh[k+1]) -1.0*(Cdc20 * pApcCdh[k+1]) +1.0*(Cdc20 * pApcCdh[k+1]) -1.0*(pApcCdh * Ca[k+1]) -1.0*(pApcCdh * Cb[k+1]) +1.0*(ApcCdh * Ca[k+1]) +1.0*(ApcCdh * Cb[k+1]) ) * ( t_sim[k+1] - t_sim[k] ) )
@NLconstraint(m, [k in 1:length(t_sim)-1],
    Cdh[k+1] == Cdh[k] + ( -1.0*(Cdh * Apc - kDiACdh[k+1] * ApcCdh[k+1]) -1.0*(pApc * Cdh - kDiACdh[k+1] * pApcCdh[k+1]) -1.0*(Cdh * Ce[k+1]) -1.0*(Cdh * Ca[k+1]) -1.0*(Cdh * Cb[k+1]) +1.0*(pCdh) ) * ( t_sim[k+1] - t_sim[k] ) )
@NLconstraint(m, [k in 1:length(t_sim)-1],
    pCdh[k+1] == pCdh[k] + ( +1.0*(ApcCdh * Ca[k+1]) +1.0*(ACE * Ce[k+1] * 0.8[k+1]) +1.0*(ACE * Ca[k+1] * 0.8[k+1]) +1.0*(ACE * Cb[k+1] * 0.8[k+1]) +1.0*(pACE * Ce[k+1] * 0.8[k+1]) +1.0*(pACE * Ca[k+1] * 0.8[k+1]) +1.0*(pACE * Cb[k+1] * 0.8[k+1]) +1.0*(pApcCdh * Ce[k+1]) +1.0*(Cdh * Ce[k+1]) +1.0*(Cdh * Ca[k+1]) +1.0*(Cdh * Cb[k+1]) -1.0*(pCdh) +1.0*(pApcCdh * Ca[k+1]) +1.0*(pApcCdh * Cb[k+1]) +1.0*(ApcCdh * Ce[k+1]) +1.0*(ApcCdh * Cb[k+1]) ) * ( t_sim[k+1] - t_sim[k] ) )
@NLconstraint(m, [k in 1:length(t_sim)-1],
    Ensa[k+1] == Ensa[k] + ( -1.0*(Ensa * pGw[k+1]) +1.0*(pEB55) ) * ( t_sim[k+1] - t_sim[k] ) )
@NLconstraint(m, [k in 1:length(t_sim)-1],
    pEB55[k+1] == pEB55[k] + ( +1.0*(pEnsa * B55 - kDipEB55[k+1] * pEB55[k+1]) -1.0*(pEB55) ) * ( t_sim[k+1] - t_sim[k] ) )
@NLconstraint(m, [k in 1:length(t_sim)-1],
    pEnsa[k+1] == pEnsa[k] + ( +1.0*(Ensa * pGw[k+1]) -1.0*(pEnsa * B55 - kDipEB55[k+1] * pEB55[k+1]) ) * ( t_sim[k+1] - t_sim[k] ) )
@NLconstraint(m, [k in 1:length(t_sim)-1],
    B55[k+1] == B55[k] + ( -1.0*(pE2fRb * B55[k+1]) +1.0*(pE2fRb * B55[k+1]) -1.0*(pApc * B55[k+1]) +1.0*(pApc * B55[k+1]) -1.0*(pCdc25 * B55[k+1]) +1.0*(pCdc25 * B55[k+1]) -1.0*(pWee * B55[k+1]) +1.0*(pWee * B55[k+1]) -1.0*(pGw * B55[k+1]) +1.0*(pGw * B55[k+1]) -1.0*(pEnsa * B55 - kDipEB55[k+1] * pEB55[k+1]) +1.0*(pEB55) -1.0*(pE2f * B55[k+1]) +1.0*(pE2f * B55[k+1]) ) * ( t_sim[k+1] - t_sim[k] ) )
@NLconstraint(m, [k in 1:length(t_sim)-1],
    Gw[k+1] == Gw[k] + ( -1.0*(Gw * Cb[k+1]) +1.0*(pGw) +1.0*(pGw * B55[k+1]) ) * ( t_sim[k+1] - t_sim[k] ) )
@NLconstraint(m, [k in 1:length(t_sim)-1],
    pGw[k+1] == pGw[k] + ( +1.0*(Gw * Cb[k+1]) -1.0*(pGw) -1.0*(pGw * B55[k+1]) -1.0*(Ensa * pGw[k+1]) +1.0*(Ensa * pGw[k+1]) ) * ( t_sim[k+1] - t_sim[k] ) )
@NLconstraint(m, [k in 1:length(t_sim)-1],
    Cdc25[k+1] == Cdc25[k] + ( -1.0*(pCb * Cdc25[k+1]) +1.0*(pCb * Cdc25[k+1]) -1.0*(Cdc25 * Cb[k+1]) -1.0*(Cdc25 * Ca[k+1]) +1.0*(pCdc25 * B55[k+1]) ) * ( t_sim[k+1] - t_sim[k] ) )
@NLconstraint(m, [k in 1:length(t_sim)-1],
    pCdc25[k+1] == pCdc25[k] + ( -1.0*(pCb * pCdc25[k+1]) +1.0*(pCb * pCdc25[k+1]) +1.0*(Cdc25 * Cb[k+1]) +1.0*(Cdc25 * Ca[k+1]) -1.0*(pCdc25 * B55[k+1]) ) * ( t_sim[k+1] - t_sim[k] ) )
@NLconstraint(m, [k in 1:length(t_sim)-1],
    pWee[k+1] == pWee[k] + ( -1.0*(Cb * pWee[k+1]) +1.0*(Cb * pWee[k+1]) -1.0*(pWee * B55[k+1]) +1.0*(Wee * Ca[k+1]) +1.0*(Wee * Cb[k+1]) ) * ( t_sim[k+1] - t_sim[k] ) )
@NLconstraint(m, [k in 1:length(t_sim)-1],
    Wee[k+1] == Wee[k] + ( -1.0*(Cb * Wee[k+1]) +1.0*(Cb * Wee[k+1]) +1.0*(pWee * B55[k+1]) -1.0*(Wee * Ca[k+1]) -1.0*(Wee * Cb[k+1]) ) * ( t_sim[k+1] - t_sim[k] ) )
@NLconstraint(m, [k in 1:length(t_sim)-1],
    E2fPx[k+1] == E2fPx[k] + ( -1.0*(E2fPx) +1.0*(E2fPx) -1.0*(E2fPx) -1.0*(E2fPx) +1.0*(E2fPx) -1.0*(E2fPx) +1.0*(E2fPx) -1.0*(E2fPx) +1.0*(E2fPx) +1.0*(E2f * Px[k+1]) -1.0*(E2fPx) +1.0*(E2fPx) ) * ( t_sim[k+1] - t_sim[k] ) )
@NLconstraint(m, [k in 1:length(t_sim)-1],
    Px[k+1] == Px[k] + ( +1.0*(E2fPx) -1.0*(E2f * Px[k+1]) ) * ( t_sim[k+1] - t_sim[k] ) )
@NLconstraint(m, [k in 1:length(t_sim)-1],
    Pcdc[k+1] == Pcdc[k] + ( -1.0*(pFox * Pcdc[k+1]) +1.0*(pFoxPcdc) ) * ( t_sim[k+1] - t_sim[k] ) )
@NLconstraint(m, [k in 1:length(t_sim)-1],
    pFoxPcdc[k+1] == pFoxPcdc[k] + ( -1.0*(pFoxPcdc) +1.0*(pFoxPcdc) +1.0*(pFox * Pcdc[k+1]) -1.0*(pFoxPcdc) ) * ( t_sim[k+1] - t_sim[k] ) )
@NLconstraint(m, [k in 1:length(t_sim)-1],
    Cb[k+1] == Cb[k] + ( -1.0*(E2fRb * Cb[k+1]) +1.0*(E2fRb * Cb[k+1]) -1.0*(E2f * Cb[k+1]) +1.0*(E2f * Cb[k+1]) -1.0*(pE2fRb * Cb[k+1]) +1.0*(pE2fRb * Cb[k+1]) -1.0*(E2fRb * Cb[k+1]) +1.0*(E2fRb * Cb[k+1]) -1.0*(Rb * Cb[k+1]) +1.0*(Rb * Cb[k+1]) -1.0*(Apc * Cb[k+1]) +1.0*(Apc * Cb[k+1]) -1.0*(Emi * Cb[k+1]) +1.0*(Emi * Cb[k+1]) -1.0*(Cdh * Cb[k+1]) +1.0*(Cdh * Cb[k+1]) -1.0*(Fox * Cb[k+1]) +1.0*(Fox * Cb[k+1]) +1.0*() +1.0*(pFoxPcb) -1.0*(Cb) -1.0*(Cb * ApcCdh[k+1]) -1.0*(Cb * pApcCdh[k+1]) -1.0*(Cb * pApcCdc[k+1]) -1.0*(Cb * Wee[k+1]) -1.0*(Cb * pWee[k+1]) +1.0*(pCb * pCdc25[k+1]) +1.0*(pCb * Cdc25[k+1]) -1.0*(Cdc25 * Cb[k+1]) +1.0*(Cdc25 * Cb[k+1]) -1.0*(Wee * Cb[k+1]) +1.0*(Wee * Cb[k+1]) -1.0*(Gw * Cb[k+1]) +1.0*(Gw * Cb[k+1]) ) * ( t_sim[k+1] - t_sim[k] ) )
@NLconstraint(m, [k in 1:length(t_sim)-1],
    pCb[k+1] == pCb[k] + ( +1.0*(Cb * Wee[k+1]) +1.0*(Cb * pWee[k+1]) -1.0*(pCb * pCdc25[k+1]) -1.0*(pCb * Cdc25[k+1]) -1.0*(pCb) -1.0*(pCb * ApcCdh[k+1]) -1.0*(pCb * pApcCdh[k+1]) -1.0*(pCb * pApcCdc[k+1]) ) * ( t_sim[k+1] - t_sim[k] ) )
@NLconstraint(m, [k in 1:length(t_sim)-1],
    Emi[k+1] == Emi[k] + ( +1.0*(ACE * Ce[k+1] * 0.8[k+1]) +1.0*(ACE * Ca[k+1] * 0.8[k+1]) +1.0*(ACE * Cb[k+1] * 0.8[k+1]) +1.0*(pACE * Ce[k+1] * 0.8[k+1]) +1.0*(pACE * Ca[k+1] * 0.8[k+1]) +1.0*(pACE * Cb[k+1] * 0.8[k+1]) -1.0*(Emi * pApcCdh - kDiACE[k+1] * pACE[k+1]) -1.0*(Emi * Ca[k+1]) -1.0*(Emi * Cb[k+1]) +1.0*(pEmi) +1.0*(E2fPx) -1.0*(ApcCdh * Emi - kDiACE[k+1] * ACE[k+1]) +1.0*() -1.0*(Emi) ) * ( t_sim[k+1] - t_sim[k] ) )
@NLconstraint(m, [k in 1:length(t_sim)-1],
    pEmi[k+1] == pEmi[k] + ( +1.0*(pACE * Ca[k+1]) +1.0*(pACE * Cb[k+1]) +1.0*(Emi * Ca[k+1]) +1.0*(Emi * Cb[k+1]) -1.0*(pEmi) +1.0*(ACE * Ca[k+1]) +1.0*(ACE * Cb[k+1]) -1.0*(pEmi) -1.0*(pEmi) ) * ( t_sim[k+1] - t_sim[k] ) )
@NLconstraint(m, [k in 1:length(t_sim)-1],
    Cdc20[k+1] == Cdc20[k] + ( +1.0*(pApcCdc * B55[k+1]) -1.0*(pApc * Cdc20 - kDipACdc[k+1] * pApcCdc[k+1]) +1.0*() +1.0*(pFoxPcdc) -1.0*(Cdc20) -1.0*(Cdc20 * ApcCdh[k+1]) -1.0*(Cdc20 * pApcCdh[k+1]) ) * ( t_sim[k+1] - t_sim[k] ) )
@NLconstraint(m, [k in 1:length(t_sim)-1],
    E2f[k+1] == E2f[k] + ( +1.0*() +1.0*(E2fPx) -1.0*(E2f) -1.0*(E2f * Rb - kDiE2fRb[k+1] * E2fRb[k+1]) +1.0*(E2fRb * Cd[k+1]) +1.0*(E2fRb * Ce[k+1]) +1.0*(E2fRb * Ca[k+1]) +1.0*(E2fRb * Cb[k+1]) -1.0*(E2f * Ca[k+1]) -1.0*(E2f * Cb[k+1]) +1.0*(pE2f) -1.0*(E2f * Px[k+1]) +1.0*(E2f * Px[k+1]) +1.0*(pE2f * B55[k+1]) ) * ( t_sim[k+1] - t_sim[k] ) )
@NLconstraint(m, [k in 1:length(t_sim)-1],
    pFox[k+1] == pFox[k] + ( -1.0*(pFox * Pcb[k+1]) +1.0*(pFox * Pcb[k+1]) -1.0*(pFox * pApcCdh[k+1]) -1.0*(pFox) -1.0*(pFox * ApcCdh[k+1]) +1.0*(Fox * Ce[k+1]) +1.0*(Fox * Ca[k+1]) +1.0*(Fox * Cb[k+1]) -1.0*(pFox) -1.0*(pFox * Pcdc[k+1]) +1.0*(pFox * Pcdc[k+1]) ) * ( t_sim[k+1] - t_sim[k] ) )
@NLconstraint(m, [k in 1:length(t_sim)-1],
    Ce[k+1] == Ce[k] + ( -1.0*(E2fRb * Ce[k+1]) +1.0*(E2fRb * Ce[k+1]) -1.0*(pE2fRb * Ce[k+1]) +1.0*(pE2fRb * Ce[k+1]) -1.0*(Rb * Ce[k+1]) +1.0*(Rb * Ce[k+1]) +1.0*() -1.0*(Ce) -2.0*(pow(Ce, 2)) +1.0*(pow(Ce, 2)) -1.0*(Cdh * Ce[k+1]) +1.0*(Cdh * Ce[k+1]) -1.0*(Ce * Ca[k+1]) -1.0*(Fox * Ce[k+1]) +1.0*(Fox * Ce[k+1]) +1.0*(E2fPx) ) * ( t_sim[k+1] - t_sim[k] ) )
@NLconstraint(m, [k in 1:length(t_sim)-1],
    Ca[k+1] == Ca[k] + ( -1.0*(E2fRb * Ca[k+1]) +1.0*(E2fRb * Ca[k+1]) -1.0*(E2f * Ca[k+1]) +1.0*(E2f * Ca[k+1]) -1.0*(pE2fRb * Ca[k+1]) +1.0*(pE2fRb * Ca[k+1]) -1.0*(E2fRb * Ca[k+1]) +1.0*(E2fRb * Ca[k+1]) -1.0*(Rb * Ca[k+1]) +1.0*(Rb * Ca[k+1]) +1.0*(E2fPx) -1.0*(Apc * Ca[k+1]) +1.0*(Apc * Ca[k+1]) -1.0*(Emi * Ca[k+1]) +1.0*(Emi * Ca[k+1]) -1.0*(Cdh * Ca[k+1]) +1.0*(Cdh * Ca[k+1]) +1.0*() -1.0*(Ca) -1.0*(Ce * Ca[k+1]) +1.0*(Ce * Ca[k+1]) -1.0*(Ca * ApcCdh[k+1]) -1.0*(Ca * pApcCdh[k+1]) -1.0*(Ca * pApcCdc[k+1]) -1.0*(Fox * Ca[k+1]) +1.0*(Fox * Ca[k+1]) -1.0*(Cdc25 * Ca[k+1]) +1.0*(Cdc25 * Ca[k+1]) -1.0*(Wee * Ca[k+1]) +1.0*(Wee * Ca[k+1]) ) * ( t_sim[k+1] - t_sim[k] ) )
@NLconstraint(m, [k in 1:length(t_sim)-1],
    pE2f[k+1] == pE2f[k] + ( +1.0*(E2f * Ca[k+1]) +1.0*(E2f * Cb[k+1]) -1.0*(pE2f) -1.0*(pE2f) -1.0*(pE2f) -1.0*(pE2f * Rb - kDiE2fRb[k+1] * pE2fRb[k+1]) +1.0*(pE2fRb * Cd[k+1]) +1.0*(pE2fRb * Ca[k+1]) +1.0*(pE2fRb * Cb[k+1]) +1.0*(pE2fRb * Ce[k+1]) -1.0*(pE2f * B55[k+1]) ) * ( t_sim[k+1] - t_sim[k] ) )
@NLconstraint(m, [k in 1:length(t_sim)-1],
    Fox[k+1] == Fox[k] + ( +1.0*() +1.0*(E2fPx) -1.0*(Fox) -1.0*(Fox * ApcCdh[k+1]) -1.0*(Fox * pApcCdh[k+1]) -1.0*(Fox * Ce[k+1]) -1.0*(Fox * Ca[k+1]) -1.0*(Fox * Cb[k+1]) +1.0*(pFox) ) * ( t_sim[k+1] - t_sim[k] ) )

# Define objective
println("Defining objective ...")
@variable(m, aux)
numerical_error_penalty = sum(weight*(Rb[k] - Rb[k-1])^pow for k in 2:length(t_sim))
    + sum(weight*(pRb[k] - pRb[k-1])^pow for k in 2:length(t_sim))
    + sum(weight*(pE2fRb[k] - pE2fRb[k-1])^pow for k in 2:length(t_sim))
    + sum(weight*(E2fRb[k] - E2fRb[k-1])^pow for k in 2:length(t_sim))
    + sum(weight*(Pcb[k] - Pcb[k-1])^pow for k in 2:length(t_sim))
    + sum(weight*(pFoxPcb[k] - pFoxPcb[k-1])^pow for k in 2:length(t_sim))
    + sum(weight*(ACE[k] - ACE[k-1])^pow for k in 2:length(t_sim))
    + sum(weight*(Apc[k] - Apc[k-1])^pow for k in 2:length(t_sim))
    + sum(weight*(ApcCdh[k] - ApcCdh[k-1])^pow for k in 2:length(t_sim))
    + sum(weight*(pACE[k] - pACE[k-1])^pow for k in 2:length(t_sim))
    + sum(weight*(pApc[k] - pApc[k-1])^pow for k in 2:length(t_sim))
    + sum(weight*(pApcCdc[k] - pApcCdc[k-1])^pow for k in 2:length(t_sim))
    + sum(weight*(pApcCdh[k] - pApcCdh[k-1])^pow for k in 2:length(t_sim))
    + sum(weight*(Cdh[k] - Cdh[k-1])^pow for k in 2:length(t_sim))
    + sum(weight*(pCdh[k] - pCdh[k-1])^pow for k in 2:length(t_sim))
    + sum(weight*(Ensa[k] - Ensa[k-1])^pow for k in 2:length(t_sim))
    + sum(weight*(pEB55[k] - pEB55[k-1])^pow for k in 2:length(t_sim))
    + sum(weight*(pEnsa[k] - pEnsa[k-1])^pow for k in 2:length(t_sim))
    + sum(weight*(B55[k] - B55[k-1])^pow for k in 2:length(t_sim))
    + sum(weight*(Gw[k] - Gw[k-1])^pow for k in 2:length(t_sim))
    + sum(weight*(pGw[k] - pGw[k-1])^pow for k in 2:length(t_sim))
    + sum(weight*(Cdc25[k] - Cdc25[k-1])^pow for k in 2:length(t_sim))
    + sum(weight*(pCdc25[k] - pCdc25[k-1])^pow for k in 2:length(t_sim))
    + sum(weight*(pWee[k] - pWee[k-1])^pow for k in 2:length(t_sim))
    + sum(weight*(Wee[k] - Wee[k-1])^pow for k in 2:length(t_sim))
    + sum(weight*(E2fPx[k] - E2fPx[k-1])^pow for k in 2:length(t_sim))
    + sum(weight*(Px[k] - Px[k-1])^pow for k in 2:length(t_sim))
    + sum(weight*(Pcdc[k] - Pcdc[k-1])^pow for k in 2:length(t_sim))
    + sum(weight*(pFoxPcdc[k] - pFoxPcdc[k-1])^pow for k in 2:length(t_sim))
    + sum(weight*(Cb[k] - Cb[k-1])^pow for k in 2:length(t_sim))
    + sum(weight*(pCb[k] - pCb[k-1])^pow for k in 2:length(t_sim))
    + sum(weight*(Emi[k] - Emi[k-1])^pow for k in 2:length(t_sim))
    + sum(weight*(pEmi[k] - pEmi[k-1])^pow for k in 2:length(t_sim))
    + sum(weight*(Cdc20[k] - Cdc20[k-1])^pow for k in 2:length(t_sim))
    + sum(weight*(E2f[k] - E2f[k-1])^pow for k in 2:length(t_sim))
    + sum(weight*(pFox[k] - pFox[k-1])^pow for k in 2:length(t_sim))
    + sum(weight*(Ce[k] - Ce[k-1])^pow for k in 2:length(t_sim))
    + sum(weight*(Ca[k] - Ca[k-1])^pow for k in 2:length(t_sim))
    + sum(weight*(pE2f[k] - pE2f[k-1])^pow for k in 2:length(t_sim))
    + sum(weight*(Fox[k] - Fox[k-1])^pow for k in 2:length(t_sim))
@constraint(m, aux == numerical_error_penalty)

@NLobjective(m, Min,sum((Rb[(k-1)*t_ratio+1]-df[k, :Rb])^2 for k in 1:length(t_exp))
    + sum((pRb[(k-1)*t_ratio+1]-df[k, :pRb])^2 for k in 1:length(t_exp))
    + sum((pE2fRb[(k-1)*t_ratio+1]-df[k, :pE2fRb])^2 for k in 1:length(t_exp))
    + sum((E2fRb[(k-1)*t_ratio+1]-df[k, :E2fRb])^2 for k in 1:length(t_exp))
    + sum((Pcb[(k-1)*t_ratio+1]-df[k, :Pcb])^2 for k in 1:length(t_exp))
    + sum((pFoxPcb[(k-1)*t_ratio+1]-df[k, :pFoxPcb])^2 for k in 1:length(t_exp))
    + sum((ACE[(k-1)*t_ratio+1]-df[k, :ACE])^2 for k in 1:length(t_exp))
    + sum((Apc[(k-1)*t_ratio+1]-df[k, :Apc])^2 for k in 1:length(t_exp))
    + sum((ApcCdh[(k-1)*t_ratio+1]-df[k, :ApcCdh])^2 for k in 1:length(t_exp))
    + sum((pACE[(k-1)*t_ratio+1]-df[k, :pACE])^2 for k in 1:length(t_exp))
    + sum((pApc[(k-1)*t_ratio+1]-df[k, :pApc])^2 for k in 1:length(t_exp))
    + sum((pApcCdc[(k-1)*t_ratio+1]-df[k, :pApcCdc])^2 for k in 1:length(t_exp))
    + sum((pApcCdh[(k-1)*t_ratio+1]-df[k, :pApcCdh])^2 for k in 1:length(t_exp))
    + sum((Cdh[(k-1)*t_ratio+1]-df[k, :Cdh])^2 for k in 1:length(t_exp))
    + sum((pCdh[(k-1)*t_ratio+1]-df[k, :pCdh])^2 for k in 1:length(t_exp))
    + sum((Ensa[(k-1)*t_ratio+1]-df[k, :Ensa])^2 for k in 1:length(t_exp))
    + sum((pEB55[(k-1)*t_ratio+1]-df[k, :pEB55])^2 for k in 1:length(t_exp))
    + sum((pEnsa[(k-1)*t_ratio+1]-df[k, :pEnsa])^2 for k in 1:length(t_exp))
    + sum((B55[(k-1)*t_ratio+1]-df[k, :B55])^2 for k in 1:length(t_exp))
    + sum((Gw[(k-1)*t_ratio+1]-df[k, :Gw])^2 for k in 1:length(t_exp))
    + sum((pGw[(k-1)*t_ratio+1]-df[k, :pGw])^2 for k in 1:length(t_exp))
    + sum((Cdc25[(k-1)*t_ratio+1]-df[k, :Cdc25])^2 for k in 1:length(t_exp))
    + sum((pCdc25[(k-1)*t_ratio+1]-df[k, :pCdc25])^2 for k in 1:length(t_exp))
    + sum((pWee[(k-1)*t_ratio+1]-df[k, :pWee])^2 for k in 1:length(t_exp))
    + sum((Wee[(k-1)*t_ratio+1]-df[k, :Wee])^2 for k in 1:length(t_exp))
    + sum((E2fPx[(k-1)*t_ratio+1]-df[k, :E2fPx])^2 for k in 1:length(t_exp))
    + sum((Px[(k-1)*t_ratio+1]-df[k, :Px])^2 for k in 1:length(t_exp))
    + sum((Pcdc[(k-1)*t_ratio+1]-df[k, :Pcdc])^2 for k in 1:length(t_exp))
    + sum((pFoxPcdc[(k-1)*t_ratio+1]-df[k, :pFoxPcdc])^2 for k in 1:length(t_exp))
    + sum((Cb[(k-1)*t_ratio+1]-df[k, :Cb])^2 for k in 1:length(t_exp))
    + sum((pCb[(k-1)*t_ratio+1]-df[k, :pCb])^2 for k in 1:length(t_exp))
    + sum((Emi[(k-1)*t_ratio+1]-df[k, :Emi])^2 for k in 1:length(t_exp))
    + sum((pEmi[(k-1)*t_ratio+1]-df[k, :pEmi])^2 for k in 1:length(t_exp))
    + sum((Cdc20[(k-1)*t_ratio+1]-df[k, :Cdc20])^2 for k in 1:length(t_exp))
    + sum((E2f[(k-1)*t_ratio+1]-df[k, :E2f])^2 for k in 1:length(t_exp))
    + sum((pFox[(k-1)*t_ratio+1]-df[k, :pFox])^2 for k in 1:length(t_exp))
    + sum((Ce[(k-1)*t_ratio+1]-df[k, :Ce])^2 for k in 1:length(t_exp))
    + sum((Ca[(k-1)*t_ratio+1]-df[k, :Ca])^2 for k in 1:length(t_exp))
    + sum((pE2f[(k-1)*t_ratio+1]-df[k, :pE2f])^2 for k in 1:length(t_exp))
    + sum((Fox[(k-1)*t_ratio+1]-df[k, :Fox])^2 for k in 1:length(t_exp))
    + aux
    )

println("Optimizing...")
optimize!(m)


println("# Obtain the solution")
println("Retreiving solution...")
species_to_plot = []
params = [tRb, tPcb, tApc, t_pApc, tCdh, tEnsa, tB55, tGw, tCdc25, tWee, tPx, tPcdc, tCb, tEmi, tCdc20, kSyE2f1, kAsFPcb, kSyE2f2, kDeE2f1, kDiEPx, kAsE2fRb, kDiE2fRb, kPhRbD, Cd, kPhRbE, kPhRb, kPhE2f, kDpE2f1, kDpE2f2, kDeE2f2, kSyCe1, kSyCa2, kDeCe, kPhCeE, kDpRb, kDpApc, kPhApcA, kPhApcB, kAsACdh, kDiACdh, kPhCdhA, kPhCdhE, kPhCdhB, kDeEmi1, kDeEmi2, kPhEmiA, kPhEmiB, kAsACE, kDiACE, kDpEmi, kDeEmi3, kSyFox1, kSyEmi2, kDpCdh, kSyCa1, kDeCa1, kPhCeA, kDeCa2, kDeCa3, kSyFox2, kDeFox1, kDeFox2, kPhFoxE, kPhFoxA, kPhFoxB, kDpFox, kSyCb1, kSyCb2, kAsEPx, kDiFPcb, kDeCb1, kDeCb2, kDeCb3, kAspACdc, kDipACdc, kDeCdc_1, kDeCdc_2, kSyCdc_1, kSyCdc_2, kAsFPcdc, kDiFPcdc, kWee2, kWee1, kCdc25_2, kCdc25_1, kPhC25B, kPhC25A, kDpCdc25, kDpWee, kPhWeeA, kPhWeeB, kPhGw, kDpGw1, kDpGw2, kPhEnsa, kAspEB55, kDipEB55, kDpEnsa, kSyEmi1, kSyCe2]
paramvalues = Dict()
for param in params
    paramvalues[param] = JuMP.value.(param)
end

variables = [ :Rb
 :pRb
 :pE2fRb
 :E2fRb
 :Pcb
 :pFoxPcb
 :ACE
 :Apc
 :ApcCdh
 :pACE
 :pApc
 :pApcCdc
 :pApcCdh
 :Cdh
 :pCdh
 :Ensa
 :pEB55
 :pEnsa
 :B55
 :Gw
 :pGw
 :Cdc25
 :pCdc25
 :pWee
 :Wee
 :E2fPx
 :Px
 :Pcdc
 :pFoxPcdc
 :Cb
 :pCb
 :Emi
 :pEmi
 :Cdc20
 :E2f
 :pFox
 :Ce
 :Ca
 :pE2f
 :Fox
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