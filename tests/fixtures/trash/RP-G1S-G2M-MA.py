from numpy import *
from matplotlib.pylab import *
from matplotlib.pyplot import *
from scipy.integrate import odeint 

LABELS = {'Ce': [], 'E2f': [], 'pE2f': [], 'Ca': [], 'Cb': [], 'E2fRb': [], 'pE2fRb': [], 'Rb': [], 'pRb': [], 'B55': [], 'E2fPx': [], 'Px': [], 'Apc': [], 'Cdc20': [], 'pApcCdc': [], 'pApc': [], 'ApcCdh': [], 'Cdh': [], 'pCdh': [], 'ACE': [], 'Emi': [], 'pACE': [], 'pApcCdh': [], 'pEmi': [], 'Fox': [], 'pFox': [], 'pFoxPcb': [], 'Pcb': [], 'Pcdc': [], 'pFoxPcdc': [], 'Wee': [], 'pCb': [], 'pWee': [], 'pCdc25': [], 'Cdc25': [], 'Gw': [], 'pGw': [], 'Ensa': [], 'pEnsa': [], 'pEB55': []}

def simulateModel(t0, tend, numPoints):
  
  #compartments
  compartment = 1.0
  
  #global parameters
  kPhRbD = 0.8
  kPhRbE = 1.0
  kDpRb = 0.2
  kPhRb = 3.0
  kSyCe1 = 0.002
  kSyCe2 = 0.06
  kDeCe = 0.02
  kDiE2fRb = 0.05
  kAsE2fRb = 50.0
  kSyE2f1 = 0.001
  kSyE2f2 = 0.01
  kDeE2f1 = 0.01
  kDeE2f2 = 0.2
  kPhE2f = 0.2
  kDpE2f1 = 0.0
  kDpE2f2 = 10.0
  kSyEmi1 = 0.005
  kSyEmi2 = 0.05
  kDeEmi1 = 0.05
  kDeEmi2 = 0.05
  kDeEmi3 = 0.5
  kPhEmiA = 0.005
  kPhEmiB = 0.025
  kDpEmi = 0.0005
  kSyCa1 = 0.001
  kSyCa2 = 0.0075
  kDeCa1 = 0.001
  kDeCa2 = 0.4
  kDeCa3 = 0.1
  kDpCdh = 0.5
  kPhCdhA = 2.5
  kPhCdhE = 0.1
  kPhCdhB = 25.0
  kDiACE = 0.05
  kAsACE = 50.0
  kDiACdh = 0.05
  kAsACdh = 50.0
  kPhCeA = 1.0
  kPhCeE = 0.01
  kPhFoxE = 0.1
  kPhFoxA = 0.1
  kPhFoxB = 0.5
  kDpFox = 0.005
  kSyCb1 = 0.001
  kSyCb2 = 0.0075
  kDeCb1 = 0.001
  kDeCb2 = 0.125
  kDeCb3 = 0.125
  kPhEnsa = 0.1
  kDpEnsa = 0.05
  kPhGw = 1.0
  kDpGw1 = 0.25
  kDpGw2 = 10.0
  kWee1 = 0.01
  kWee2 = 1.0
  kPhWeeA = 0.1
  kPhWeeB = 1.0
  kDpWee = 10.0
  kCdc25_1 = 0.1
  kCdc25_2 = 1.0
  kPhC25A = 0.0
  kDpCdc25 = 10.0
  kDipEB55 = 0.0068
  kAspEB55 = 57.0
  kSyCdc_1 = 0.001
  kSyCdc_2 = 0.01
  kDeCdc_1 = 0.004
  kDeCdc_2 = 0.04
  kDipACdc = 0.125
  kAspACdc = 1.0
  kPhApcA = 0.01
  kPhApcB = 0.01
  kDpApc = 5.0
  kSyFox1 = 0.001
  kSyFox2 = 0.0075
  kDeFox1 = 0.001
  kDeFox2 = 0.5
  kAsEPx = 50.0
  kDiEPx = 5.0
  kAsFPcb = 50.0
  kDiFPcb = 10.0
  kAsFPcdc = 50.0
  kDiFPcdc = 25.0
  f1 = 1.0
  f2 = 1.0
  f3 = 1.0
  f4 = 0.8
  kPhC25B = 1.0
  
  #boundary species
  Cd = 1.0
  tRb = 1.0
  tPcb = 1.0
  tApc = 5.0
  t_pApc = 0.0
  tCdh = 1.0
  tEnsa = 1.0
  tB55 = 0.25
  tGw = 1.0
  tCdc25 = 1.0
  tWee = 1.0
  tPx = 1.0
  tPcdc = 1.0
  tCb = 0.0
  tEmi = 0.0
  tCdc20 = 0.0
  
  def ode_fun(__Y__, t):
    Ce = __Y__[0]
    E2f = __Y__[1]
    pE2f = __Y__[2]
    Ca = __Y__[3]
    Cb = __Y__[4]
    E2fRb = __Y__[5]
    pE2fRb = __Y__[6]
    Rb = __Y__[7]
    pRb = __Y__[8]
    B55 = __Y__[9]
    E2fPx = __Y__[10]
    Px = __Y__[11]
    Apc = __Y__[12]
    Cdc20 = __Y__[13]
    pApcCdc = __Y__[14]
    pApc = __Y__[15]
    ApcCdh = __Y__[16]
    Cdh = __Y__[17]
    pCdh = __Y__[18]
    ACE = __Y__[19]
    Emi = __Y__[20]
    pACE = __Y__[21]
    pApcCdh = __Y__[22]
    pEmi = __Y__[23]
    Fox = __Y__[24]
    pFox = __Y__[25]
    pFoxPcb = __Y__[26]
    Pcb = __Y__[27]
    Pcdc = __Y__[28]
    pFoxPcdc = __Y__[29]
    Wee = __Y__[30]
    pCb = __Y__[31]
    pWee = __Y__[32]
    pCdc25 = __Y__[33]
    Cdc25 = __Y__[34]
    Gw = __Y__[35]
    pGw = __Y__[36]
    Ensa = __Y__[37]
    pEnsa = __Y__[38]
    pEB55 = __Y__[39]

    SyE2f1 = compartment * kSyE2f1
    AspFoxPcb = compartment * kAsFPcb * pFox * Pcb
    SyE2f2 = compartment * kSyE2f2 * E2fPx
    DeE2f = compartment * kDeE2f1 * E2f
    DiE2fPx = compartment * kDiEPx * E2fPx
    AsE2fRb = compartment * (kAsE2fRb * E2f * Rb - kDiE2fRb * E2fRb)
    PhRbE2fByCd = compartment * kPhRbD * E2fRb * Cd
    PhRbE2fByCe = compartment * kPhRbE * E2fRb * Ce
    PhRbE2fByCa = compartment * kPhRb * E2fRb * Ca
    PhRbE2fByCb = compartment * kPhRb * E2fRb * Cb
    PhE2fByCa = compartment * kPhE2f * E2f * Ca
    PhE2fByCb = compartment * kPhE2f * E2f * Cb
    DppE2f1 = compartment * kDpE2f1 * pE2f
    DppE2fRb2 = compartment * kDpE2f2 * pE2fRb * B55
    DepE2f1 = compartment * kDeE2f1 * pE2f
    DepE2f2 = compartment * kDeE2f2 * pE2f
    AspE2fRb = compartment * (kAsE2fRb * pE2f * Rb - kDiE2fRb * pE2fRb)
    PhRbpE2fByCd = compartment * kPhRbD * pE2fRb * Cd
    PhRbpE2fByCa = compartment * kPhRb * pE2fRb * Ca
    PhRbpE2fByCb = compartment * kPhRb * pE2fRb * Cb
    PhRbpE2fByCe = compartment * kPhRbE * pE2fRb * Ce
    PhE2fRbByCa = compartment * kPhE2f * E2fRb * Ca
    PhE2fRbByCb = compartment * kPhE2f * E2fRb * Cb
    DppE2fRb1 = compartment * kDpE2f1 * pE2fRb
    DeE2fRb = compartment * kDeE2f1 * E2fRb
    DepE2fRb1 = compartment * kDeE2f1 * pE2fRb
    DepE2fRb2 = compartment * kDeE2f2 * pE2fRb
    PhRbByCd = compartment * kPhRbD * Rb * Cd
    PhRbByCe = compartment * kPhRbE * Rb * Ce
    PhRbByCa = compartment * kPhRb * Rb * Ca
    PhRbByCb = compartment * kPhRb * Rb * Cb
    SyCe1 = compartment * kSyCe1
    SyCa2 = compartment * kSyCa2 * E2fPx
    DeCe1 = compartment * kDeCe * Ce
    DeCe2 = compartment * kPhCeE * pow(Ce, 2)
    DpRb = compartment * kDpRb * pRb
    DppApcCdc20 = compartment * kDpApc * f1 * pApcCdc * B55
    PhApcByCa = compartment * kPhApcA * Apc * Ca
    PhApcByCb = compartment * kPhApcB * Apc * Cb
    DppApc = compartment * kDpApc * pApc * B55
    AsApcCdh = compartment * (kAsACdh * Cdh * Apc - kDiACdh * ApcCdh)
    PhCdhApcByCa = compartment * kPhCdhA * f2 * ApcCdh * Ca
    PhCAEByCe = compartment * kPhCdhE * f4 * ACE * Ce
    PhCAEByCa = compartment * kPhCdhA * f4 * ACE * Ca
    PhCAEByCb = compartment * kPhCdhB * f4 * ACE * Cb
    PhCpAEByCe = compartment * kPhCdhE * f4 * pACE * Ce
    PhCpAEByCa = compartment * kPhCdhA * f4 * pACE * Ca
    PhCpAEByCb = compartment * kPhCdhB * f4 * pACE * Cb
    PhCdhpApcByCe = compartment * kPhCdhE * f2 * pApcCdh * Ce
    DppApcCdh = compartment * kDpApc * f1 * pApcCdh * B55
    AspApcCdh = compartment * (kAsACdh * pApc * Cdh - kDiACdh * pApcCdh)
    DeEpAC1 = compartment * kDeEmi1 * pACE
    DeEpAC2 = compartment * kDeEmi2 * pACE
    PhEpACByCa = compartment * kPhEmiA * f3 * pACE * Ca
    PhEpACByCb = compartment * kPhEmiB * f3 * pACE * Cb
    AspACE = compartment * (kAsACE * Emi * pApcCdh - kDiACE * pACE)
    PhACEByCa = compartment * kPhApcA * f1 * ACE * Ca
    PhACEByCb = compartment * kPhApcB * f1 * ACE * Cb
    DppACE = compartment * kDpApc * f1 * pACE * B55
    PhEmiByCa = compartment * kPhEmiA * Emi * Ca
    PhEmiByCb = compartment * kPhEmiB * Emi * Cb
    DppEmi = compartment * kDpEmi * pEmi
    PhEACByCa = compartment * kPhEmiA * f3 * ACE * Ca
    PhEACByCb = compartment * kPhEmiB * f3 * ACE * Cb
    DepEmi1 = compartment * kDeEmi1 * pEmi
    DepEmi2 = compartment * kDeEmi3 * pEmi
    SyFox1 = compartment * kSyFox1
    SyEmi2 = compartment * kSyEmi2 * E2fPx
    AsACE = compartment * (kAsACE * ApcCdh * Emi - kDiACE * ACE)
    PhCdhByCe = compartment * kPhCdhE * Cdh * Ce
    PhCdhByCa = compartment * kPhCdhA * Cdh * Ca
    PhCdhByCb = compartment * kPhCdhB * Cdh * Cb
    DppCdh = compartment * kDpCdh * pCdh
    SyCa1 = compartment * kSyCa1
    DeCa1 = compartment * kDeCa1 * Ca
    DeCe3 = compartment * kPhCeA * Ce * Ca
    DeCa2_1 = compartment * kDeCa2 * Ca * ApcCdh
    DeCa2_2 = compartment * kDeCa2 * Ca * pApcCdh
    DeCa3 = compartment * kDeCa3 * Ca * pApcCdc
    SyFox2 = compartment * kSyFox2 * E2fPx
    DeFox1 = compartment * kDeFox1 * Fox
    DeFox2_1 = compartment * kDeFox2 * Fox * ApcCdh
    DepFox2_2 = compartment * kDeFox2 * pFox * pApcCdh
    DepFox1 = compartment * kDeFox1 * pFox
    DepFox2_1 = compartment * kDeFox2 * pFox * ApcCdh
    DeFox2_2 = compartment * kDeFox2 * Fox * pApcCdh
    PhFoxByCe = compartment * kPhFoxE * Fox * Ce
    PhFoxByCa = compartment * kPhFoxA * Fox * Ca
    PhFoxByCb = compartment * kPhFoxB * Fox * Cb
    DppFox = compartment * kDpFox * pFox
    SyCb1 = compartment * kSyCb1
    SyCb2 = compartment * kSyCb2 * pFoxPcb
    AsE2fPx = compartment * kAsEPx * E2f * Px
    DipFoxPcb = compartment * kDiFPcb * pFoxPcb
    DeCb1 = compartment * kDeCb1 * Cb
    DeCb2_1 = compartment * kDeCb2 * Cb * ApcCdh
    DeCb2_2 = compartment * kDeCb2 * Cb * pApcCdh
    DeCb3 = compartment * kDeCb3 * Cb * pApcCdc
    AspApcCdc20 = compartment * (kAspACdc * pApc * Cdc20 - kDipACdc * pApcCdc)
    DeCdcpApc1 = compartment * kDeCdc_1 * pApcCdc
    DeCdcpApc2_1 = compartment * kDeCdc_2 * pApcCdc * ApcCdh
    DeCdcpApc2_2 = compartment * kDeCdc_2 * pApcCdc * pApcCdh
    SyCdc20_1 = compartment * kSyCdc_1
    SyCdc20_2 = compartment * kSyCdc_2 * pFoxPcdc
    AspFoxPcdc = compartment * kAsFPcdc * pFox * Pcdc
    DipFoxPcdc = compartment * kDiFPcdc * pFoxPcdc
    PhCbByWee = compartment * kWee2 * Cb * Wee
    PhCbBypWee = compartment * kWee1 * Cb * pWee
    DppCbBypCdc25 = compartment * kCdc25_2 * pCb * pCdc25
    DppCbByCdc25 = compartment * kCdc25_1 * pCb * Cdc25
    DepCb1 = compartment * kDeCb1 * pCb
    DepCb2_1 = compartment * kDeCb2 * pCb * ApcCdh
    DepCb2_2 = compartment * kDeCb2 * pCb * pApcCdh
    DepCb3 = compartment * kDeCb3 * pCb * pApcCdc
    PhCdc25ByCb = compartment * kPhC25B * Cdc25 * Cb
    PhCdc25ByCa = compartment * kPhC25A * Cdc25 * Ca
    DppCdc25 = compartment * kDpCdc25 * pCdc25 * B55
    DppWee = compartment * kDpWee * pWee * B55
    PhWeeByCa = compartment * kPhWeeA * Wee * Ca
    PhWeeByCb = compartment * kPhWeeB * Wee * Cb
    PhGw = compartment * kPhGw * Gw * Cb
    DppGw1 = compartment * kDpGw1 * pGw
    DppGw2 = compartment * kDpGw2 * pGw * B55
    PhEnsa = compartment * kPhEnsa * Ensa * pGw
    AspEB55 = compartment * (kAspEB55 * pEnsa * B55 - kDipEB55 * pEB55)
    DppEB55 = compartment * kDpEnsa * pEB55
    SyEmi1 = compartment * kSyEmi1
    SyCe2 = compartment * kSyCe2 * E2fPx
    DeCdc20_1 = compartment * kDeCdc_1 * Cdc20
    DeCdc20_2_1 = compartment * kDeCdc_2 * Cdc20 * ApcCdh
    DeCdc20_2_2 = compartment * kDeCdc_2 * Cdc20 * pApcCdh
    DppE2f2 = compartment * kDpE2f2 * pE2f * B55
    PhCdhpApcByCa = compartment * kPhCdhA * f2 * pApcCdh * Ca
    PhCdhpApcByCb = compartment * kPhCdhB * f2 * pApcCdh * Cb
    PhCdhApcByCe = compartment * kPhCdhE * f2 * ApcCdh * Ce
    PhCdhApcByCb = compartment * kPhCdhB * f2 * ApcCdh * Cb
    DeEAC1 = compartment * kDeEmi1 * ACE
    DeEAC2 = compartment * kDeEmi2 * ACE
    PhApcCdhByCa = compartment * kPhApcA * f1 * ApcCdh * Ca
    PhApcCdhByCb = compartment * kPhApcB * f1 * ApcCdh * Cb
    DeEmi1 = compartment * kDeEmi1 * Emi

    return array([ + (-PhRbE2fByCe) + (PhRbE2fByCe) + (-PhRbpE2fByCe) + (PhRbpE2fByCe) + (-PhRbByCe) + (PhRbByCe) + (SyCe1) + (-DeCe1) + (-(2.0)*DeCe2) + (DeCe2) + (-PhCdhByCe) + (PhCdhByCe) + (-DeCe3) + (-PhFoxByCe) + (PhFoxByCe) + (SyCe2),
       + (SyE2f1) + (SyE2f2) + (-DeE2f) + (-AsE2fRb) + (PhRbE2fByCd) + (PhRbE2fByCe) + (PhRbE2fByCa) + (PhRbE2fByCb) + (-PhE2fByCa) + (-PhE2fByCb) + (DppE2f1) + (-AsE2fPx) + (AsE2fPx) + (DppE2f2),
       + (PhE2fByCa) + (PhE2fByCb) + (-DppE2f1) + (-DepE2f1) + (-DepE2f2) + (-AspE2fRb) + (PhRbpE2fByCd) + (PhRbpE2fByCa) + (PhRbpE2fByCb) + (PhRbpE2fByCe) + (-DppE2f2),
       + (-PhRbE2fByCa) + (PhRbE2fByCa) + (-PhE2fByCa) + (PhE2fByCa) + (-PhRbpE2fByCa) + (PhRbpE2fByCa) + (-PhE2fRbByCa) + (PhE2fRbByCa) + (-PhRbByCa) + (PhRbByCa) + (SyCa2) + (-PhApcByCa) + (PhApcByCa) + (-PhEmiByCa) + (PhEmiByCa) + (-PhCdhByCa) + (PhCdhByCa) + (SyCa1) + (-DeCa1) + (-DeCe3) + (DeCe3) + (-DeCa2_1) + (-DeCa2_2) + (-DeCa3) + (-PhFoxByCa) + (PhFoxByCa) + (-PhCdc25ByCa) + (PhCdc25ByCa) + (-PhWeeByCa) + (PhWeeByCa),
       + (-PhRbE2fByCb) + (PhRbE2fByCb) + (-PhE2fByCb) + (PhE2fByCb) + (-PhRbpE2fByCb) + (PhRbpE2fByCb) + (-PhE2fRbByCb) + (PhE2fRbByCb) + (-PhRbByCb) + (PhRbByCb) + (-PhApcByCb) + (PhApcByCb) + (-PhEmiByCb) + (PhEmiByCb) + (-PhCdhByCb) + (PhCdhByCb) + (-PhFoxByCb) + (PhFoxByCb) + (SyCb1) + (SyCb2) + (-DeCb1) + (-DeCb2_1) + (-DeCb2_2) + (-DeCb3) + (-PhCbByWee) + (-PhCbBypWee) + (DppCbBypCdc25) + (DppCbByCdc25) + (-PhCdc25ByCb) + (PhCdc25ByCb) + (-PhWeeByCb) + (PhWeeByCb) + (-PhGw) + (PhGw),
       + (AsE2fRb) + (-PhRbE2fByCd) + (-PhRbE2fByCe) + (-PhRbE2fByCa) + (-PhRbE2fByCb) + (DppE2fRb2) + (-PhE2fRbByCa) + (-PhE2fRbByCb) + (DppE2fRb1) + (-DeE2fRb),
       + (-DppE2fRb2) + (AspE2fRb) + (-PhRbpE2fByCd) + (-PhRbpE2fByCa) + (-PhRbpE2fByCb) + (-PhRbpE2fByCe) + (PhE2fRbByCa) + (PhE2fRbByCb) + (-DppE2fRb1) + (-DepE2fRb1) + (-DepE2fRb2),
       + (-AsE2fRb) + (-AspE2fRb) + (DeE2fRb) + (DepE2fRb1) + (DepE2fRb2) + (-PhRbByCd) + (-PhRbByCe) + (-PhRbByCa) + (-PhRbByCb) + (DpRb),
       + (PhRbE2fByCd) + (PhRbE2fByCe) + (PhRbE2fByCa) + (PhRbE2fByCb) + (PhRbpE2fByCd) + (PhRbpE2fByCa) + (PhRbpE2fByCb) + (PhRbpE2fByCe) + (PhRbByCd) + (PhRbByCe) + (PhRbByCa) + (PhRbByCb) + (-DpRb),
       + (-DppE2fRb2) + (DppE2fRb2) + (-DppApc) + (DppApc) + (-DppCdc25) + (DppCdc25) + (-DppWee) + (DppWee) + (-DppGw2) + (DppGw2) + (-AspEB55) + (DppEB55) + (-DppE2f2) + (DppE2f2),
       + (-SyE2f2) + (SyE2f2) + (-DiE2fPx) + (-SyCa2) + (SyCa2) + (-SyEmi2) + (SyEmi2) + (-SyFox2) + (SyFox2) + (AsE2fPx) + (-SyCe2) + (SyCe2),
       + (DiE2fPx) + (-AsE2fPx),
       + (DppApcCdc20) + (-PhApcByCa) + (-PhApcByCb) + (DppApc) + (-AsApcCdh) + (PhCdhApcByCa) + (PhCAEByCe) + (PhCAEByCa) + (PhCAEByCb) + (PhCdhApcByCe) + (PhCdhApcByCb),
       + (DppApcCdc20) + (-AspApcCdc20) + (SyCdc20_1) + (SyCdc20_2) + (-DeCdc20_1) + (-DeCdc20_2_1) + (-DeCdc20_2_2),
       + (-DppApcCdc20) + (-DeCa3) + (DeCa3) + (-DeCb3) + (DeCb3) + (AspApcCdc20) + (-DeCdcpApc1) + (-DeCdcpApc2_1) + (-DeCdcpApc2_2) + (-DepCb3) + (DepCb3),
       + (PhApcByCa) + (PhApcByCb) + (-DppApc) + (PhCpAEByCe) + (PhCpAEByCa) + (PhCpAEByCb) + (PhCdhpApcByCe) + (-AspApcCdh) + (-AspApcCdc20) + (DeCdcpApc1) + (DeCdcpApc2_1) + (DeCdcpApc2_2) + (PhCdhpApcByCa) + (PhCdhpApcByCb),
       + (AsApcCdh) + (-PhCdhApcByCa) + (DppApcCdh) + (PhEACByCa) + (PhEACByCb) + (-AsACE) + (-DeCa2_1) + (DeCa2_1) + (-DeFox2_1) + (DeFox2_1) + (-DepFox2_1) + (DepFox2_1) + (-DeCb2_1) + (DeCb2_1) + (-DeCdcpApc2_1) + (DeCdcpApc2_1) + (-DepCb2_1) + (DepCb2_1) + (-DeCdc20_2_1) + (DeCdc20_2_1) + (-PhCdhApcByCe) + (-PhCdhApcByCb) + (DeEAC1) + (DeEAC2) + (-PhApcCdhByCa) + (-PhApcCdhByCb),
       + (-AsApcCdh) + (-AspApcCdh) + (-PhCdhByCe) + (-PhCdhByCa) + (-PhCdhByCb) + (DppCdh),
       + (PhCdhApcByCa) + (PhCAEByCe) + (PhCAEByCa) + (PhCAEByCb) + (PhCpAEByCe) + (PhCpAEByCa) + (PhCpAEByCb) + (PhCdhpApcByCe) + (PhCdhByCe) + (PhCdhByCa) + (PhCdhByCb) + (-DppCdh) + (PhCdhpApcByCa) + (PhCdhpApcByCb) + (PhCdhApcByCe) + (PhCdhApcByCb),
       + (-PhCAEByCe) + (-PhCAEByCa) + (-PhCAEByCb) + (-PhACEByCa) + (-PhACEByCb) + (DppACE) + (-PhEACByCa) + (-PhEACByCb) + (AsACE) + (-DeEAC1) + (-DeEAC2),
       + (PhCAEByCe) + (PhCAEByCa) + (PhCAEByCb) + (PhCpAEByCe) + (PhCpAEByCa) + (PhCpAEByCb) + (-AspACE) + (-PhEmiByCa) + (-PhEmiByCb) + (DppEmi) + (SyEmi2) + (-AsACE) + (SyEmi1) + (-DeEmi1),
       + (-PhCpAEByCe) + (-PhCpAEByCa) + (-PhCpAEByCb) + (-DeEpAC1) + (-DeEpAC2) + (-PhEpACByCa) + (-PhEpACByCb) + (AspACE) + (PhACEByCa) + (PhACEByCb) + (-DppACE),
       + (-PhCdhpApcByCe) + (-DppApcCdh) + (AspApcCdh) + (DeEpAC1) + (DeEpAC2) + (PhEpACByCa) + (PhEpACByCb) + (-AspACE) + (-DeCa2_2) + (DeCa2_2) + (-DepFox2_2) + (DepFox2_2) + (-DeFox2_2) + (DeFox2_2) + (-DeCb2_2) + (DeCb2_2) + (-DeCdcpApc2_2) + (DeCdcpApc2_2) + (-DepCb2_2) + (DepCb2_2) + (-DeCdc20_2_2) + (DeCdc20_2_2) + (-PhCdhpApcByCa) + (-PhCdhpApcByCb) + (PhApcCdhByCa) + (PhApcCdhByCb),
       + (PhEpACByCa) + (PhEpACByCb) + (PhEmiByCa) + (PhEmiByCb) + (-DppEmi) + (PhEACByCa) + (PhEACByCb) + (-DepEmi1) + (-DepEmi2),
       + (SyFox1) + (SyFox2) + (-DeFox1) + (-DeFox2_1) + (-DeFox2_2) + (-PhFoxByCe) + (-PhFoxByCa) + (-PhFoxByCb) + (DppFox),
       + (-AspFoxPcb) + (AspFoxPcb) + (-DepFox2_2) + (-DepFox1) + (-DepFox2_1) + (PhFoxByCe) + (PhFoxByCa) + (PhFoxByCb) + (-DppFox) + (-AspFoxPcdc) + (AspFoxPcdc),
       + (AspFoxPcb) + (-SyCb2) + (SyCb2) + (-DipFoxPcb),
       + (-AspFoxPcb) + (DipFoxPcb),
       + (-AspFoxPcdc) + (DipFoxPcdc),
       + (-SyCdc20_2) + (SyCdc20_2) + (AspFoxPcdc) + (-DipFoxPcdc),
       + (-PhCbByWee) + (PhCbByWee) + (DppWee) + (-PhWeeByCa) + (-PhWeeByCb),
       + (PhCbByWee) + (PhCbBypWee) + (-DppCbBypCdc25) + (-DppCbByCdc25) + (-DepCb1) + (-DepCb2_1) + (-DepCb2_2) + (-DepCb3),
       + (-PhCbBypWee) + (PhCbBypWee) + (-DppWee) + (PhWeeByCa) + (PhWeeByCb),
       + (-DppCbBypCdc25) + (DppCbBypCdc25) + (PhCdc25ByCb) + (PhCdc25ByCa) + (-DppCdc25),
       + (-DppCbByCdc25) + (DppCbByCdc25) + (-PhCdc25ByCb) + (-PhCdc25ByCa) + (DppCdc25),
       + (-PhGw) + (DppGw1) + (DppGw2),
       + (PhGw) + (-DppGw1) + (-DppGw2) + (-PhEnsa) + (PhEnsa),
       + (-PhEnsa) + (DppEB55),
       + (PhEnsa) + (-AspEB55),
       + (AspEB55) + (-DppEB55)    ])

  time = linspace(t0, tend, numPoints)
  yinit= array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.25, 0.0, 1.0, 5.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 1.0, 0.0, 0.0])
  
  y = odeint(ode_fun, yinit, time)

  return time, y


time, result = simulateModel(0, 1000, 1001)

fig = figure()
ax = subplot(111)
plot(time,result[:,0], label="Ce", lw=1.5)
plot(time,result[:,1], label="E2f", lw=1.5)
plot(time,result[:,2], label="pE2f", lw=1.5)
plot(time,result[:,3], label="Ca", lw=1.5)
plot(time,result[:,4], label="Cb", lw=1.5)
plot(time,result[:,5], label="E2fRb", lw=1.5)
plot(time,result[:,6], label="pE2fRb", lw=1.5)
plot(time,result[:,7], label="Rb", lw=1.5)
plot(time,result[:,8], label="pRb", lw=1.5)
plot(time,result[:,9], label="B55", lw=1.5)
plot(time,result[:,10], label="E2fPx", lw=1.5)
plot(time,result[:,11], label="Px", lw=1.5)
plot(time,result[:,12], label="Apc", lw=1.5)
plot(time,result[:,13], label="Cdc20", lw=1.5)
plot(time,result[:,14], label="pApcCdc", lw=1.5)
plot(time,result[:,15], label="pApc", lw=1.5)
plot(time,result[:,16], label="ApcCdh", lw=1.5)
plot(time,result[:,17], label="Cdh", lw=1.5)
plot(time,result[:,18], label="pCdh", lw=1.5)
plot(time,result[:,19], label="ACE", lw=1.5)
plot(time,result[:,20], label="Emi", lw=1.5)
plot(time,result[:,21], label="pACE", lw=1.5)
plot(time,result[:,22], label="pApcCdh", lw=1.5)
plot(time,result[:,23], label="pEmi", lw=1.5)
plot(time,result[:,24], label="Fox", lw=1.5)
plot(time,result[:,25], label="pFox", lw=1.5)
plot(time,result[:,26], label="pFoxPcb", lw=1.5)
plot(time,result[:,27], label="Pcb", lw=1.5)
plot(time,result[:,28], label="Pcdc", lw=1.5)
plot(time,result[:,29], label="pFoxPcdc", lw=1.5)
plot(time,result[:,30], label="Wee", lw=1.5)
plot(time,result[:,31], label="pCb", lw=1.5)
plot(time,result[:,32], label="pWee", lw=1.5)
plot(time,result[:,33], label="pCdc25", lw=1.5)
plot(time,result[:,34], label="Cdc25", lw=1.5)
plot(time,result[:,35], label="Gw", lw=1.5)
plot(time,result[:,36], label="pGw", lw=1.5)
plot(time,result[:,37], label="Ensa", lw=1.5)
plot(time,result[:,38], label="pEnsa", lw=1.5)
plot(time,result[:,39], label="pEB55", lw=1.5)
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.7, box.height])
xlabel("time")
ylabel("concentration")
legend(loc="center left", bbox_to_anchor=(1, 0.5))
show()
