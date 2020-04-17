from numpy import *
from matplotlib.pylab import *
from matplotlib.pyplot import *
from scipy.integrate import odeint 

class PythonModel(object):
  def __init__(self):
        self._params = [0.8, 1.0, 0.2, 3.0, 0.002, 0.06, 0.02, 0.05, 50.0, 0.001, 0.01, 0.01, 0.2, 0.2, 0.0, 10.0, 0.005, 0.05, 0.05, 0.05, 0.5, 0.005, 0.025, 0.0005, 0.001, 0.0075, 0.001, 0.4, 0.1, 0.5, 2.5, 0.1, 25.0, 0.05, 50.0, 0.05, 50.0, 1.0, 0.01, 0.1, 0.1, 0.5, 0.005, 0.001, 0.0075, 0.001, 0.125, 0.125, 0.1, 0.05, 1.0, 0.25, 10.0, 0.01, 1.0, 0.1, 1.0, 10.0, 0.1, 1.0, 0.0, 10.0, 0.0068, 57.0, 0.001, 0.01, 0.004, 0.04, 0.125, 1.0, 0.01, 0.01, 5.0, 0.001, 0.0075, 0.001, 0.5, 50.0, 5.0, 50.0, 10.0, 50.0, 25.0, 1.0, 1.0, 1.0, 0.8, 1.0]
        self._x_0 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.25, 0.0, 1.0, 5.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 1.0, 0.0, 0.0]
        self._labels = ['Ce', 'E2f', 'pE2f', 'Ca', 'Cb', 'E2fRb', 'pE2fRb', 'Rb', 'pRb', 'B55', 'E2fPx', 'Px', 'Apc', 'Cdc20', 'pApcCdc', 'pApc', 'ApcCdh', 'Cdh', 'pCdh', 'ACE', 'Emi', 'pACE', 'pApcCdh', 'pEmi', 'Fox', 'pFox', 'pFoxPcb', 'Pcb', 'Pcdc', 'pFoxPcdc', 'Wee', 'pCb', 'pWee', 'pCdc25', 'Cdc25', 'Gw', 'pGw', 'Ensa', 'pEnsa', 'pEB55']
        self.simulate()

  def dx_dt(self, x, t, par='NOTHING'):
  
  #compartments
    compartment = 1.0
    
    #global parameters
    if par == 'NOTHING':
      par = self.params
    
    Ce = x[0]
    E2f = x[1]
    pE2f = x[2]
    Ca = x[3]
    Cb = x[4]
    E2fRb = x[5]
    pE2fRb = x[6]
    Rb = x[7]
    pRb = x[8]
    B55 = x[9]
    E2fPx = x[10]
    Px = x[11]
    Apc = x[12]
    Cdc20 = x[13]
    pApcCdc = x[14]
    pApc = x[15]
    ApcCdh = x[16]
    Cdh = x[17]
    pCdh = x[18]
    ACE = x[19]
    Emi = x[20]
    pACE = x[21]
    pApcCdh = x[22]
    pEmi = x[23]
    Fox = x[24]
    pFox = x[25]
    pFoxPcb = x[26]
    Pcb = x[27]
    Pcdc = x[28]
    pFoxPcdc = x[29]
    Wee = x[30]
    pCb = x[31]
    pWee = x[32]
    pCdc25 = x[33]
    Cdc25 = x[34]
    Gw = x[35]
    pGw = x[36]
    Ensa = x[37]
    pEnsa = x[38]
    pEB55 = x[39]
  
    kPhRbD = par[0]
    kPhRbE = par[1]
    kDpRb = par[2]
    kPhRb = par[3]
    kSyCe1 = par[4]
    kSyCe2 = par[5]
    kDeCe = par[6]
    kDiE2fRb = par[7]
    kAsE2fRb = par[8]
    kSyE2f1 = par[9]
    kSyE2f2 = par[10]
    kDeE2f1 = par[11]
    kDeE2f2 = par[12]
    kPhE2f = par[13]
    kDpE2f1 = par[14]
    kDpE2f2 = par[15]
    kSyEmi1 = par[16]
    kSyEmi2 = par[17]
    kDeEmi1 = par[18]
    kDeEmi2 = par[19]
    kDeEmi3 = par[20]
    kPhEmiA = par[21]
    kPhEmiB = par[22]
    kDpEmi = par[23]
    kSyCa1 = par[24]
    kSyCa2 = par[25]
    kDeCa1 = par[26]
    kDeCa2 = par[27]
    kDeCa3 = par[28]
    kDpCdh = par[29]
    kPhCdhA = par[30]
    kPhCdhE = par[31]
    kPhCdhB = par[32]
    kDiACE = par[33]
    kAsACE = par[34]
    kDiACdh = par[35]
    kAsACdh = par[36]
    kPhCeA = par[37]
    kPhCeE = par[38]
    kPhFoxE = par[39]
    kPhFoxA = par[40]
    kPhFoxB = par[41]
    kDpFox = par[42]
    kSyCb1 = par[43]
    kSyCb2 = par[44]
    kDeCb1 = par[45]
    kDeCb2 = par[46]
    kDeCb3 = par[47]
    kPhEnsa = par[48]
    kDpEnsa = par[49]
    kPhGw = par[50]
    kDpGw1 = par[51]
    kDpGw2 = par[52]
    kWee1 = par[53]
    kWee2 = par[54]
    kPhWeeA = par[55]
    kPhWeeB = par[56]
    kDpWee = par[57]
    kCdc25_1 = par[58]
    kCdc25_2 = par[59]
    kPhC25A = par[60]
    kDpCdc25 = par[61]
    kDipEB55 = par[62]
    kAspEB55 = par[63]
    kSyCdc_1 = par[64]
    kSyCdc_2 = par[65]
    kDeCdc_1 = par[66]
    kDeCdc_2 = par[67]
    kDipACdc = par[68]
    kAspACdc = par[69]
    kPhApcA = par[70]
    kPhApcB = par[71]
    kDpApc = par[72]
    kSyFox1 = par[73]
    kSyFox2 = par[74]
    kDeFox1 = par[75]
    kDeFox2 = par[76]
    kAsEPx = par[77]
    kDiEPx = par[78]
    kAsFPcb = par[79]
    kDiFPcb = par[80]
    kAsFPcdc = par[81]
    kDiFPcdc = par[82]
    f1 = par[83]
    f2 = par[84]
    f3 = par[85]
    f4 = par[86]
    kPhC25B = par[87]
    
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

  def simulate(self, t=linspace(0, 10, 11)):
    self._result = odeint(self.dx_dt, self._x_0, t)
    self._t = t
    return (self._t, self._result, self._labels)

  @property
  def result(self):
    return (self._t, self._result, self._labels)

  @property
  def labels(self):
    return self._labels

  @property
  def params(self):
    return self._params

  @params.setter
  def params(self, value):
    if len(value) != len(self._params) or not isinstance(value, list):
      raise ValueError('params must be a list of length {} but is {} of length {}'.format(len(self._params), type(value), len(value)))
    self._params = value
    self.simulate()

  @property
  def x_0(self):
    return self._x_0

  @x_0.setter
  def x_0(self, value):
    if len(value) != len(self._x_0) or not isinstance(value, list):
      raise ValueError('x_0 must be a list of length {} but is {} of length {}'.format(len(self._x_0), type(value), len(value)))
    self._x_0 = value
    self.simulate()
