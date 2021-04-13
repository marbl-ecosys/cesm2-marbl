import numpy as np
from numba import jit

def calc_cfc12sol(S, T):
    return _calc_cfcsol(S, T, 12)

def calc_cfc11sol(S,T):
    return _calc_cfcsol(S, T, 11)


#@jit(nopython=True)
def _calc_cfcsol(PS, PT, kn):  
    """
    Compute CFC11 and CFC12 Solubilities in seawater.
        Reference: Warner & Weiss (1985) , Deep Sea Research, vol32
    
    Translated from cfc11_mod.F90 (MCL, 2011)
    INPUT:
    PT: temperature (degree Celsius)
    PS: salinity
    kn: 11 = CFC11, 12 = CFC12
    OUTPUT:
    SOLUBILITY_CFC: returned value in mol/m3/pptv
    1 pptv = 1 part per trillion = 10^-12 atm = 1 picoatm
    """
    T0_Kelvin = 273.16 # this is indeed what POP uses
    c1000 = 1000.0
   
    assert kn in [11, 12], 'kn must be either 11 or 12'
    if kn == 11:
        a1 = -229.9261
        a2 = 319.6552
        a3 = 119.4471
        a4 = -1.39165
        b1 = -0.142382
        b2 = 0.091459
        b3 = -0.0157274
    elif kn == 12:
        a1 = -218.0971
        a2 = 298.9702
        a3 = 113.8049
        a4 = -1.39165
        b1 = -0.143566
        b2 = 0.091015
        b3 = -0.0153924
    
    WORK = ((PT + T0_Kelvin) * 0.01)

    #  coefficient for solubility in  mol/l/atm
    SOLUBILITY_CFC = np.exp( a1 + a2 / WORK + a3 * np.log ( WORK )
                           + a4 * WORK * WORK
                           + PS * ( ( b3 * WORK + b2 ) * WORK + b1 ) )

    #  conversion from mol/(l * atm) to mol/(m^3 * atm) to mol/(m3 * pptv)
    SOLUBILITY_CFC = c1000 * SOLUBILITY_CFC
    SOLUBILITY_CFC = 1.0e-12 * SOLUBILITY_CFC

    return SOLUBILITY_CFC # mol/m^3/patm

#-------------------------------------------------------------------------------
#--- main
#-------------------------------------------------------------------------------
if __name__ == '__main__':
    print(calc_cfc12sol(35., 15.))
    print(calc_cfc11sol(35., 15.))
