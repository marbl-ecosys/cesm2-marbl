import os
import numpy as np
import xarray as xr


os.environ['CESMDATAROOT'] = '/glade/scratch/mclong/inputdata'
import pop_tools

def calc_cfc12sol(S, T):
    return _calc_cfcsol(S, T, 12)

def calc_cfc11sol(S, T):
    return _calc_cfcsol(S, T, 11)


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


def mld_dsigma(SALT, TEMP, KMT, dsigma, z_t=None):
    dims_in = SALT.dims
    assert dims_in == TEMP.dims, 'dimension mismatch'
    if z_t is None:
        try:
            z_t = SALT.z_t if 'z_t' in SALT.coords else TEMP.z_t
        except:
            raise ValueError('z_t not found in SALT or TEMP coords; pass it in explicitly')
    
    lateral_dims = dims_in[-2:]
    vertical_dim = dims_in[-3]
    other_dims = set(dims_in) - set(lateral_dims) - {vertical_dim}    

    rho = pop_tools.eos(SALT, TEMP, depth=z_t * 1e-2)
    rho = rho - 1000.

    dims_out_ordered = [d for d in dims_in if d != vertical_dim]
    if not other_dims:
        rho = rho.expand_dims({'dummy': 1})
        other_dims = ('dummy',)
        dims_out_ordered = [d for d in dims_in if d != vertical_dim] + ['dummy']

    rho_stack = rho.stack(other_dims=other_dims, lateral_dims=lateral_dims,)

    kmt_stack = KMT.stack(lateral_dims=lateral_dims)
    mld_stack = xr.full_like(rho_stack.isel({vertical_dim: 0}), 
                             fill_value=np.nan).drop(vertical_dim)
    
    mld_stack.name = f'MLD_{dsigma:g}'.replace('0.', '')

    z_t_data = z_t.values * 1e-2 # convert from cm to m 
    for l in range(len(rho_stack.other_dims)):
        for ij in range(len(rho_stack.lateral_dims)):
            kmt_ij = kmt_stack[ij].values
            if kmt_ij == 0:
                continue
            rho_ijl = rho_stack[:kmt_ij, l, ij]
            z = z_t_data[:kmt_ij]
            rho_ijl_sort, I = np.unique(rho_ijl, return_index=True)
            z = z[I]
            if len(I) >= 3:
                mld_stack[l, ij] = np.interp(rho_ijl_sort[0] + dsigma, 
                                             rho_ijl_sort, z)
                
    mld = mld_stack.unstack().transpose(*dims_out_ordered)
    if 'dummy' in mld.dims:
        mld = mld.isel(dummy=0).drop('dummy')
    mld.attrs['long_name'] = 'MLD'
    mld.attrs['definition'] = f'$\Delta\sigma = {dsigma:0.3f}$'
    mld.attrs['units'] = 'm'
    return mld