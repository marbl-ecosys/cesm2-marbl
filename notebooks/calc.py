import os

from scipy import interpolate
import numpy as np
import xarray as xr

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




def mld_dsigma(SALT, TEMP, dsigma=0.03, rho_chunks={'nlat': 16, 'nlon': 16}):
    """
    Compute MLD based on ∆σ criterion. Uses xarray.map_blocks.
    
    Parameters
    ----------
    
    SALT : xarray.DataArray
      Salinity
    TEMP : xarray.DataArray
      Potential temperature
    dsigma : float, optional
      The value for ∆σ.
      
    Returns
    -------
    
    MLD : xarray.DataArray
      The MLD (m) defined as the point in the water column where
      density exceeds rho[0] + dsigma.      
    """
    
    # determine dimensionality
    dims_in = SALT.dims
    assert dims_in == TEMP.dims, 'dimension mismatch'
    assert 'z_t' in SALT.coords, 'z_t not found in SALT coords'

   
    # drop ancillary coordinates (this may not be necessary)
    SALT = SALT.reset_coords(drop=True)
    TEMP = TEMP.reset_coords(drop=True)
    
    # compute density
    rho = pop_tools.eos(SALT.chunk({'z_t': 10}), 
                        TEMP.chunk({'z_t': 10}), 
                        depth=SALT.z_t * 0.).compute()
    
    if 'nlat' in rho.dims:
        rho = rho.assign_coords({
            'nlat': xr.DataArray(np.arange(len(SALT.nlat)), dims=('nlat')),
            'nlon': xr.DataArray(np.arange(len(SALT.nlon)), dims=('nlon')),
        })
    rho = rho.chunk(rho_chunks).persist()
    
    # compute and return MLD
    template = rho.isel(z_t=0).drop('z_t')
    template.attrs['long_name'] = 'MLD'
    template.attrs['units'] = SALT.z_t.attrs['units']
    template.name = 'MLD'
    
    return xr.map_blocks(
        _interp_mld, rho,
        kwargs=dict(dsigma=dsigma), 
        template=template,
    )
    
    
def _interp_mld(rho_in, dsigma=0.03):
    """compute MLD at point using interpolation"""
    
    non_vertical_dims = [d for d in rho_in.dims if d not in ['z_t']]
    rho_stack = rho_in.stack(non_vertical_dims=non_vertical_dims) 
    mld_stack = xr.full_like(rho_stack.isel(z_t=slice(0, 1)), fill_value=np.nan)
    mld_stack.name = 'MLD'
    z_t = rho_in.z_t
    
    for i in range(len(rho_stack.non_vertical_dims)):
        # TODO: get more specific about extrapolation rules
        #       desired behavior: all rho < rho[0] + dsigma, 
        #       return z_t[np.argmax(rho)]?
        if np.isnan(rho_stack[:, i]).all():
            continue

        if (rho_stack[:, i] < rho_stack[0, i] + dsigma).all():
            k = np.where(~np.isnan(rho_stack[:, i]))[0]
            mld_stack[:, i] = z_t[k[-1]]
        else: 
            f = interpolate.interp1d(
                rho_stack[:, i], z_t, 
                assume_sorted=False,
            )
            mld_stack[:, i] = f(rho_stack[0, i] + dsigma)

    return mld_stack.unstack().isel(z_t=0, drop=True).transpose(*non_vertical_dims)



# def mld_dsigma(SALT, TEMP, dsigma=0.03):
#     """
#     Compute MLD based on ∆σ criterion. Uses xarray.map_blocks.
    
#     Parameters
#     ----------
    
#     SALT : xarray.DataArray
#       Salinity
#     TEMP : xarray.DataArray
#       Potential temperature
#     dsigma : float, optional
#       The value for ∆σ.
      
#     Returns
#     -------
    
#     MLD : xarray.DataArray
#       The MLD (m) defined as the point in the water column where
#       density exceeds rho[0] + dsigma.      
#     """
    
#     # determine dimensionality
#     dims_in = SALT.dims
#     assert dims_in == TEMP.dims, 'dimension mismatch'
#     assert 'z_t' in SALT.coords, 'z_t not found in SALT coords'

#     # assume vertical dimension is called "z_t"
#     non_vertical_dims = set(dims_in) - {'z_t'}
    
#     # drop ancillary coordinates (this may not be necessary)
#     SALT = SALT.reset_coords(drop=True)
#     TEMP = TEMP.reset_coords(drop=True)
    
#     # define chunks in each non-vertical dimension
#     chunk_dict = {k: 1 for k in non_vertical_dims}
    
#     # compute density
#     rho = pop_tools.eos(SALT, 
#                         TEMP, 
#                         depth=SALT.z_t * 1e-2).compute()
        
#     print('\trho computation complete')
#     # compute and return MLD
#     rho = rho.chunk(chunk_dict)
#     return xr.map_blocks(
#         _interp_mld, rho,
#         kwargs=dict(dsigma=dsigma), 
#         template=rho.isel(z_t=0).drop('z_t').reset_coords(drop=True),
#     ).to_dataset()
    
    
# def _interp_mld(rho_1d, dsigma=0.03):
#     """compute MLD at point using interpolation"""
#     mld = np.empty(rho_1d.isel(z_t=0).shape)
#     if np.isnan(rho_1d).all():
#         mld[:] = np.nan
#     else:
#         # TODO: get more specific about extrapolation rules
#         #       desired behavior: all rho < rho[0] + dsigma, 
#         #       return z_t[np.argmax(rho)]?
#         f = interpolate.interp1d(
#             rho_1d.squeeze(), rho_1d.z_t*1e-2, 
#             assume_sorted=False,
#         )
#         mld[:] = f(rho_1d.isel(z_t=0) + dsigma)
    
#     return xr.DataArray(
#         mld, 
#         dims=[k for k in rho_1d.dims if k != 'z_t'], 
#         coords={k: v for k, v in rho_1d.coords.items() if k != 'z_t'},
#         attrs={'long_name': 'MLD', 'units': 'm', 'note': f'Dsigma = {dsigma:g}'}
#     )    