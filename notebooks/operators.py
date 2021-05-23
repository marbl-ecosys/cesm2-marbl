import numpy as np
import xarray as xr

import funnel as fn

import variable_defs

nmols_to_PgCyr = 1e-9 * 86400. * 365. * 12e-15


def resample_ann(ds):
    """compute the annual mean of an xarray.Dataset"""
    
    ds = ds.copy()
    
    # compute temporal weights using time_bound attr
    assert 'bounds' in ds.time.attrs, 'missing "bounds" attr on time'
    tb_name = ds.time.attrs['bounds']        
    assert tb_name in ds, f'missing "{tb_name}"'
    
    dim = ds[tb_name].dims[-1]
    ds['time'] = ds[tb_name].compute().mean(dim).squeeze()   
    
    # compute weigths from diff of time_bound
    weights = ds[tb_name].compute().diff(dim).squeeze()
    weights = weights.groupby('time.year') / weights.groupby('time.year').sum()
   
    # ensure they all add to one
    # TODO: build support for situations when they don't, i.e. define min coverage threshold
    nyr = len(weights.groupby('time.year'))
    np.testing.assert_allclose(weights.groupby('time.year').sum().values, np.ones(nyr))
        
    # ascertain which variables have time and which don't
    time_vars = [v for v in ds.data_vars if 'time' in ds[v].dims and v != tb_name]
    other_vars = list(set(ds.variables) - set(time_vars) - {tb_name, 'time'} )

    # compute
    with xr.set_options(keep_attrs=True):        
        return xr.merge((
            ds[other_vars],         
            (ds[time_vars] * weights).groupby('time.year').sum(dim='time'),
        )).rename({'year': 'time'})    


def global_mean(ds, normalize=True):
    """
    Compute the global mean on a POP dataset. 
    Return computed quantity in conventional units.
    """
    masked_area = ds.TAREA.where(ds.REGION_MASK > 0).fillna(0.)
    compute_vars = [
        v for v in ds 
        if 'time' in ds[v].dims and ('nlat', 'nlon') == ds[v].dims[-2:]
    ]
    other_vars = list(set(ds.variables) - set(compute_vars))

    with xr.set_options(keep_attrs=True):
        
        dso = (ds[compute_vars] * masked_area).sum(['nlat', 'nlon'])    
        if normalize:
            dso = dso[compute_vars] / masked_area.sum(['nlat', 'nlon'])
        else:
            for v in compute_vars:
                if v in variable_defs.C_flux_vars:
                    dso[v] = dso[v] * nmols_to_PgCyr
                    dso[v].attrs['units'] = 'Pg C yr$^{-1}$'
                
        return xr.merge([dso, ds[other_vars]])
    
