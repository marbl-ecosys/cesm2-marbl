import pop_tools

from . import util


def get_grid(model):
    """returnt the model grid"""
    if model in ['cesm1']:
        grid = pop_tools.get_grid('POP_gx1v6')
    elif model in ['cesm2']:
        grid = pop_tools.get_grid('POP_gx1v7')
    else:
        raise ValueError(f'unknown model: {model}') 
        
    grid['dz_150m'] = grid.dz.isel(z_t=slice(0, 15)).rename(z_t='z_t_150m')

    return util.fix_units(grid)


def compute_chl_surf(ds, drop=True):
    """compute surface chl"""

    ds = ds.copy(deep=True)
    ds['Chl_surf'] = (ds.diatChl + ds.spChl + ds.diazChl).isel(z_t_150m=0)
    ds.Chl_surf.attrs = ds.diatChl.attrs
    ds.Chl_surf.attrs['long_name'] = 'Surface chlorophyll'

    if drop:
        ds = ds.drop(['diatChl', 'spChl', 'diazChl'])

    return ds