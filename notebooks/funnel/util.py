import numpy as np
import xarray as xr
import cftime

from . defaults import unit_str_replacements


def compute_time(ds, time_coord_name='time'):
    """reset time in a dataset to be the mean of the time_bounds variable """
    
    ds = ds.copy(deep=True)

    time_bnd_name = get_time_bnd_name(ds, time_coord_name)            
    time_bnds_dim2 = ds[time_bnd_name].dims[1]    

    time_attrs = ds[time_coord_name].attrs
    time_encoding = ds[time_coord_name].encoding
    time_encoding.update(
        dict(
            units=time_attrs.pop('units'),
            calendar=time_attrs.pop('calendar'),
            _FillValue=None,
        )
    )

    ds[time_coord_name] = xr.DataArray(
        cftime.num2date(
            ds[time_bnd_name].mean(dim=time_bnds_dim2), 
            units=time_encoding['units'], 
            calendar=time_encoding['calendar'],
        ), 
        dims=('time'),
        attrs=time_attrs,
    )
    ds[time_coord_name].encoding = time_encoding
    
    return ds


def fix_units(ds):
    """fix the units"""
    
    ds = ds.copy(deep=True)
    
    for v in ds.variables:
        if 'units' in ds[v].attrs:
            old_units = ds[v].attrs['units']            
            if old_units in unit_str_replacements:
                ds[v].attrs['units'] = unit_str_replacements[old_units]

    return ds


def calc_ann_mean(ds, time_coord_name='time'):
    """compute annual means of a dataset"""
    
    time_bnd_name = get_time_bnd_name(ds, time_coord_name)            
    tb_dim2 = ds[time_bnd_name].dims[1]

    group_by_year = time_coord_name+'.year'
    rename = {'year': time_coord_name}
    
    # compute time weights as time-bound diff
    time_wgt = ds[time_bnd_name].diff(dim=tb_dim2)
    time_wgt_grouped = time_wgt.groupby(group_by_year)
    time_wgt = time_wgt_grouped / time_wgt_grouped.sum(dim=xr.ALL_DIMS)

    nyr = len(time_wgt_grouped.groups)
    time_wgt = time_wgt.squeeze()

    # ensure that weights sum to 1
    np.testing.assert_almost_equal(time_wgt.groupby(group_by_year).sum(dim=xr.ALL_DIMS), 
                                   np.ones(nyr))

    # set non-time related vars to coords to avoid xarray adding a time-dim
    nontime_vars = set([v for v in ds.variables if 'time' not in ds[v].dims]) - set(ds.coords)
    dsop = ds.set_coords(nontime_vars)

    # compute the annual means
    with xr.set_options(keep_attrs=True):
        ds_ann = (dsop * time_wgt).groupby(group_by_year).sum(dim='time')

    # construct outgoing time bounds
    time_bnd_min = ds[time_bnd_name].groupby(group_by_year).min()[:, 0]
    time_bnd_max = ds[time_bnd_name].groupby(group_by_year).max()[:, 1]    
    ds_ann[time_bnd_name] = xr.concat(
        (time_bnd_min, time_bnd_max), dim=tb_dim2
    ).transpose('year', tb_dim2)
    
    # rename time and put back the coords variable
    return ds_ann.reset_coords(nontime_vars).rename(rename)


def calc_time_mean(ds, time_coord_name='time'):

    time_bnd_name = get_time_bnd_name(ds, time_coord_name)            
        
    # compute time weights as time-bound diff
    time_wgt = ds[time_bnd_name].diff(dim=ds[time_bnd_name].dims[1]).squeeze()
    time_wgt = time_wgt / time_wgt.sum(dim=xr.ALL_DIMS)

    # set non-time related vars to coords to avoid xarray adding a time-dim
    nontime_vars = set([v for v in ds.variables if 'time' not in ds[v].dims]) - set(ds.coords)
    dsop = ds.set_coords(nontime_vars).drop(time_bnd_name)

    # compute the annual means
    with xr.set_options(keep_attrs=True):
        ds_mean = (dsop * time_wgt).sum(dim='time')

    # put back the coords variable
    return ds_mean.reset_coords(nontime_vars)


def get_time_bnd_name(ds, time_coord_name):
    """return the name of the time bounds variable"""
    try:
        return ds[time_coord_name].attrs['bounds']
    except:      
        print('-'*80)
        print('problem with the following dataset:')
        ds.info()
        raise ValueError('could not determine time bounds variable name')
        
        
