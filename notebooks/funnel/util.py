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
        
        
def test_dataset():
    """Generate a simple test dataset"""
    
    start_date = np.array([0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334], dtype=np.float64)
    start_date = np.append(start_date, start_date + 365)
    end_date = np.array([31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365], dtype=np.float64)
    end_date = np.append(end_date, end_date + 365)

    ds = xr.Dataset(coords={'time': 24, 'lat': 2, 'lon': 2, 'd2': 2})
    ds['time'] = xr.DataArray(end_date, dims='time')
    ds['lat'] = xr.DataArray([0, 1], dims='lat')
    ds['lon'] = xr.DataArray([0, 1], dims='lon')
    ds['d2'] = xr.DataArray([0, 1], dims='d2')
    ds['time_bound'] = xr.DataArray(
        np.array([start_date, end_date]).transpose(), dims=['time', 'd2']
    )
    ds['variable_1'] = xr.DataArray(
        np.append(
            np.zeros([12, 2, 2], dtype='float32'), np.ones([12, 2, 2], dtype='float32'), axis=0
        ),
        dims=['time', 'lat', 'lon'],
    )
    ds.variable_1.attrs['description'] = 'All zeroes for year 1, all ones for year 2'
    
    ds['variable_2'] = xr.DataArray(
        np.append(
            np.ones([12, 2, 2], dtype='float32'), np.zeros([12, 2, 2], dtype='float32'), axis=0
        ),
        dims=['time', 'lat', 'lon'],
    )
    ds.variable_2.attrs['description'] = 'All ones for year 1, all zeroes for year 2'
    
    ds['non_time_variable_1'] = xr.DataArray(np.ones((2, 2)), dims=['lat', 'lon'])
    
    ds.time.attrs['units'] = 'days since 0000-01-01 00:00:00'
    ds.time.attrs['calendar'] = 'noleap'
    ds.time.attrs['bounds'] = 'time_bound'

    return ds        