def read_CESM_var(time_slice, variable, mean_dims=None, postprocess=None):
    import intake
    import intake_esm
    import xarray as xr
    
    is_derived = variable in derived_vars_defined

    if is_derived:
        dv = derived_var(variable)
        varlist = dv.dependent_vars
    else:
        varlist = [variable]
        
    catalog = intake.open_esm_datastore(
        'data/campaign-cesm2-cmip6-timeseries.json'
    )

    dq = (catalog
          .search(experiment='historical', 
                  component='ocn', 
                  variable=varlist
                 )
          .to_dataset_dict(
              cdf_kwargs={'chunks': {'time': 4}},
          ))

    keep_vars = ['REGION_MASK', 'z_t', 'dz', 'TAREA', 'TLONG', 'KMT', 
                 'TLAT', 'time', 'time_bound', 'member_id', 
                 'ctrl_member_id',] + varlist

    dset = dq['ocn.historical.pop.h']
    dset = dset[[v for v in keep_vars if v in dset]].sel(time=time_slice)

    if is_derived:
        dset = dv.compute(dset)
        if variable == 'Cant_v1':
            dset = dset.rename({'Cant': 'Cant_v1'})
            
    if mean_dims is not None:
        with xr.set_options(keep_attrs=True):
            dset = dset.mean(dim=mean_dims)
    
    if postprocess is not None:
        dset = postprocess(dset)
        
    for v in dset.variables:
        if '_FillValue' not in dset[v].encoding:
            dset[v].encoding['_FillValue'] = None
        
    return dset.compute()


def read_obs(src, variable=None, freq='monthly'):
    import os
    import xarray as xr

    if src not in ['WOA', 'SeaWiFS']:
        raise ValueError(f'{src} is not a valid source')

    if freq in ['month', 'monthly']:
        freq = 'mon'
    if freq == 'annual':
        freq = 'ann'
    
    xr_kwargs = dict()
    if src == 'WOA':
        if freq not in ['mon']:
            raise ValueError(f'{freq} is not a valid frequency of data for {src}')
        # varmap is used to define the filename
        varmap = {'NO3' : 'n', 'PO4' : 'p', 'SiO3' : 'i'}
        if variable not in varmap:
            raise ValueError(f'{variable} is not a valid variable for WOA')
        root_dir=os.path.join(os.path.sep, 'glade', 'p', 'cgd', 'oce', 'projects',
                              'cesm2-marbl', 'woa2018-data', 'POP_gx1v7', 'annual')
        filename=f'woa18_all_{varmap[variable]}00_gx1v7.nc'
        xr_kwargs['decode_times'] = False

    if src == 'SeaWiFS':
        if freq not in ['mon', 'ann', 'JJA', 'DJF']:
            raise ValueError(f'{freq} is not a valid frequency of data for {src}')
        root_dir=os.path.join(os.path.sep, 'glade', 'p', 'cgd', 'oce', 'projects',
                              'cesm2-marbl', 'seaWIFS-data')
        filename=f'seaWIFS.chl_gsm.{freq}_climo.Sep1997_Dec2010.nc'

    # Read in dataset, use first time dim from WOA
    ds=xr.open_dataset(os.path.join(root_dir, filename), **xr_kwargs)
    if src == 'WOA':
        ds = ds.isel(time=0)
    return(ds)


def get_fesedflux_forcing():
    import xarray as xr
    µmolm2d_to_mmolm2yr = 1e-3 * 365.
    
    file_fesedflux = '/glade/p/cesmdata/cseg/inputdata/ocn/pop/gx1v6/forcing/fesedfluxTot_gx1v6_cesm2_2018_c180618.nc'
    file_feventflux = '/glade/p/cesmdata/cseg/inputdata/ocn/pop/gx1v6/forcing/feventflux_gx1v6_5gmol_cesm1_97_2017.nc'

    dsi = xr.merge((
        xr.open_dataset(file_feventflux).rename({'FESEDFLUXIN': 'Fe_ventflux'}) * µmolm2d_to_mmolm2yr,
        xr.open_dataset(file_fesedflux).rename({'FESEDFLUXIN': 'Fe_sedflux'}) * µmolm2d_to_mmolm2yr,
    )).rename({'z': 'z_t', 'y': 'nlat', 'x': 'nlon'})

    dsi['Fe_ventflux'] = dsi.Fe_ventflux.sum('z_t') # since units are already /m^2, we can just sum to integrate in z
    dsi['Fe_sedflux'] = dsi.Fe_sedflux.sum('z_t')

    dsi.Fe_ventflux.attrs['units'] = 'mmol m$^{-2}$ yr$^{-1}$'
    dsi.Fe_sedflux.attrs['units'] = 'mmol m$^{-2}$ yr$^{-1}$'
    return dsi

def plot_surface_vals(variable, ds, da, da_obs, obs_src='obs', levels=None, bias_levels=None, force_units=None):
    import matplotlib.pyplot as plt
    from matplotlib.colors import BoundaryNorm
    from matplotlib.ticker import MaxNLocator

    import cartopy
    import cartopy.crs as ccrs

    if type(ds) == dict:
        ds = ds[variable]
    if type(da) == dict:
        da = da[variable]
    if type(da_obs) == dict:
        da_obs = da_obs[variable]
    TLONG = ds.TLONG
    TLAT = ds.TLAT
    computed = da
    observed = da_obs
    bias=computed-observed
    # If user provided specified units, assume these are pint objects
    # Otherwise treat them as dataarrays
    if force_units:
        computed = computed.to(force_units).magnitude
        observed = observed.to(force_units).magnitude
        bias = bias.to(force_units).magnitude
    else:
        computed = computed.data
        observed = observed.data
        bias = bias.data

    # Determine contours
    if levels is None:
#         pass
# #         levels = MaxNLocator(nbins=len(levels)-1).tick_values(levels)
#     else:
        min_lev = 0
        if variable == 'NO3':
            max_lev = 42
        if variable == 'PO4':
            max_lev = 3.2
        if variable == 'SiO3':
            max_lev = 180
        if variable == 'totChl':
            max_lev = 20
        levels = MaxNLocator(nbins=15).tick_values(min_lev, max_lev)
    cmap = plt.get_cmap('rainbow')
    norm = BoundaryNorm(levels, ncolors=cmap.N)
    bias_kwargs=dict()
    if bias_levels is not None:
        bias_kwargs['norm'] = BoundaryNorm(bias_levels, ncolors=cmap.N)


    fig = plt.figure(figsize=(8, 14))
    ax = plt.subplot(3, 1, 1, projection=ccrs.Robinson(central_longitude=305.0))

    pc = ax.pcolormesh(TLONG,
                       TLAT,
                       computed,
                       cmap=cmap,
                       transform=ccrs.PlateCarree(),
                       norm=norm)

    # ax.add_feature(cartopy.feature.NaturalEarthFeature('physical','land','110m',
    #                                                    edgecolor='face',
    #                                                    facecolor='lightgray'))
    ax.set_global()
    ax.coastlines(linewidth=0.5)
    cb = plt.colorbar(pc, shrink=0.6)
    ax.set_title(f'{variable} from POP run');
    if force_units:
        cb.set_label(force_units)
    else:
        cb.set_label(da.attrs['units'])

    ax = plt.subplot(3, 1, 2, projection=ccrs.Robinson(central_longitude=305.0))

    pc = ax.pcolormesh(TLONG,
                       TLAT,
                       observed,
                       cmap=cmap,
                       transform=ccrs.PlateCarree(),
                       norm=norm)

    # ax.add_feature(cartopy.feature.NaturalEarthFeature('physical','land','110m',
    #                                                    edgecolor='face',
    #                                                    facecolor='lightgray'))
    ax.set_global()
    ax.coastlines(linewidth=0.5)
    cb = plt.colorbar(pc, shrink=0.6)
    ax.set_title(f'{variable} from {obs_src}');
    if force_units:
        cb.set_label(force_units)
    else:
        cb.set_label(da_obs.attrs['units'])


    ax = plt.subplot(3, 1, 3, projection=ccrs.Robinson(central_longitude=305.0))

    pc = ax.pcolormesh(TLONG,
                       TLAT,
                       bias,
                       cmap=plt.get_cmap('bwr'),
                       transform=ccrs.PlateCarree(), **bias_kwargs)

    # ax.add_feature(cartopy.feature.NaturalEarthFeature('physical','land','110m',
    #                                                    edgecolor='face',
    #                                                    facecolor='lightgray'))
    ax.set_global()
    ax.coastlines(linewidth=0.5)
    cb = plt.colorbar(pc, shrink=0.6)
    ax.set_title(f'Bias in {variable}');
    if force_units:
        cb.set_label(force_units)
    else:
        cb.set_label(da.attrs['units'])


def plot_global_profile(variables, units, ds, da, obs):
    import matplotlib.pyplot as plt

    plt_cnt = len(variables)
    z = ds[variables[0]].z_t.data * units[ds[variables[0]]['z_t'].attrs['units']]

    fig = plt.figure(figsize=(4*plt_cnt, 4))
    for n, variable in enumerate(variables):
        computed = da[variable]
        observed = obs[variable]

        ax = fig.add_subplot(1, plt_cnt, n+1)

        ax.plot(computed.magnitude, (z.to('m')).magnitude, 'b-', label='CESM2', linewidth=2)
        ax.plot(observed.magnitude, (z.to('m')).magnitude, 'r:', label='WOA2018', linewidth=2)

        ax.set_title(f'Global Profile of {variable}')
        ax.set(xlabel='concentration (mmol m$^{-3}$)')
        if n == 0:
            ax.set(ylabel='depth (m)')
        else:
            ax.set_yticklabels('')
        plt.gca().invert_yaxis()
        ax.legend()


def return_magnitude_in_units(pint_obj, units):
    return((pint_obj.to('mmol/m^3')).magnitude)


def plot_zonal_averages_by_region(variables, region, da, obs, lat, z):
    import matplotlib.pyplot as plt
    import numpy as np

    if type(variables) != list:
        variables = [variables]

    # For now use levels defined in Kristen's notebook
    levels = dict()
    levels['NO3'] = np.arange(0, 48, 4) # 0, 4, 8, ..., 40, 44
    levels['PO4'] = np.arange(0, 4.4, 0.4) # 0, 0.4, 0.8, ..., 3.6, 4.0
    levels['SiO3'] = np.concatenate((np.arange(0, 100, 10), np.arange(100, 300, 20))) # 0, 10, ... 90, 100, 120, 140, ..., 260, 280

    bias_levels = dict()
    bias_levels['NO3'] = np.arange(-18, 20, 3) # -18, -15, ..., 15, 18
    bias_levels['PO4'] = np.arange(-0.9, 1, 0.1) # -0.9, -0.8, ..., 0.8, 0.9
    bias_levels['SiO3'] = np.arange(-80, 85, 5) # -80, -75, -70, ..., 75, 80

    fig = plt.figure(figsize=(8*len(variables),12))
    plt.suptitle(f'Zonal Means ({region})', fontsize=14)

    # TOP: CESM

    for n, variable in enumerate(variables):
        cesm_out = return_magnitude_in_units(da[region][variable], 'mmol/m^3')
        woa_vals = return_magnitude_in_units(obs[region][variable], 'mmol/m^3')
        bias = cesm_out - woa_vals

        ax = fig.add_subplot(3, len(variables), n + 1)
        ax.set_title(f' {variable}\nCESM2')
        pc=ax.contourf(lat, z, cesm_out, levels=levels[variable], cmap='rainbow', extend='both')
        pc2 = ax.contour(lat, z, cesm_out, levels[variable], colors='k')
        ax.clabel(pc2, colors = 'k', fmt = '%2.1f', fontsize=10)
        if n==0:
            ax.set(ylabel='depth (m)')
        else:
            ax.set(ylabel='', yticklabels='')
        ax.set(xlabel='', xticklabels='')
        ax.invert_yaxis()

        #MIDDLE: WOA

        ax = fig.add_subplot(3, len(variables), (n + 1) + len(variables))
        ax.set_title(f'WOA')
        pc=ax.contourf(lat, z, woa_vals, levels=levels[variable], cmap='rainbow', extend='both')
        pc2 = ax.contour(lat, z, woa_vals, levels[variable], colors='k')
        ax.clabel(pc2, colors = 'k', fmt = '%2.1f', fontsize=10)
        if n==0:
            ax.set(ylabel='depth (m)')
        else:
            ax.set(ylabel='', yticklabels='')
        ax.set(xlabel='', xticklabels='')
        ax.invert_yaxis()

        #BOTTOM: Bias

        ax = fig.add_subplot(3, len(variables), (n + 1) + 2*len(variables))
        ax.set_title(f'CESM2 - WOA')
        pc=ax.contourf(lat, z, bias, levels=bias_levels[variable], cmap='bwr', extend='both')
        pc2 = ax.contour(lat, z, bias, bias_levels[variable], colors='k')
        ax.clabel(pc2, colors = 'k', fmt = '%2.1f', fontsize=10)
        if n==0:
            ax.set(ylabel='depth (m)')
        else:
            ax.set(ylabel='', yticklabels='')
        ax.set(xlabel='Latitude')
        ax.invert_yaxis()

    
def _ensure_variables(ds, req_var):
    """ensure that required variables are present"""
    missing_var_error = False
    for v in req_var:
        if v not in ds:
            print('ERROR: Missing required variable: {v}')
            missing_var_error = True
    if missing_var_error:
        raise ValueError('Variables missing')

                
def derive_var_pCFC11(ds, drop_derivedfrom_vars=True):
    """compute pCFC11"""
    from calc import calc_cfc11sol

    ds['pCFC11'] = ds['CFC11'] * 1e-9 / calc_cfc11sol(ds.SALT, ds.TEMP)
    ds.pCFC11.attrs['long_name'] = 'pCFC-11'
    ds.pCFC11.attrs['units'] = 'patm'
    if 'coordinates' in ds.TEMP.attrs:
        ds.pCFC11.attrs['coordinates'] = ds.TEMP.attrs['coordinates']
    ds.pCFC11.encoding = ds.TEMP.encoding

    if drop_derivedfrom_vars:
        ds = ds.drop(['CFC11', 'TEMP', 'SALT'])

    return ds


def derive_var_pCFC12(ds, drop_derivedfrom_vars=True):
    """compute pCFC12"""
    from calc import calc_cfc12sol

    ds['pCFC12'] = ds['CFC12'] * 1e-9 / calc_cfc12sol(ds['SALT'],ds['TEMP'])
    ds.pCFC12.attrs['long_name'] = 'pCFC-12'
    ds.pCFC12.attrs['units'] = 'patm'
    if 'coordinates' in ds.TEMP.attrs:
        ds.pCFC12.attrs['coordinates'] = ds.TEMP.attrs['coordinates']
    ds.pCFC12.encoding = ds.TEMP.encoding

    if drop_derivedfrom_vars:
        ds = ds.drop(['CFC12', 'TEMP', 'SALT'])

    return ds        

def derive_var_Cant(ds, drop_derivedfrom_vars=True):
    """compute Cant"""

    ds['Cant'] = ds['DIC'] - ds['DIC_ALT_CO2']
    ds.Cant.attrs = ds.DIC.attrs
    ds.Cant.attrs['long_name'] = 'Anthropogenic CO$_2$'

    if 'coordinates' in ds.DIC.attrs:
        ds.Cant.attrs['coordinates'] = ds.DIC.attrs['coordinates']
    ds.Cant.encoding = ds.DIC.encoding

    if drop_derivedfrom_vars:
        ds = ds.drop(['DIC', 'DIC_ALT_CO2'])

    return ds     

def derive_var_Del14C(ds, drop_derivedfrom_vars=True):
    """compute Del14C"""

    ds['Del14C'] = 1000. * (ds['ABIO_DIC14'] / ds['ABIO_DIC'] - 1.)
    ds.Del14C.attrs = ds.ABIO_DIC14.attrs
    ds.Del14C.attrs['long_name'] = '$\Delta^{14}$C'
    ds.Del14C.attrs['units'] = 'permille'    

    if 'coordinates' in ds.ABIO_DIC14.attrs:
        ds.Del14C.attrs['coordinates'] = ds.ABIO_DIC14.attrs['coordinates']
    ds.Del14C.encoding = ds.ABIO_DIC14.encoding

    if drop_derivedfrom_vars:
        ds = ds.drop(['ABIO_DIC14', 'ABIO_DIC'])

    return ds     


derived_vars_defined = dict(
    pCFC11=dict(
        dependent_vars=['CFC11', 'TEMP', 'SALT'],
        method=derive_var_pCFC11,
    ),
    pCFC12=dict(
        dependent_vars=['CFC12', 'TEMP', 'SALT'],
        method=derive_var_pCFC12,
    ),
    Cant=dict(
        dependent_vars=['DIC', 'DIC_ALT_CO2'],
        method=derive_var_Cant,        
    ),
    Cant_v1=dict(
        dependent_vars=['DIC', 'DIC_ALT_CO2'],
        method=derive_var_Cant,        
    ),    
    Del14C=dict(
        dependent_vars=['ABIO_DIC14', 'ABIO_DIC'],
        method=derive_var_Del14C,        
    ),        
)
class derived_var(object):      
    def __init__(self, varname):
        assert varname in derived_vars_defined
        self.varname = varname
        self.dependent_vars = derived_vars_defined[self.varname]['dependent_vars']
        self._callable = derived_vars_defined[self.varname]['method']
        
    def compute(self, ds, **kwargs):
        _ensure_variables(ds, self.dependent_vars)
        return self._callable(ds, **kwargs)
    
    
    