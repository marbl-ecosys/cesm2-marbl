def read_CESM_var(time_slice, variable, mean_dims=None):
    import intake
    import intake_esm

    # Define catalog
    catalog = intake.open_esm_datastore('data/campaign-cesm2-cmip6-timeseries.json')
    dq = catalog.search(experiment='historical', component='ocn', variable=variable).to_dataset_dict(cdf_kwargs={'chunks':{'time': 4}})

    mean_kwargs = dict()
    if mean_dims:
        mean_kwargs['dim'] = mean_dims

    # Define datasets
    dataset = dq['ocn.historical.pop.h']

    keep_vars = ['REGION_MASK', 'z_t', 'dz', 'TAREA', 'TLONG', 'TLAT', 'time', 'time_bound', 'member_id', 'ctrl_member_id', variable]
    dataset = dataset.drop([v for v in dataset.variables if v not in keep_vars]).sel(time=time_slice).mean(**mean_kwargs).compute()

    return(dataset)

##################################################

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

##################################################

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


##################################################

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

##################################################

def return_magnitude_in_units(pint_obj, units):
    return((pint_obj.to('mmol/m^3')).magnitude)

##################################################

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

        
        