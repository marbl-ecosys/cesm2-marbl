import os
import subprocess
import tempfile
import xarray as xr
import numpy as np

def zonal_mean_via_fortran(ds, var, grid=None, rmask_file=None):
    """
    Write ds to a temporary netCDF file, compute zonal mean for
    a given variable based on Keith L's fortran program, read
    resulting netcdf file, and return the new xarray dataset

    If three_ocean_regions=True, use a region mask that extends the
    Pacific, Indian, and Atlantic to the coast of Antarctica (and does
    not provide separate Arctic Ocean, Lab Sea, etc regions)
    """

    # xarray doesn't require the ".nc" suffix, but it's useful to know what the file is for
    ds_in_file = tempfile.NamedTemporaryFile(suffix='.nc')
    ds_out_file = tempfile.NamedTemporaryFile(suffix='.nc')
    ds.to_netcdf(ds_in_file.name)

    # Set up location of the zonal average executable
    za_exe = os.path.join(os.path.sep,
                          'glade',
                          'u',
                          'home',
                          'klindsay',
                          'bin',
                          'zon_avg',
                          'za')
    if grid:
        # Point to a file that contains all necessary grid variables
        # I think that is KMT, TLAT, TLONG, TAREA, ULAT, ULONG, UAREA, and REGION_MASK
        if grid == 'gx1v7':
            grid_file = os.path.join(os.path.sep,
                                     'glade',
                                     'collections',
                                     'cdg',
                                     'timeseries-cmip6',
                                     'b.e21.B1850.f09_g17.CMIP6-piControl.001',
                                     'ocn',
                                     'proc',
                                     'tseries',
                                     'month_1',
                                     'b.e21.B1850.f09_g17.CMIP6-piControl.001.pop.h.NO3.000101-009912.nc')
        else:
            print(f'WARNING: no grid file for {grid}, using xarray dataset for grid vars')
            grid_file = ds_in_file.name
    else:
        # Assume xarray dataset contains all needed fields
        grid_file = ds_in_file.name
    if rmask_file:
        region_mask = ['-rmask_file', rmask_file]
    else:
        region_mask = []

    # Set up the call to za with correct options
    za_call = [za_exe,
               '-v', var] + \
              region_mask + \
              [
               '-grid_file', grid_file,
               '-kmt_file', grid_file,
               '-O', '-o', ds_out_file.name, # -O overwrites existing file, -o gives file name
               ds_in_file.name]

    # Use subprocess to call za, allows us to capture stdout and print it
    proc = subprocess.Popen(za_call, stdout=subprocess.PIPE)
    (out, err) = proc.communicate()
    if not out:
        # Read in the newly-generated file
        print('za ran successfully, writing netcdf output')
        ds_out = xr.open_dataset(ds_out_file.name)
    else:
        print(f'za reported an error:\n{out.decode("utf-8")}')

    # Delete the temporary files and return the new xarray dataset
    ds_in_file.close()
    ds_out_file.close()
    if not out:
        return(ds_out)
    return(None)

def pop_add_cyclic(ds):
    
    nj = ds.TLAT.shape[0]
    ni = ds.TLONG.shape[1]

    xL = int(ni/2 - 1)
    xR = int(xL + ni)

    tlon = ds.TLONG.data
    tlat = ds.TLAT.data
    
    tlon = np.where(np.greater_equal(tlon, min(tlon[:,0])), tlon-360., tlon)    
    lon  = np.concatenate((tlon, tlon + 360.), 1)
    lon = lon[:, xL:xR]

    if ni == 320:
        lon[367:-3, 0] = lon[367:-3, 0] + 360.        
    lon = lon - 360.
    
    lon = np.hstack((lon, lon[:, 0:1] + 360.))
    if ni == 320:
        lon[367:, -1] = lon[367:, -1] - 360.

    #-- trick cartopy into doing the right thing:
    #   it gets confused when the cyclic coords are identical
    lon[:, 0] = lon[:, 0] - 1e-8

    #-- periodicity
    lat = np.concatenate((tlat, tlat), 1)
    lat = lat[:, xL:xR]
    lat = np.hstack((lat, lat[:,0:1]))

    TLAT = xr.DataArray(lat, dims=('nlat', 'nlon'))
    TLONG = xr.DataArray(lon, dims=('nlat', 'nlon'))
    
    dso = xr.Dataset({'TLAT': TLAT, 'TLONG': TLONG})

    # copy vars
    varlist = [v for v in ds.data_vars if v not in ['TLAT', 'TLONG']]
    for v in varlist:
        v_dims = ds[v].dims
        if not ('nlat' in v_dims and 'nlon' in v_dims):
            dso[v] = ds[v]
        else:
            # determine and sort other dimensions
            other_dims = set(v_dims) - {'nlat', 'nlon'}
            other_dims = tuple([d for d in v_dims if d in other_dims])
            lon_dim = ds[v].dims.index('nlon')
            field = ds[v].data
            field = np.concatenate((field, field), lon_dim)
            field = field[..., :, xL:xR]
            field = np.concatenate((field, field[..., :, 0:1]), lon_dim)       
            dso[v] = xr.DataArray(field, dims=other_dims+('nlat', 'nlon'), 
                                  attrs=ds[v].attrs)


    # copy coords
    for v, da in ds.coords.items():
        if not ('nlat' in da.dims and 'nlon' in da.dims):
            dso = dso.assign_coords(**{v: da})
                
            
    return dso

def label_plots(fig, axs, xoff=-0.04, yoff=0.02):
    alp = [chr(i).upper() for i in range(97,97+26)]
    for i, ax in enumerate(axs):    
        p = ax.get_position()
        x = p.x0 + xoff
        y = p.y1 + yoff
        fig.text(
            x, y , f'{alp[i]}',
            fontsize=14,
            fontweight='semibold'
        ) 