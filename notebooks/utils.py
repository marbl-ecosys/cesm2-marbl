import os
import subprocess
import tempfile
import xarray as xr

def zonal_mean_via_fortran(ds, var, grid=None, three_ocean_regions=False):
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
    if three_ocean_regions:
        rm_file = os.path.join(os.path.sep,
                               'glade',
                               'p',
                               'cgd',
                               'oce',
                               'people',
                               'klindsay',
                               'oldcgdoce',
                               'REGION_MASK',
                               'new_REGION_MASK_gx1v6.nc')
        region_mask = ['-rmask_file', rm_file]
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

