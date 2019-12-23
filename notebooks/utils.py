import os
import subprocess
import tempfile
import xarray as xr

def zonal_mean_via_fortran(ds, var):
    ds_in_file = tempfile.NamedTemporaryFile()
    ds_out_file = tempfile.NamedTemporaryFile()
    ds.to_netcdf(ds_in_file.name)
    ds.to_netcdf('/glade/scratch/mlevy/xarray_tonc.nc')
    za_exe = os.path.join(os.path.sep,
                          'glade',
                          'u',
                          'home',
                          'klindsay',
                          'bin',
                          'zon_avg',
                          'za')
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
    za_call = [za_exe,
               '-v', var,
               '-grid_file', grid_file,
               '-kmt_file', grid_file,
               '-O', '-o', ds_out_file.name, # -O overwrites existing file, -o gives file name
               ds_in_file.name]
    proc = subprocess.Popen(za_call, stdout=subprocess.PIPE)
    (out, err) = proc.communicate()
    print(out.decode('utf-8'))
    ds_out = xr.open_dataset(ds_out_file.name)
    ds_in_file.close()
    ds_out_file.close()
    return(ds_out)

