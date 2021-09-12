import os
import shutil

import subprocess
import tempfile

import yaml

import matplotlib.pyplot as plt

import xarray as xr
import numpy as np

import pop_tools


def _gen_plotname(plot_name):
    """generate a name for plotting"""
    plot_basename, ext = os.path.splitext(plot_name)

    fig_key = None    
    if os.path.exists('_figure-order.yml'):
        with open('_figure-order.yml') as fid:
            fig_map = yaml.safe_load(fid)
        assert len(set(fig_map.values())) == len(fig_map.values()), (
            'non-unique figure names found in _figure-order.yml'
        )        
        for key, value in fig_map.items():
            if plot_basename == value:
                fig_key = key 
                break
                
    if 'TEST_PLOT_FIGURE_DIR' in os.environ:
        dirout_main = os.environ['TEST_PLOT_FIGURE_DIR']
    else:
        dirout_main = 'figures'    

    if fig_key is not None:
        plot_basename = f'{fig_key}-{plot_basename}'
        dirout = dirout_main
    else:
        dirout = f'{dirout_main}/misc'        
    os.makedirs(dirout, exist_ok=True)    
        
    return dirout, plot_basename


def savefig(plot_name):
    """Write figure"""
    
    dirout, plot_basename = _gen_plotname(plot_name)
    
    for ext in ['pdf', 'png']:
        plt.savefig(f'{dirout}/{plot_basename}.{ext}', 
                    dpi=300, 
                    bbox_inches='tight', 
                    metadata={'CreationDate': None})
        
        
def write_ds_out(dso, file_out):
    file_out = os.path.realpath(file_out)

    os.makedirs(os.path.dirname(file_out), exist_ok=True)

    if os.path.exists(file_out):
        shutil.rmtree(file_out)
    print('-'*30)
    print(f'Writing {file_out}')
    dso.info()
    print()
    dso.to_zarr(file_out);


def compute_kmt(da):
    """compute KMT based on missing values"""
    nk = len(da.z_t)
    
    KMT = np.zeros(da.shape[-2:]).astype(int)

    #-- where surface is missing, KMT = 0, else full depth
    KMT = np.where(np.isnan(da.data[0, :, :]), 0, nk)

    #   where level k is missing: KMT = k, i.e. the level above in 1-based indexing
    for k in range(1, nk):
        KMT = np.where(np.isnan(da.data[k, :, :]) & (KMT > k), k, KMT)

    return xr.DataArray(KMT, dims=('nlat', 'nlon'), name='KMT')


def zonal_mean_via_fortran(ds, var=None, grid=None, region_mask=None, replace_kmt=False):
    """
    Write ds to a temporary netCDF file, compute zonal mean for
    a given variable based on Keith L's fortran program, read
    resulting netcdf file, and return the new xarray dataset

    If three_ocean_regions=True, use a region mask that extends the
    Pacific, Indian, and Atlantic to the coast of Antarctica (and does
    not provide separate Arctic Ocean, Lab Sea, etc regions)
    """
    if replace_kmt and (var is None or ',' in var):
        raise ValueError('if "replace_kmt" is True, a single "var" must be specified.')
        
    ds_in_file = tempfile.NamedTemporaryFile(suffix='.nc')
    ds_out_file = tempfile.NamedTemporaryFile(suffix='.nc')
    
    ds = ds.copy()
    ds.attrs = {} # for some reason, za does not like file attrs---perhaps "coordinates"?
    ds.to_netcdf(ds_in_file.name)

    za_exe = '/glade/u/home/klindsay/bin/zon_avg/za'

    grid_file = None
    rmask_file = None
    
    if grid is not None:
        grid = pop_tools.get_grid(grid)
        
        grid_file = tempfile.NamedTemporaryFile(suffix='.nc')
        grid_file_name = grid_file.name
        if replace_kmt:
            grid['KMT'] = compute_kmt(ds[var])
        
        grid.to_netcdf(grid_file_name)
        
    else:
        # Assume xarray dataset contains all needed fields
        grid_file_name = ds_in_file.name

    # Set up the call to za with correct options
    za_call = [za_exe]
    if var is not None:
        za_call += ['-v', var] 
            
    if region_mask is not None:
        rmask_file = tempfile.NamedTemporaryFile(suffix='.nc')
        region_mask.to_netcdf(rmask_file.name)
        za_call += ['-rmask_file', rmask_file.name]

    za_call += [
        '-grid_file', grid_file_name,
        '-kmt_file', grid_file_name,
        '-O', '-o', ds_out_file.name, # -O overwrites existing file, -o gives file name
        ds_in_file.name
    ]
    # Use subprocess to call za, allows us to capture stdout and print it
    proc = subprocess.Popen(za_call, stdout=subprocess.PIPE)
    (out, err) = proc.communicate()

    subprocess.check_call(['cp', '-v', ds_in_file.name, f'{os.environ["TMPDIR"]}/za-in.nc'])
    subprocess.check_call(['cp', '-v', grid_file_name, f'{os.environ["TMPDIR"]}/za-grid.nc'])
    subprocess.check_call(['cp', '-v', rmask_file.name, f'{os.environ["TMPDIR"]}/za-rmask.nc'])        
    
    if not out:
        # Read in the newly-generated file
        print('za ran successfully, writing netcdf output')
        ds_out = xr.open_dataset(ds_out_file.name)
    else:
        print(f'za reported an error:\n{out.decode("utf-8")}')
        print(za_call)
        return 
    
    # clean up
    ds_in_file.close()
    ds_out_file.close()        
    if grid_file is not None:
        grid_file.close()
    if rmask_file is not None:
        rmask_file.close()
        
    return ds_out    


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

def label_plots(fig, axs, xoff=-0.04, yoff=0.02, start=0):
    
    alp = [chr(i).upper() for i in range(97,97+26)][start:]
    
    for i, ax in enumerate(axs):    
        p = ax.get_position()
        x = p.x0 + xoff
        y = p.y1 + yoff
        fig.text(
            x, y , f'{alp[i]}',
            fontsize=14,
            fontweight='semibold'
        ) 
        
        
        
def get_pop_region_mask_za(mask_type='3d', grid_name='POP_gx1v7',):
    """return a region mask for zonal averaging"""
    mask3d = pop_tools.region_mask_3d(grid_name, mask_name='Pacific-Indian-Atlantic')
    nregion = len(mask3d.region)
    
    if mask_type.lower() == '3d':
        return mask3d

    elif mask_type.lower() == '2d':
        mask2d = xr.full_like(
            mask3d.isel(region=0), fill_value=0, dtype=np.int32
        )
        for i in range(1, nregion): # skip first index because "za" puts the global field in there
            mask2d = xr.where(mask3d.isel(region=i)==1, i, mask2d)
        mask2d.name = 'REGION_MASK'
        return mask2d
    raise ValueError(f'unknown mask type: {mask_type}\nexpecting either "2d" or "3d"')
    
    

def subplot_col_labels(axs, col_labels, xoff=0.):
    assert len(axs) == len(col_labels)
    for ax, col_label in zip(axs, col_labels):
        ax.annotate(col_label,
                    xy=(np.mean(ax.get_xlim())+xoff, np.max(ax.get_ylim()) - np.diff(ax.get_ylim())*0.12), 
                    xytext=(np.mean(ax.get_xlim())+xoff, np.max(ax.get_ylim()) + np.diff(ax.get_ylim())*0.12), 
                    fontsize='14', fontweight='bold', ha='center', va='center')

def subplot_row_labels(axs, row_labels, yoff=0., xoff=0.):    
    assert len(axs) == len(row_labels)
    for ax, row_label in zip(axs, row_labels):
        ax.annotate(row_label, xy=(0+xoff, 0.5+yoff), xytext=(-ax.yaxis.labelpad-12+xoff, 0+yoff),
                    xycoords=ax.yaxis.label, textcoords='offset points',
                    rotation=90,
                    fontsize='14', 
                    fontweight='bold', ha='center', va='center')    
        

def get_ClusterClient():
    import dask
    from dask_jobqueue import PBSCluster
    from dask.distributed import Client
    
    USER = os.environ['USER']
    
    cluster = PBSCluster(
        cores=1,
        memory='25GB',
        processes=1,
        queue='casper',
        local_directory=f'/glade/scratch/{USER}/dask-workers',
        log_directory=f'/glade/scratch/{USER}/dask-workers',
        resource_spec='select=1:ncpus=1:mem=25GB',
        project='NCGD0011',
        walltime='06:00:00',
        interface='ib0',)

    dask.config.set({
        'distributed.dashboard.link':
        'https://jupyterhub.hpc.ucar.edu/stable/user/{USER}/proxy/{port}/status'
    })
    client = Client(cluster)
    return cluster, client

