import os
from glob import glob

import numpy as np
import xarray as xr
import urllib
import netCDF4 as nc

os.environ['CESMDATAROOT'] = '/glade/scratch/mclong/inputdata'
import pop_tools

import regrid_tools


cache_dir = f'/glade/p/cgd/oce/projects/cesm2-marbl/glodap-data'

known_products = [
    'GLODAPv1',
    'GLODAPv2.2016b_MappedClimatologies',
    # 'GLODAPv2_Mapped_Climatologies', # this appears to be on pressure surfaces
]

def _ensure_datafiles(product_name='GLODAPv2.2016b_MappedClimatologies'):
    """
    get data files from website and return dictionary
    
    product_name='GLODAPv2.2016b_MappedClimatologies'
    Variables returned = 'Cant', 'NO3', 'OmegaA', 'OmegaC', 'PI_TCO2', 'PO4', 
                          'TAlk', 'TCO2', 'oxygen', 'pHts25p0', 'pHtsinsitutp', 
                          'salinity', 'silicate', 'temperature', 
    
    Alternative to default:
    product_name='GLODAPv2_Mapped_Climatologies'
    Variables returned = 'OmegaAinsitu',  'OmegaCinsitu',  'nitrate',  'oxygen',  
                          'pHts25p0',  'pHtsinsitu',  'phosphate',  'salinity', 
                          'silicate',  'talk',  'tco2',  'theta',    
    """
    
    url = 'https://www.nodc.noaa.gov/archive/arc0107/0162565/2.2/data/0-data/mapped'    

    filename = (
        'GLODAPv2_Mapped_Climatology.tar.gz'
        if product_name == 'GLODAPv2_Mapped_Climatologies' 
        else f'{product_name}.tar.gz'
    )

    files = sorted(glob(f'{cache_dir}/{product_name}/*.nc'))
    if not files:
        os.makedirs(cache_dir, exist_ok=True)
        local_file = f'{cache_dir}/{filename}' 
        urllib.request.urlretrieve(f'{url}/{filename}', local_file)
        check_call(['gunzip', local_file])
        check_call(['tar', '-xvf', local_file.replace('.gz', ''), '-C', cache_dir])
        files = sorted(glob(f'{cache_dir}/{product_name}/*.nc'))

    return {f.split('.')[-2]: f for f in files}
    

def _ensure_datafiles_v1():
    """get data friles for GLODAPv1 (txt files, arrrgh) and return dictionary"""
    
    url = 'https://www.ncei.noaa.gov/data/oceans/ncei/ocads/data/0001644/GLODAP_gridded.data'
    subdir = [
        'Alk.data',
        'AnthCO2.data',
        'BkgC14.data',
        'BombC14.data',
        'C14.data',
        'CFC.data',
        'TCO2.data',
    ]
    product_name = 'GLODAPv1'
    os.makedirs(f'{cache_dir}/{product_name}', exist_ok=True)
    for d in subdir:
        if d == 'CFC.data':
            these_files = ['CFC-11.data.txt', 'CFC-12.data.txt', 'pCFC-11.data.txt', 'pCFC-12.data.txt']
        else:
            these_files = [f'{d}.txt']

        for f in these_files:
            local_file = f'{cache_dir}/{product_name}/{f}'
            if not os.path.exists(local_file):
                print(f'downloading {url}/{d}/{f}')
                urllib.request.urlretrieve(f'{url}/{d}/{f}', local_file)

        files = sorted(glob(f'{cache_dir}/{product_name}/*.data.txt'))

    return {os.path.basename(f).split('.')[-3].replace('-', ''): f for f in files}    


def _gen_v1_dataset(clobber=True):
    """convert *.data.txt to dataset"""
    
    netcdf_file = f'{cache_dir}/GLODAPv1/GLODAPv1.nc'
    if os.path.exists(netcdf_file) and clobber:
        os.remove(netcdf_file)
        
    if os.path.exists(netcdf_file):
        return xr.open_dataset(netcdf_file)
    else:
        obs_files = _ensure_datafiles_v1()

        attrs = dict(
            pCFC11=dict(long_name='pCFC-11', units='patm'),
            pCFC12=dict(long_name='pCFC-12', units='patm'),            
            Cant=dict(long_name='C$_{ant}$', units='µmol kg$^{-1}$'),
            DIC=dict(long_name='DIC', units='µmol kg$^{-1}$'),
            ALK=dict(long_name='Alkalinity', units='µmol kg$^{-1}$'),            
            Del14C=dict(long_name='$^{14}$C', units='permille'),                        
        )
        
        depth_edges = xr.DataArray(
            np.array([0, 5, 15, 25, 40, 60, 85, 110, 135, 175, 
                      225, 275, 350, 450, 550, 650, 750, 850, 950, 
                       1050, 1150, 1250, 1350, 1450, 1600, 1850, 2250, 
                       2750, 3250, 3750, 4250, 4750, 5250, 6500]).astype(np.float), 
            dims=('depth_edges',),
        )

        depth = xr.DataArray(
            np.array([0, 10, 20, 30, 50, 75, 100, 125, 150, 200,250,300,400,
                      500,600,700,800,900,1000,1100,1200,1300,1400,1500,1750,
                      2000,2500,3000,3500,4000,4500,5000,5500]).astype(np.float), 
            dims=('depth',),
        )

        dx = 1
        dy = 1
        lon0 = -180.
        lat = xr.DataArray(np.arange(-90. + dy / 2., 90., dy), dims=('lat'))
        lon = xr.DataArray(np.arange(lon0 + dx / 2., lon0 + 360., dx), dims=('lon'))
        ds = xr.Dataset(dict(
            lon=lon,
            lat=lat,
            depth=depth,
            depth_edges=depth_edges,
        ))

        def convert_str(string):
            if string == '-999':
                return np.nan
            else:
                try:
                    value = np.float(string)
                except:
                    print('conversion error')
                    print(string)
                    raise
                return value

        nk, nj, ni = len(depth), len(lat), len(lon)
        
        for varname, file_in in obs_files.items():
            if varname in ['BombC14']:
                continue

            data = np.ones((nk, nj, ni))
            with open(file_in, 'r') as fid:
                n = 0
                for k in range(nk): 
                    for i in range(ni):
                        line = fid.readline()
                        data[k, :, i] = [convert_str(v.strip()) for v in line.split(',')]                   
                        n += 1               

            ds[varname] = xr.DataArray(data, dims=('depth', 'lat', 'lon'))
        ds = ds.rename({
            'AnthCO2': 'Cant_v1', 
            'TCO2': 'DIC', 
            'Alk': 
            'ALK', 
            'C14': 'Del14C',
        })
        for v in ds.data_vars:
            if v in attrs:
                ds[v].attrs = attrs[v]
                
        ds.to_netcdf(netcdf_file)
        return ds


def open_glodap(product='GLODAPv1'):
    """return GLODAP dataset"""
    assert product in known_products
    if product == 'GLODAPv1':
        return _gen_v1_dataset()
    else:
        obs_files = _ensure_datafiles(product)
        ds_list = []
        for varname, file_in in obs_files.items():
            ds = xr.open_dataset(file_in)
            depth = 'Depth' if 'Depth' in ds else 'depth'
            ds_list.append(ds[[depth, varname]])
        ds = xr.merge(ds_list)
        ds = ds.rename({'Depth': 'depth'})
        ds = ds.rename({'depth_surface': 'depth'}).set_coords('depth')
        ds = ds.rename({'TAlk': 'ALK', 'TCO2': 'DIC', 'oxygen': 'O2',})            
        return ds
    

def open_glodap_pop_grid(product='GLODAPv1', model_grid='POP_gx1v7', method='bilinear'):
    """return GLODAP dataset"""
    assert product in known_products
    
    ds_src = open_glodap(product)
    dst_grid = regrid_tools.grid(model_grid, clobber=False)
    
    if product == 'GLODAPv1':
        src_grid = regrid_tools.grid('latlon_glodapv1', nx=360, ny=180, lon0=-180.)#, grid_imask=grid_imask)
    else:
        src_grid = regrid_tools.grid('latlon_glodapv2', nx=360, ny=180, lon0=20.5)#, grid_imask=grid_imask)
        
    regridder = regrid_tools.regridder(src_grid, dst_grid, method, clobber=False)
    
    spatial_vars = [v for v in ds_src.data_vars if ('lat', 'lon') == ds_src[v].dims[-2:]]
    non_spatial_vars = set(ds_src.data_vars) - set(spatial_vars)
    
    ds_dst_xy = regridder.regrid(
        ds_src[spatial_vars], renormalize=True, apply_mask=False,
    )
    
    for v in spatial_vars:
        ds_dst_xy[v].attrs = ds_src[v].attrs
        
    for v in non_spatial_vars:
        ds_dst_xy[v] = ds_src[v]
        
    return add_coords_regrid_vertical(ds_dst_xy)


def add_coords_regrid_vertical(ds_dst_xy, pop_grid='POP_gx1v7'):
    """perform vertical regridding"""    
    ydim = 'lat'
    xdim = 'lon'
    zdim = 'depth'
    
    ds_dst_xy = ds_dst_xy.assign_coords({zdim: ds_dst_xy[zdim] * 1e2}) # m --> cm
    
    grid_vars = ['TLONG', 'TLAT', 'TAREA', 'z_t', 'dz', 'KMT']
    ds_pop = pop_tools.get_grid(pop_grid)
    ds_dst = ds_pop[grid_vars].set_coords(grid_vars)
        
    for v in ds_dst_xy.data_vars:
        da = ds_dst_xy[v]
        if zdim not in da.dims:
            continue
        else:
            with xr.set_options(keep_attrs=True):
                da_out = da.interp(
                    coords={zdim: ds_pop.z_t},
                    method='linear',
                    assume_sorted=True,
                    kwargs={'bounds_error': False}
                )
                da_out = xr.where(np.isclose(da_out, 0., atol=1e-10), 0., da_out)
                da_out.encoding['_FillValue'] = nc.default_fillvals['f8']   
                da_out.encoding['coordinates'] = 'TLONG TLAT z_t'
        ds_dst[v] = da_out
        ds_dst[v].attrs = da.attrs
    
    for v in grid_vars:
        ds_dst[v].encoding['_FillValue'] = None
        
    return ds_dst.drop([zdim]).rename({ydim: 'nlat', xdim: 'nlon'})    