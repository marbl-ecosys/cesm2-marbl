import os

import numpy as np
import pandas as pd
import xarray as xr

import ESMF

USER = os.environ['USER']

os.environ['CESMDATAROOT'] = f'/glade/scratch/{USER}/inputdata'
import pop_tools 

path_to_here = os.path.dirname(os.path.realpath(__file__))

regrid_dir = f'/glade/work/{USER}/adhoc-regridding'
os.makedirs(regrid_dir, exist_ok=True)


def _ensure_grid_file(grid_name, clobber):
    """ensure that grid file exists"""
    
    grid_file = f'{regrid_dir}/{grid_name}.nc' 
    if os.path.exists(grid_file) and not clobber:
        return grid_file
        
    # generate file if needed
    if grid_name in ['POP_gx1v6', 'POP_gx1v7', 'POP_gx3v7',]:
        dso = pop_tools.get_grid(grid_name, scrip=True)            
    else:
        raise ValueError('unknown grid')   
    
    dso.to_netcdf(grid_file)
    return grid_file


def _esmf_pop_grid(grid_name, clobber=False):
    """instantiate an ESMF grid object"""
    return ESMF.Grid(
        filename=_ensure_grid_file(grid_name, clobber),
        filetype=ESMF.api.constants.FileFormat.SCRIP,
        add_corner_stagger=True,
        pole_kind=[ESMF.api.constants.PoleKind.NONE, ESMF.api.constants.PoleKind.NONE], # is this the right choice?
    )


def _esmf_locstream(lon, lat):
    """instantiate an ESMF locstream object"""
    locstream = ESMF.LocStream(
        len(lon), coord_sys=ESMF.CoordSys.SPH_DEG,
    )
    locstream["ESMF:Lon"] = lon.astype(np.dtype('f8'))
    locstream["ESMF:Lat"] = lat.astype(np.dtype('f8'))

    return locstream

       
def open_datastream(obs_name):
    """open raw dataset"""
    
    filename_dict = dict(
        dFe=f'{path_to_here}/dFe-database-2021-04-20.csv',
        DOM=f'{path_to_here}/DOMobs.csv',
        test=f'{path_to_here}/dfe-test.csv',
    )

    try:
        filename = filename_dict[obs_name]
    except:
        raise ValueError(f'unknown obs name {obs_name}')

    return pd.read_csv(filename, na_values=-999.).dropna(axis=0, how='all')


@pd.api.extensions.register_dataframe_accessor('obs_stream')
class obs_datastream:
    def __init__(self, pandas_obj):
        self._validate(pandas_obj)
        self._obj = pandas_obj

    @staticmethod
    def _validate(obj):
        """verify the requried columns are present"""
        for field in ['lat', 'lon', 'depth']:
            if field not in obj.columns:
                raise AttributeError(f"Must have '{field}' column.")

    def add_model_field(self, da_in, model_grid=None, field_name=None, method='bilinear'):
        """return a DataFrame with obs and model"""
        
        # determine dimensions
        if da_in.dims == ('z_t', 'nlat', 'nlon'):
            nk, nj, ni = da_in.shape
        elif da_in.dims == ('nlat', 'nlon'):
            nk = 0
            nj, ni = da_in.shape
        else:
            raise ValueError(f'dimensions not supported: {da_in.dims}')

         
        # get model grid
        if model_grid is None:
            if (nj, ni) == (116, 100):
                model_grid = 'POP_gx3v7'
            elif (nj, ni) == (384, 320):
                model_grid = 'POP_gx1v7'
            else:
                raise ValueError(f'cannot infer model grid: {da_in.dims}')
        
        grid = _esmf_pop_grid(model_grid)

        # define locstream
        df = self._obj
        n_obs = len(df)
        locstream = _esmf_locstream(df.lon.values, df.lat.values)

        # set up remapping TODO: precompute and save regrid?
        srcfield = ESMF.Field(grid, name='srcfield')
        dstfield = ESMF.Field(locstream, name='dstfield')

        method_dict = dict(
            bilinear=ESMF.RegridMethod.BILINEAR,
            nearest=ESMF.RegridMethod.NEAREST_STOD,
        )
        try:
            ESMF_RegridMethod = method_dict[method]
        except:
            raise ValueError(f'unkown method {method}')

        regrid = ESMF.Regrid(
            srcfield, dstfield,
            regrid_method=ESMF_RegridMethod,
            unmapped_action=ESMF.UnmappedAction.ERROR,
        )
        
        if field_name is None:
            field_name = da_in.name
            i = 1
            while field_name in df:
                field_name = f'{da_in.name}_{i}'
                i += 1

        # 2D field
        if nk == 0:
            dstfield.data[...] = np.nan            
            srcfield.data[...] = da_in.data[:, :].T
            dstfield = regrid(srcfield, dstfield, zero_region=ESMF.Region.SELECT)
            df[field_name] = dstfield.data
        else:
            
            # TODO: this is a little clunky, would be better to simply do 3D interpolation
            da_out_columns = np.ones((nk, n_obs)) * np.nan 
            for k in range(nk):
                dstfield.data[...] = np.nan                
                srcfield.data[...] = da_in.data[k, :, :].T
                dstfield = regrid(srcfield, dstfield, zero_region=ESMF.Region.SELECT)
                da_out_columns[k, :] = dstfield.data

            dstfield_z = np.ones((n_obs)) * np.nan
            for n in range(n_obs):
                dstfield_z[n] = np.interp(df.depth.values[n]*1e2, da_in.z_t, da_out_columns[:, n])            
            
            df[field_name] = dstfield_z