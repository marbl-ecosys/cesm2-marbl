import os
import shutil

import pickle
from toolz import curry

import dask
import xarray as xr
import xpersist as xp

import intake
import intake_esm

from . config import path_to_here, project_tmpdir
from . import units
from . import util
from . import pop


class Component(object):
    """A class to enable provenance tracking"""
    
    def __init__(
        self, 
        model, 
        experiment, 
        variable,
        component='ocn', 
        stream='pop.h',
        cache_path=None,        
        cache_format='zarr',
        cdf_kwargs=None,
        **kwargs,
    ):
        # TODO: add support beyond POP        
        assert component == 'ocn'
        
        # TODO: determine the right catalog file
        self.model = model
        self.catalog_file = f'{path_to_here}/catalogs/campaign-cesm2-cmip6-timeseries.json'        
        
        # TODO: construct a name for this that can track catalog/query
        self.name = '-'.join([
            model,
            experiment,
            component,
            stream,
            variable,
        ])
        
        if cache_path is None:
            self.cache_path = project_tmpdir
        else:
            self.cache_path = cache_path
        self.cache_format = cache_format
        
        if cdf_kwargs is None:
            self.cdf_kwargs = dict(
                chunks=dict(time=180),
                decode_coords=False,
                decode_times=False,
            )
        
        # apply the query and get a dataset
        catalog = intake.open_esm_datastore(self.catalog_file, sep=':')

        if not isinstance(variable, list):
            variable = [variable]

        self.data_vars = variable
            
        # manage derived variables
        query_variables, derived_variables = manage_var_dep(
            variable_list=variable,
            defined_model_variables=list(catalog.df.variable.unique()),
        )
        
        self._derived_var_func = None
        if derived_variables:
            self._derived_var_func = defined_derived_variables[variable]['function']
        
        # set the query
        self.query = dict(
            experiment=experiment,
            component=component,
            stream=stream,
            variable=query_variables,
        )               
        self.catalog = catalog.search(**self.query)   
        self._ds_key = ':'.join([component, experiment, stream])
        self._ds = None

    def __repr__(self):
        """return string representation of object"""
        return (f'{self.name}:\n'
                f'  catalog: {self.catalog_file}\n'
                f'  query: {self.query}\n'
                f'  assets: {self.catalog}'
               )
    
    @property
    def ds(self):
        """get a dataset from catalog and query, but only if necessary"""        

        if self._ds is None:               
            dsets = self.catalog.to_dataset_dict(cdf_kwargs=self.cdf_kwargs)
            
            # TODO: should we support concatentation logic?
            assert len(dsets.keys()) == 1

            ds = dsets[self._ds_key]

            # TODO: need test
            # apply derivation function
            if self._derived_var_func is not None:
                ds = self._derived_var_func(ds)
                
            # prune and conform
            keep_vars = self.data_vars + ['time', ds.time.bounds]
            ds = ds[keep_vars]
            ds = util.fix_units(ds)
            ds = util.compute_time(ds)
            self._ds = ds.assign_coords(self._grid.variables) 
            
        return self._ds
        
    @property    
    def _grid(self):
        """get the grid"""        
        # TODO: extend to other components
        return pop.get_grid(self.model)    

    @curry
    def _persist_ds(
        compute_func, 
        dataset_name, 
    ):
        """call operator within a wrapper to persist a dataset with generated name"""
        
        # define the wrapper function
        def compute_wrapper(self, **kwargs):
            # get variables meant for local control
            persist = kwargs.pop('persist', True)
            clobber = kwargs.pop('clobber', False)

            if persist:
                
                # TODO: some of this probably belongs in xpersist
                # there should be a data registry that keeps track of the
                # files. Files should have nominally human-readable
                # names. The catalog, query and dask task graph 
                # should be saved in the registry and compared 
                # to determine if recompute is necessary. 
                
                # generate cache file name               
                token = dask.base.tokenize(dataset_name, kwargs)
                cache_name = f'{self.name}-{dataset_name}-{token}'                
                cache_path = f'{self.cache_path}/{self.name}'
                cache_file = f'{cache_path}/{cache_name}.{self.cache_format}'
                cache_db_entry = {
                    cache_file: dict(
                        method=dataset_name,
                        kwargs=kwargs,
                    )
                }                             
                
                # register dataset
                cache_db = {}
                if os.path.exists(cache_database_file):
                    with open(cache_database_file, 'rb') as fid:
                        cache_db = pickle.load(fid)

                if cache_file in cache_db:                   
                    # TODO: ensure that catalog/query/kwargs/task-graph match
                    pass
                else:
                    cache_db.update(cache_db_entry)
                    with open(cache_database_file, 'wb') as fid:
                        pickle.dump(cache_db, fid)
                
                # I think there is a bug in xpersist            
                if clobber and os.path.exists(cache_file):
                    shutil.rmtree(cache_file)
                
                # generate xpersist partial            
                xp_persist_ds = xp.persist_ds(
                    name=cache_name, 
                    path=cache_path, 
                    trust_cache=True,
                    format=self.cache_format,                
                )
                
                # call the function
                return xp_persist_ds(compute_func)(self, **kwargs)
            
            else:
                return compute_func(self, **kwargs)

        return compute_wrapper    
    
    @_persist_ds(dataset_name='area-mean')
    def area_mean(self, normalize=True, region_mask_name=None):
        """compute area-weighted average"""

        # TODO: extend with region mask definitions
        if region_mask_name is None:
            masked_area = self.ds.TAREA.where(self.ds.REGION_MASK>0)
        else:
            raise NotImplemented('region mask selection not available yet')

        # compute 
        xy_vars = [
            v for v in self.ds.data_vars 
            if not {'nlat', 'nlon'} - set(self.ds[v].dims)
        ]

        other_vars = set(self.ds.data_vars) - set(xy_vars) 

        with xr.set_options(keep_attrs=True):
            ds_xy_mean = (self.ds[xy_vars] * masked_area).sum(['nlat', 'nlon'])
            if normalize:
                ds_xy_mean = ds_xy_mean / masked_area.sum(['nlat', 'nlon'])

        # TODO: convert to output units                    
        # update units
        if not normalize:
            for v in xy_vars:
                ds_xy_mean[v].attrs['units'] = ' '.join([
                    self.ds[v].attrs['units'], masked_area.units
                ])

        # copy back other variables
        for v in other_vars:
            ds_xy_mean[v] = self.ds[v]

        # copy back coords
        for v in self.ds.coords:
            if {'nlat', 'nlon'} - set(self.ds[v].dims):
                ds_xy_mean[v] = self.ds[v]

        return ds_xy_mean.convert_units(self.data_vars).compute()

    @_persist_ds(dataset_name='timeseries-ann')
    def timeseries_ann(self, normalize=True, region_mask_name=None):
        """compute an annual mean timeseries"""
        ds_xy_mean = self.area_mean(persist=False,
                                    normalize=normalize,                                    
                                    region_mask_name=region_mask_name)
        return util.calc_ann_mean(ds_xy_mean)
    
    @_persist_ds(dataset_name='time-mean')
    def time_mean(self, time_slice=None):
        """compute the time-mean"""
        if time_slice is None:
            ds = util.calc_time_mean(ds)
        else:
            ds = util.calc_time_mean(self.ds.sel(time=time_slice))

        return ds.convert_units(self.data_vars).compute()

    
def manage_var_dep(variable_list, defined_model_variables):
    """Determine if a variable is written directly by 
       the model (and therefore just read-in)
       or derived from dependencies.
    """
    query_variables = []    
    derived_variables = []
    for v in variable_list:    
        if v in defined_model_variables:
            query_variables.append(v)
            
        elif v in defined_derived_variables:
            q, d = manage_var_dep(
                defined_derived_variables[v]['dependencies']
            )
            derived_variables.append(v)
            query_variables.extend(q)
        else:    
            raise ValueError(f'unknown variable {v}')    

    return (
        sorted(list(set(query_variables))), 
        sorted(list(set(derived_variables)))
    )


defined_derived_variables = {
    'Chl_surf': {
        'dependencies': ['diatChl', 'spChl', 'diazChl'],
        'function': pop.compute_chl_surf,
    },
}
        
    