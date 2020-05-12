import os
import shutil

import xarray as xr
import xpersist as xp

import intake
import intake_esm

from . config import path_to_here, project_tmpdir
from . import util
from . defaults import cdf_kwargs as cdf_kwargs_defaults
from . import pop


cache_rootdir_local = f'{path_to_here}/data'
cache_rootdir_tmpdir = project_tmpdir


class Component(object):
    """A class to enable provenance tracking"""
    
    def __init__(
        self, 
        model, 
        experiment, 
        variable,
        component='ocn', 
        stream='pop.h',
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
        return (
            f'{self.name}:\n'
            f'  catalog: {self.catalog_file}\n'
            f'  query: {self.query}'
        )
    
    @property
    def ds(self, cdf_kwargs=None):
        """get a dataset from catalog and query, but only if necessary"""        

        if self._ds is None:

            # should we do a `update` of cdf_kwargs_defaults or replacement?
            # here I have opted for replacement
            if cdf_kwargs is None:
                cdf_kwargs = cdf_kwargs_defaults
                
            dsets = self.catalog.to_dataset_dict(cdf_kwargs=cdf_kwargs)
            
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
        
    @property
    def _dataset_methods(self):
        """This is a dictionary mapping named operations to functions"""
        
        # TODO: make this user settable/extensible?
        #       select appropriate methods for different subclasses
        return dict(
            area_mean=dict(
                function=self._compute_area_mean,
                cache_rootdir=cache_rootdir_tmpdir,
                cache_format='zarr',
            ),
            timeseries_ann=dict(
                function=self._compute_timeseries_ann,
                cache_rootdir=cache_rootdir_tmpdir,
                cache_format='nc',
            ),
            time_mean=dict(
                function=self._compute_time_mean,
                cache_rootdir=cache_rootdir_tmpdir,
                cache_format='nc',
            ),
        )
    
    def get_dataset(self, dataset_name, persist=True,
                     clobber=False, **kwargs):
        """return a computed dataset"""
        
        dataset_name = dataset_name.replace('-', '_')        
        assert dataset_name in self._dataset_methods, (
            f'unrecognized dataset: {dataset_name}' 
        )
        compute_func = self._dataset_methods[dataset_name]['function']
        
        if persist:
            # cache file
            cache_path = f'{self._dataset_methods[dataset_name]["cache_rootdir"]}/{self.name}'
            cache_format = self._dataset_methods[dataset_name]["cache_format"]

            # TODO: test for kwargs that should not be serialized here, like arrays, etc.
            cache_name = '.'.join(
                [dataset_name] + [
                    f'{k}={v}'.replace('.', '_')
                    for k, v in kwargs.items()
                ]
            )
            
            cache_name = dataset_name
            cache_file = f'{cache_path}/{cache_name}.{cache_format}'

            # I think there is a bug in xpersist            
            if clobber and os.path.exists(cache_file):
                shutil.rmtree(cache_file)

            # generate xpersist partial            
            persist_ds = xp.persist_ds(
                name=cache_name, 
                path=cache_path, 
                trust_cache=True,
                format=cache_format,                
            )

            return persist_ds(compute_func)(**kwargs)

        else:
            return compute_func(**kwargs)
        
    def _compute_area_mean(self, normalize=True, region_mask_name=None):
        """compute area-weighted average"""

        ds = self.ds
        
        # TODO: extend with region mask definitions
        if region_mask_name is None:
            masked_area = ds.TAREA.where(self.ds.REGION_MASK>0)
        else:
            raise NotImplemented('region mask selection not available yet')

        # compute 
        xy_vars = [
            v for v in ds.data_vars 
            if not {'nlat', 'nlon'} - set(ds[v].dims)
        ]

        other_vars = set(ds.data_vars) - set(xy_vars) 

        with xr.set_options(keep_attrs=True):
            ds_xy_mean = (ds[xy_vars] * masked_area).sum(['nlat', 'nlon'])
            if normalize:
                ds_xy_mean = ds_xy_mean / masked_area.sum(['nlat', 'nlon'])

        # TODO: convert to output units                    
        # update units
        if not normalize:
            for v in xy_vars:
                ds_xy_mean[v].attrs['units'] = ' '.join([
                    ds[v].attrs['units'], masked_area.units
                ])

        # copy back other variables
        for v in other_vars:
            ds_xy_mean[v] = ds[v]

        # copy back coords
        for v in ds.coords:
            if {'nlat', 'nlon'} - set(ds[v].dims):
                ds_xy_mean[v] = ds[v]

        return ds_xy_mean.compute()

    def _compute_timeseries_ann(self, normalize=True, region_mask_name=None):
        """compute an annual mean timeseries"""
        ds_xy_mean = self._compute_area_mean(normalize=normalize, 
                                             region_mask_name=region_mask_name)        
        return util.calc_ann_mean(ds_xy_mean)
    
    def _compute_time_mean(self, time_slice=None):
        """compute the time-mean"""
        if time_slice is None:
            ds = util.calc_time_mean(ds)
        else:
            print('computing time mean over:')
            print(self.ds.sel(time=time_slice).time)
            print()
            ds = util.calc_time_mean(self.ds.sel(time=time_slice))
        # TODO: convert to output units            
        return ds

    
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
        