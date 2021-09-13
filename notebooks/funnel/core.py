import os
from glob import glob
import shutil
import pprint

from itertools import product

import traceback
import warnings
import json 
import yaml

import intake

import pandas as pd
import dask
import xarray as xr

from toolz import curry

from . config import cache_catalog_dir, cache_catalog_prefix, cache_format
from . registry import _DERIVED_VAR_REGISTRY, _QUERY_DEPENDENT_OP_REGISTRY

os.makedirs(cache_catalog_dir, exist_ok=True)


class Collection(object):
    """
    Catalog + Query + Operation = Unique Dataset
    
    Extend an intake-esm catalog with the ability to:
      (1) Compute derived variables combining multipe catalog entries;
      (2) Apply postprocessing operations to returned datasets;
      (3) Optionally cache the results and curate a discrete catalog of processed assets; 
          skip computation on subsequent calls for existing cataloged assets.

    Parameters
    ----------
    
    name : str
      Name of this collection.
      
    catalog : str
      JSON file defining `intake-esm` catalog
      
    query : dict
      Dictionary of keyword, value pairs defining a query subsetting the `catalog`.
          
    postproccess : callable, list of callables, optional
      Call these functions on each dataset following aggregation.
    
    postproccess_kwargs: dict, list of dict's, optional
       Keyword arguments to `postprocess` functions.
       **Note** These arguments must be serializable to yaml.
    
    persist : bool, optional
       Persist datasets to disk.
       
    cache_dir : str, optional
       Directory for storing cache files.
             
    kwargs : dict
       Defaults for keyword arguments to pass to `intake_esm.core.esm_datastore.to_dataset_dict`.
       (Note: these can be updated/overridden for each call to self.dsets below.)
    """
    def __init__(self, 
                 name,
                 esm_collection_json, 
                 query, 
                 postproccess=[], 
                 postproccess_kwargs=[], 
                 persist=False,
                 cache_dir='.', 
                 **kwargs,
                ):
       
        self.name = name
    
        # TODO: what should we do if "variable" is in the query?
        #       for now, just remove it
        variable = query.get('variable', None)
        if variable is not None:            
            warnings.warn('removed "variable" key from query')
    
        # TODO: accept multiple catalogs and concatenate them into one?
        self.catalog = intake.open_esm_datastore(esm_collection_json).search(**query)
        
        # make the assumption that the catalog.keys() provides a means of uniquely
        # identifying datasets within the catalog. This assumption is valid if
        # the user only constructs searches using the `groupby_attrs` columns,
        # so here we check to ensure this is the case.
        # TODO: an extension of this framework could strip out the `groupby_attrs`
        #       fields from the query and track the other keys of the query.
        #       For example, the user might include `date_range` in their query,
        #       which could truncate the time axis of the resulting datasets---
        #       but wouldn't change the keys: i.e. datasets with the same "key" might
        #       differ based on the fields of the `query` other than the `groupby_attrs`. 
        with open(esm_collection_json) as fid:
            catalog_def = json.load(fid)
        groupby_attrs = catalog_def['aggregation_control']['groupby_attrs']        
        assert not (set(query.keys()) - set(groupby_attrs)), (
            f'queries can only include the following fields: {groupby_attrs}'
        )
                    
        # setup cache info
        assert cache_format in ['nc', 'zarr'], f'unsupported format {cache_format}'             
        self._format = cache_format        
        self.persist = persist        
        self.cache_dir = cache_dir         
        if self.persist and not os.path.exists(cache_dir):
            os.makedirs(cache_dir)
        self._open_cache_kwargs = dict() # TODO: get this from config?
        if self._format == 'zarr' and 'consolidated' not in self._open_cache_kwargs:
            self._open_cache_kwargs['consolidated'] = True
        
        self._to_dsets_kwargs_defaults = kwargs
        
        # setup operators attrs
        if not isinstance(postproccess, list):
            postproccess = [postproccess]
        self.operators = postproccess
 
        if not postproccess_kwargs:
            self.ops_kwargs = [{} for op in self.operators]
        else:
            if not isinstance(postproccess_kwargs, list):
                postproccess_kwargs = [postproccess_kwargs]            
            assert len(postproccess_kwargs) == len(postproccess), 'mismatched ops/ops_kwargs'
            self.ops_kwargs = postproccess_kwargs
        
        # pull out "prepocess"
        preprocess_name = None if 'preprocess' not in kwargs else kwargs['preprocess'].__name__ 
        
        # build origins dict
        # TODO: save more about the function, including the source code 
        #       rather than simply the name
        self.origins_dict = dict(
            name=name,
            esm_collection=esm_collection_json, 
            preprocess=preprocess_name,
            operators=[op.__name__ for op in self.operators],
            operator_kwargs=self.ops_kwargs,
        )                

        # assemble database of existing cache's
        self._assemble_cache_db()
        
    def _assemble_cache_db(self):
        """loop over yaml files in cache dir; find ones that match"""
        self._cache_files = {}        
        if not self.persist:
            return
        
        cache_id_files = self._find_cache_id_files()
        
        if cache_id_files:
            for file in cache_id_files:
                try:
                    with open(file, 'r') as fid:
                        cache_id = yaml.load(fid, Loader=yaml.Loader)
                except: # TODO: be more informative here
                    print(f'cannot read cache: {file}...skipping')
                    print('skipping')
                    
                cache_origins_dict = {
                    k: cache_id[k] 
                    for k in self.origins_dict.keys()
                    if k in cache_id
                }
                if cache_origins_dict == self.origins_dict:
                    if not os.path.exists(cache_id['asset']):
                        print(f'missing: {cache_id["asset"]}')
                        continue
                    token = self._gen_cache_token(
                        cache_id['key'],                        
                        cache_id['variable'],
                    )
                    self._cache_files[token] = cache_id['asset']

    def _catalog_subset(self, variable, prefer_derived, refine_query):
        """apply query to catalog, determine if variable is derived"""
        catalog_subset_var = self.catalog.search(variable=variable, **refine_query)            

        # check for variable in derived registry
        is_derived = variable in _DERIVED_VAR_REGISTRY

        if len(catalog_subset_var):
            if is_derived and not prefer_derived:
                warnings.warn(
                    f'found variable "{variable}" in catalog and derived_var registry'
                )
                is_derived = False                
        else:
            if not is_derived:
                raise ValueError(f'variable not found {variable}')
        
        if is_derived:
            derived_var_obj = _DERIVED_VAR_REGISTRY[variable]
            query_vars = derived_var_obj.dependent_vars
            catalog_subset_var = self.catalog.search(variable=query_vars, **refine_query)
        
        return catalog_subset_var, is_derived
    
    def to_dataset_dict(self, variable, clobber=False, 
                        prefer_derived=False, refine_query={}, **kwargs):
        """
        Get dataset_dicts for assets.
        
        Parameters
        ----------
        
        variable : str or list
          The variable or list of variables to read and process.
          
        prefer_derived : bool, optional (default=False)
           Specify whether to use a "derived variable" regardless if 
           a variable by the same name is found in the catalog.
    
        refine_query : dict, optional
          Additional arguments to refine catalog search.
          
        clobber : bool, optional (default=False)
          Remove and recreate any existing cache.
          
        kwargs : dict, optional 
          Keyword arguments to pass to `intake_esm.core.esm_datastore.to_dataset_dict`
        
        """
       
        if isinstance(variable, list):
            dsets_list = [
                self.to_dataset_dict(v, clobber, prefer_derived, refine_query, **kwargs) 
                for v in variable
            ]
            keys = list(set([k for dsets in dsets_list for k in dsets.keys()]))
            dsets_merged = {}
            for key in keys:
                ds_list = []
                for dsets in dsets_list:                    
                    if key in dsets:
                        ds_list.append(dsets[key])
                dsets_merged[key] = xr.merge(ds_list, compat='override')
            return dsets_merged        
        
        catalog_subset_var, is_derived = self._catalog_subset(variable, prefer_derived, refine_query)        
        catalog_keys = catalog_subset_var.keys()
        key_info = intake_esm_get_keys_info(catalog_subset_var)

        dsets = {}
        for key in catalog_keys:
            # TODO: check for lock file, wait if present        
            if self._cache_exists(key, variable, clobber):        
                dsets[key] = self._read_cache(key, variable)
            else:                
                to_dsets_kwargs = kwargs.copy()
                to_dsets_kwargs.update(self._to_dsets_kwargs_defaults)            
                
                # TODO: set lock file                
                #       optionally spin up a cluster here
                dsets[key] = self._generate_dset(
                    key, variable, is_derived, catalog_subset_var, key_info, **to_dsets_kwargs
                )
                # TODO: release lock
        return dsets
    
    def _generate_dset(self, key, variable, is_derived, catalog_subset_var, key_info, **kwargs):
        """Do the computation necessary to make `dsets`"""
        
        dset = catalog_subset_var.search(**key_info[key]).to_dataset_dict(**kwargs)      
        assert len(dset.keys()) == 1, (
            f'expecting a single key ({key})\nfound the following: {dset.keys()}'
        )
        assert list(dset.keys())[0] == key, (
            f'mismatch in key:\nexpecting {key}\ngot: {dset.keys()}'
        )
        
        _, ds = dset.popitem()
        
        if is_derived:
            derived_var_obj = _DERIVED_VAR_REGISTRY[variable]
            ds = derived_var_obj(ds)

        for op, kw in zip(self.operators, self.ops_kwargs):
            if hash(op) in _QUERY_DEPENDENT_OP_REGISTRY:
                op_obj = _QUERY_DEPENDENT_OP_REGISTRY[hash(op)] 
                query_dict = dict(**key_info[key])
                query_dict['variable'] = variable
                ds = op_obj(
                    ds, 
                    query_dict=query_dict,
                    **kw,
                )
            else:
                ds = op(ds, **kw)

        if self.persist:
            self._make_cache(ds, key, variable)

        return ds
    
    def _make_cache(self, ds, key, variable):
        """write cache file"""
        
        token = self._gen_cache_token(key, variable)
        
        cache_id_dict = self.origins_dict.copy()
        cache_id_dict['key'] = key
        cache_id_dict['variable'] = variable
        cache_id_dict['asset'] = self._gen_cache_file_name(key, variable)        
        if self._format == 'nc':
            ds.to_netcdf(cache_id_dict['asset'])            
        elif self._format == 'zarr':
            ds.to_zarr(cache_id_dict['asset'], mode='w', consolidated=True)

        cache_id_file = self._gen_cache_id_file_name(key, variable)            
        with open(cache_id_file, 'w') as fid:
            yaml.dump(cache_id_dict, fid)

        self._cache_files[token] = cache_id_dict['asset']
                    
    def _read_cache(self, key, variable):            
        """read cache files"""
        token = self._gen_cache_token(key, variable)
        asset = self._cache_files[token]        
        if self._format == 'nc':
            with xr.open_dataset(asset, **self._open_cache_kwargs) as ds:
                return ds
        elif self._format == 'zarr':
            with xr.open_zarr(asset, **self._open_cache_kwargs) as ds:
                return ds

    def _cache_exists(self, key, variable, clobber):        
        """determine if cache files exist (or clobber them)"""
        # assemble database of existing cache's
        self._assemble_cache_db()        
        
        token = self._gen_cache_token(key, variable)                
        if token not in self._cache_files:
            return False
        
        asset = self._cache_files[token]            
        if os.path.exists(asset):
            if clobber:
                self._remove_asset(asset)
                return False
            else:
                return True
        else:
            return False
        
    def _remove_asset(self, asset):
        """delete asset from disk"""
        if not os.path.exists(asset):
            return
        if self._format == 'nc':
            os.remove(asset)
        elif self._format == 'zarr':
            shutil.rmtree(asset)
    
    def _gen_cache_token(self, key, variable):
        return dask.base.tokenize(self.origins_dict, key, variable)
        
    def _gen_cache_file_name(self, key, variable):
        """generate a file cache name"""
        # TODO: accept a user-provided callable to generate human-readable
        #       file name
        token_key = self._gen_cache_token(key, variable)
        return f'{self.cache_dir}/{token_key}.{self._format}'
    
    def _gen_cache_id_file_name(self, key, variable):
        """generate a unique cache file name"""
        token_key = self._gen_cache_token(key, variable)
        return f'{cache_catalog_dir}/{cache_catalog_prefix}-{token_key}.yml'
        
    def _find_cache_id_files(self):
        return sorted(glob(f'{cache_catalog_dir}/{cache_catalog_prefix}-*.yml'))        

    def __repr__(self):
        return pprint.pformat(self.origins_dict, indent=2, width=1, compact=True)
        
class derived_var(object):
    """
    Support computation of variables that depend on multiple variables, 
    i.e., "derived vars"
    """
    def __init__(self, dependent_vars, func):
        self.dependent_vars = dependent_vars
        self._callable = func                
        
    def __call__(self, ds, **kwargs):
        """call the function to compute derived var"""
        self._ensure_variables(ds)
        return self._callable(ds, **kwargs)
            
    def _ensure_variables(self, ds):
        """ensure that required variables are present"""
        missing_var = set(self.dependent_vars) - set(ds.variables)
        if missing_var:
            raise ValueError(f'Variables missing: {missing_var}')


class query_dependent_op(object):
    """
    Support calling functions that depend on the values of the query.
    """
    def __init__(self, func, query_keys):
        self._callable = func
        self._query_keys = query_keys
    
    def __call__(self, ds, query_dict, **kwargs):
        """call function with query keys added to keyword args"""
        kwargs.update({k: query_dict[k] for k in self._query_keys})
        return self._callable(ds, **kwargs)
    
    
@curry
def register_derived_var(func, varname, dependent_vars):
    """register a function for computing derived variables"""
    if varname in _DERIVED_VAR_REGISTRY:
        warnings.warn(
            f'overwriting derived variable "{varname}" definition'
        )

    _DERIVED_VAR_REGISTRY[varname] = derived_var(
        dependent_vars, func,
    )    
    return func


@curry
def register_query_dependent_op(func, query_keys):
    """register a function for computing derived variables"""
    func_hash = hash(func)
    if func_hash in _QUERY_DEPENDENT_OP_REGISTRY:
        warnings.warn(
            f'overwriting query dependent operator "{func.__name__}" definition'
        )

    _QUERY_DEPENDENT_OP_REGISTRY[func_hash] = query_dependent_op(
        func, query_keys,
    )    
    return func


def intake_esm_get_keys_info(cat):
    """return a dictionary with the values of components of the keys in an
       intake catalog
       
       Example:
         key_info = {
           'experiment': '20C',
           'component': 'ocn',
           'stream': 'pop.h',
           'member_id': 1,
         }
    """    
    
    # generate a list of lists with all possible values of each groupby_attr
    iterables = [
        cat.unique(columns=key)[key]['values'] 
        for key in cat.groupby_attrs
    ]    

    # generate a dictionary of keys with the values of its attributes
    key_info = {}    
    for values in product(*iterables):
        key = cat.sep.join([str(v) for v in values])
        if key in cat.keys():
            key_info[key] = {k: values[i] for i, k in enumerate(cat.groupby_attrs)}
    return key_info


def to_intake_esm():
    """generate an intake-esm data catalog from funnel collections"""
    
    catalog_csv_file = f'{cache_catalog_dir}/collection-summary.csv.gz'
    catalog_json_file = f'{cache_catalog_dir}/collection-summary.json'
    
    files = sorted(glob(f'{cache_catalog_dir}/*.yml'))
    data = {}
    for f in files:
        with open(f) as fid:
            data[f] = yaml.unsafe_load(fid)

    # assume that there is a *single* esm_collection
    # this could be extended, but supporting multiple collections
    # raises all sorts of questions about validation
    esm_collection = [v['esm_collection'] for v in data.values()]
    assert len(set(esm_collection)) == 1
    catalog = intake.open_esm_datastore(esm_collection[0])
    
    # generate a dictionary of the key info dictionaries for each catalog
    groupby_attrs_values = intake_esm_get_keys_info(catalog)
    first_key = list(groupby_attrs_values.keys())[0]
    columns = list(groupby_attrs_values[first_key].keys()) + ['variable', 'name', 'path']
    
    lines = []    
    for f in files:        
        column_data = dict(**groupby_attrs_values[data[f]['key']])
        column_data['variable'] = data[f]['variable']
        column_data['name'] = data[f]['name']
        column_data['path'] = data[f]['asset']
        lines.append(column_data)
    
    df = pd.DataFrame(lines)        
    assert set(df.columns) == set(columns), 'mismatch in expected columns'
        
    # modify the json
    with open(esm_collection[0]) as fid:
        catalog_def = json.load(fid)

    catalog_def['catalog_file'] = catalog_csv_file
    catalog_def['attributes'] = [{'column_name': k, 'vocabulary': ''} for k in columns]
    catalog_def['assets'] = {'column_name': 'path', 'format': cache_format}

    # ensure that all `groupby_attrs` are in the columns of the DataFrame
    assert all(
        [c in columns for c in catalog_def['aggregation_control']['groupby_attrs']]
    ), 'not all groupby attrs found in columns'
    
    # add `name` to `groupby_attrs`
    catalog_def['aggregation_control']['groupby_attrs'] += ['name']
    
    # filter the aggregations rules to ensure only existing columns are included
    catalog_def['aggregation_control']['aggregations'] = [
        d for d in catalog_def['aggregation_control']['aggregations']
        if d['attribute_name'] in columns
    ]    
    
    # persist
    df.to_csv(catalog_csv_file, index=False)
    
    with open(catalog_json_file, 'w') as fid:
        json.dump(catalog_def, fid)
        
    return catalog_json_file
    
    

    