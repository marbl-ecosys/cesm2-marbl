import os
from jinja2 import Template
import yaml

import funnel as fn
import operators as ops
import variable_defs


def get_config_dict():
    """return the configuration dictionary with environment variables replaced"""
    
    with open('_config_calc.yml') as fid:
        config_dict_in = yaml.safe_load(fid)    
    
    required_keys = ['esm_collection', 'data_collections', 'cache_dir']
    for key in required_keys:
        assert key in config_dict_in, f'config missing {key}'
        
    config_dict = {}
    for key, value in config_dict_in.items():
        if isinstance(value, str):
            t = Template(value)
            config_dict[key] = t.render(env=os.environ)
        else:
            config_dict[key] = value
    return config_dict



config_dict = get_config_dict()
_collections_settings = config_dict['data_collections']
esm_collection_json = config_dict['esm_collection']


@fn.register_query_dependent_op(
    query_keys=['experiment'],
)
def _mean_time_for_experiment(ds, experiment):
    """compute the mean over time"""
    time_range = _collections_settings['epoch_mean']['experiment'][experiment]
    sel_dict = dict(time=slice(time_range[0], time_range[1]))
    return ops.mean_time(
        ds, 
        sel_dict=sel_dict,
    )


def epoch_mean(query, name='epoch_mean', center_time=True, time_range=None):
    """Instantiate a `funnel.Collection` object for computing epoch means."""
    if center_time:
        postproccess = [ops.center_time]
        postproccess_kwargs = [{}]
    else:
        postproccess = []
        postproccess_kwargs = []    
           
    if time_range is None:
        postproccess += [_mean_time_for_experiment] 
        postproccess_kwargs += [{}]    
    else:
        sel_dict = dict(time=slice(time_range[0], time_range[1]))
        postproccess += [ops.mean_time] 
        postproccess_kwargs += [dict(sel_dict=sel_dict)]            
        
    return fn.Collection(
        name=name,
        esm_collection_json=esm_collection_json,
        postproccess=postproccess,  
        postproccess_kwargs=postproccess_kwargs,
        query=query,
        cache_dir=config_dict['cache_dir'],
        persist=True,
        cdf_kwargs=dict(chunks={'time': 4}), 
    )


def global_mean_timeseries_ann(query, name='global_mean_timeseries_ann',
                               center_time=True):
    """
    Instantiate a `funnel.Collection` object for computing 
    global mean, annual mean timeseries.
    """
    
    postproccess = [ops.global_mean, ops.resample_ann] 
    postproccess_kwargs = [dict(normalize=True, include_ms=False), {}]    
    
    if center_time:
        postproccess = [ops.center_time] + postproccess
        postproccess_kwargs = [{}] + postproccess_kwargs

    return fn.Collection(
        name=name,
        esm_collection_json=esm_collection_json,
        postproccess=postproccess,  
        postproccess_kwargs=postproccess_kwargs,
        query=query,
        cache_dir=config_dict['cache_dir'],
        persist=True,
        cdf_kwargs=dict(chunks={'time': 4}, decode_coords=False), 
    )    

def global_integral_timeseries_ann(query, name='global_integral_timeseries_ann',
                               center_time=True):
    """
    Instantiate a `funnel.Collection` object for computing 
    global mean, annual mean timeseries.
    """
    
    postproccess = [ops.global_mean, ops.resample_ann] 
    postproccess_kwargs = [dict(normalize=False, include_ms=False), {}]    
    
    if center_time:
        postproccess = [ops.center_time] + postproccess
        postproccess_kwargs = [{}] + postproccess_kwargs

    return fn.Collection(
        name=name,
        esm_collection_json=esm_collection_json,
        postproccess=postproccess,  
        postproccess_kwargs=postproccess_kwargs,
        query=query,
        cache_dir=config_dict['cache_dir'],
        persist=True,
        cdf_kwargs=dict(chunks={'time': 4}, decode_coords=False), 
    ) 
