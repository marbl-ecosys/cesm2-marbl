import yaml

import funnel as fn
import operators as ops
import variable_defs

with open('config.yml') as fid:
    config_dict = yaml.load(fid, Loader=yaml.Loader)
    
        
def epoch_mean(query, sel_dict, name='epoch_mean', center_time=True):
    """Instantiate a `funnel.Collection` object for computing epoch means."""
    postproccess = [ops.mean_time] 
    postproccess_kwargs = [dict(sel_dict=sel_dict)]    
    if center_time:
        postproccess = [ops.center_time] + postproccess
        postproccess_kwargs = [{}] + postproccess_kwargs

    return fn.Collection(
        name=name,
        esm_collection_json=config_dict['esm_collection'],
        postproccess=postproccess,  
        postproccess_kwargs=postproccess_kwargs,
        query=query,
        cache_dir=config_dict['cache_dir'],
        persist=True,
        cdf_kwargs=dict(chunks={'time': 4}, decode_coords=False), 
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
        esm_collection_json=config_dict['esm_collection'],
        postproccess=postproccess,  
        postproccess_kwargs=postproccess_kwargs,
        query=query,
        cache_dir=config_dict['cache_dir'],
        persist=True,
        cdf_kwargs=dict(chunks={'time': 4}, decode_coords=False), 
    )    

def global_integral_timeseries_ann(query, name='global_mean_timeseries_ann',
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
        esm_collection_json=config_dict['esm_collection'],
        postproccess=postproccess,  
        postproccess_kwargs=postproccess_kwargs,
        query=query,
        cache_dir=config_dict['cache_dir'],
        persist=True,
        cdf_kwargs=dict(chunks={'time': 4}, decode_coords=False), 
    ) 

def present_day_epoch_mean(experiment='historical', stream='pop.h'):
    """Instantiate a `funnel.Collection` object for computing epoch means."""
    query = dict(
        experiment=experiment,
        stream=stream,
    )
    
    sel_dict = dict(
        time=slice('1990', '2014')
    )       
    return epoch_mean(query, sel_dict, center_time=True)