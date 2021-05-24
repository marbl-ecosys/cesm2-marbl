import yaml

import funnel as fn
import operators as ops


with open('config.yml') as fid:
    config_dict = yaml.load(fid, Loader=yaml.Loader)
    
        
def epoch_mean(epoch_def, query):
    """Instantiate a `funnel.Collection` object for computing epoch means."""
    
    time_slice_kwargs = None
    if isinstance(epoch_def, slice):
        time_slice_kwargs = dict(time=epoch_def)
        
    elif isinstance(epoch_def, str):
        if epoch_def in config_dict['epochs']:
            time_slice_kwargs = config_dict['epochs'][epoch_def]
 
    if time_slice_kwargs is None:
        raise ValueError(
            f'epoch_def: {epoch_def} not recognized as pre-defined or slice object'
        )
        
    return fn.Collection(
        esm_collection_json=config_dict['esm_collection'],
        postproccess=[ops.epoch_mean],  
        postproccess_kwargs=[dict(sel_dict=time_slice_kwargs)],
        query=query,
        cache_dir=config_dict['cache_dir'],
        persist=True,
        cdf_kwargs=dict(chunks={'time': 4}, decode_coords=False), 
    )    