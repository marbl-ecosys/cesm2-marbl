import yaml
import xarray as xr

from . config import units_info_yaml

@xr.register_dataset_accessor('convert_units')
class convert_units(object):
    
    # class attributes
    from pint import UnitRegistry
    units_reg = UnitRegistry()
    
    def __init__(self, ds=None):
        """initialize class"""
        self._ds = ds

        # TODO: provide user control on this file as settable parameter
        # define attrs from yaml
        with open(units_info_yaml, 'r') as fid:
            units_info = yaml.safe_load(fid) 

        self.var_groups = units_info['variable_groups']            
        self.latex_replacements = units_info['latex_replacements']
        self.expected_units = units_info['expected_units']
            
        # Add new units to UnitRegistry
        for pint_def in units_info['pint_definitions']:
            self.units_reg.define(pint_def)
            
    def _get_pint_units_out(self, varname, units_in):
        """determine output units on the basis of varname and inputs"""

        # TODO: change all print statements to warnings
        #       or better yet would be a logging function that enables failures
        #       so we can test each condition
        
        # determine the "type" of variable looping over defined groups
        for group_name, group_var_list in self.var_groups.items():
            
            # if the variable is in group "group_name"
            if varname in group_var_list:
                
                # do we have a rubric for converting units for this group?
                if group_name in self.expected_units:                    
                    units_from_to_dict = self.expected_units[group_name]                    
                    
                    # are the units in the var attrs what we expect?
                    if units_in in units_from_to_dict:
                        return units_from_to_dict[units_in]
                    else:
                        print(f'warning: unknown units: {units_in} for {varname}')
                        return
                else:
                    print(f'warning: unit conversions not specified {group_name}')
                    return 
        
        print(f'var type not specified for {varname}\nadd to variables.yml')

    def _conversion_factor(self, units_in, units_out):
        """return a scalar conversion factor"""
        return (1. * self.units_reg[units_in]).to(units_out).magnitude
    
        
    def __call__(self, varname, inplace=False, latex=True):        
        """apply units conversion to Dataset"""
        
        if inplace:
            ds = self._ds
        else:
            ds = self._ds.copy()
        
        if isinstance(varname, str):
            varname = [varname]
        else:
            assert isinstance(varname, list), (
                'varname must be a list'
            )
        
        # loop over all variables and convert units
        for v in varname:
            if 'units' not in ds[v].attrs:
                print(f'{v} has no units attribute, skipping')
                continue
            
            # get input units and determine output
            units_in = ds[v].attrs['units']
            units_out = self._get_pint_units_out(v, units_in)            
            if units_out is None:
                print(f'{v}: could not determine target units, skipping')
                continue
            
            # get conversion factor and apply it
            factor = self._conversion_factor(units_in, units_out)        
            ds[v].data = ds[v].data * factor

            # set units string
            if latex: 
                units_out = self._latex(units_out)
            ds[v].attrs['units'] = units_out   

        return ds
    
    def _latex(self, units_out):
        """make latex strings out of pint-compatible units"""
        
        for old, new in self.latex_replacements.items():
            if old in units_out:
                units_out = units_out.replace(old, new)
        return units_out