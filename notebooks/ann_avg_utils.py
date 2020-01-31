def get_pint_units():
    from pint import UnitRegistry

    # Add new units to UnitRegistry
    units = UnitRegistry()
    units.define('gram N = mol / 14 = gN')
    units.define('gram C = mol / 12 = gC')
    units.define('year = 365 day = yr')

    # Define final units
    PgC_per_year = 'PgC/yr'
    TgN_per_year = 'TgN/yr'
    uM = 'uM'

    final_units = dict()
    final_units['photoC_TOT_zint'] = PgC_per_year
    final_units['photoC_diat_zint'] = PgC_per_year
    final_units['photoC_TOT_zint_100m'] = PgC_per_year
    final_units['photoC_diat_zint_100m'] = PgC_per_year
    final_units['POC_FLUX_100m'] = PgC_per_year
    final_units['CaCO3_FLUX_100m'] = PgC_per_year
    final_units['diaz_Nfix'] = TgN_per_year
    final_units['NOx_FLUX'] = TgN_per_year
    final_units['NHy_FLUX'] = TgN_per_year
    final_units['DENITRIF'] = TgN_per_year
    final_units['SedDenitrif'] = TgN_per_year
    final_units['DON_RIV_FLUX'] = TgN_per_year
    final_units['DONr_RIV_FLUX'] = TgN_per_year
    final_units['FG_CO2'] = PgC_per_year
    final_units['O2'] = 'uM'
    final_units['O2_under_thres'] = 'Pm * m^2'

    return units, final_units

def get_ann_means_and_units(xp_dir, vars, experiments, experiment_longnames, units):
    import os
    import xarray as xr

    cache_dir = os.path.join(os.path.sep, 'glade', 'p', 'cgd', 'oce', 'projects', 'cesm2-marbl', 'xpersist_cache', xp_dir)

    ann_avg = dict()
    cesm_units = dict()
    for variable in vars:
        ann_avg[variable] = dict()
        cesm_units[variable] = dict()
        for model_version in experiments:
            for exp in experiments[model_version]:
                filename = os.path.join(cache_dir, f'{exp}_{variable}.nc')
                # Skip files that do not exist
                if not os.path.isfile(filename):
                    continue
                ann_avg[variable][exp] = xr.open_dataset(filename)
                cesm_units[variable][exp] = units[ann_avg[variable][exp][variable].attrs['units']]

    # Build time series that combines cesm1_hist to 2004 and cesm1_RCP85 after
    if 'cesm1_hist' in experiments['cesm1'] and 'cesm1_RCP85' in experiments['cesm1']:
        for var in ann_avg:
            if 'cesm1_hist' in ann_avg[var]:
                ann_avg[var]['cesm1_hist_RCP85'] = xr.concat([ann_avg[var]['cesm1_hist'].isel(time=slice(0,-1)),
                                                              ann_avg[var]['cesm1_RCP85']],
                                                             dim='time'
                                                            )
                cesm_units[var]['cesm1_hist_RCP85'] = cesm_units[var]['cesm1_hist'].copy()

    return ann_avg, cesm_units

def print_exp_time_bounds(ann_avg_var0, time_slices):
    # Verify time bounds for each experiment
    for exp in ann_avg_var0:
        try:
            bounds = list(ann_avg_var0[exp].sel(time=time_slices[exp]).time_bound.values[ind] for ind in [(0,0), (-1,-1)])
        except:
            bounds = list(ann_avg_var0[exp].isel(time=time_slices[exp]).time_bound.values[ind] for ind in [(0,0), (-1,-1)])
        print(f'Experiment: {exp}\nRequested time bounds\n----\n{bounds}\n\n')

def global_vars():
    # all variables will be packed into a dictionary
    all_vars = dict()

    # Should the annual averages include the marginal seas or just be open ocean?
    # (USE TRUE FOR PAPER)
    all_vars['include_marg_seas'] = True
    # all_vars['include_marg_seas'] = False

    # Subdirectory in xpersist cache containing netCDF files
    all_vars['xp_dir'] = 'with_marginal_seas' if all_vars['include_marg_seas'] else 'no_marginal_seas'

    # Process for updating intake-esm catalog
    #       1. download all data from HPSS via get_ocn_cmip5_files.sh
    #       2. rm /glade/u/home/mlevy/.intake_esm/collections/CESM1-CMIP5.nc
    #       3. regenerate it via Anderson's legacy intake-esm
    #       4. re-run build intake collections notebook
    #       5. commit change to .csv.gz in /glade/work/mlevy/intake-esm-collection/csv.gz/
    # NOTE: steps 2-5 can be done with notebooks/intake-esm-collection-defs/rebuild.sh
    all_vars['vars'] = [
                        'photoC_TOT_zint_100m', 'photoC_diat_zint_100m',
                        'photoC_TOT_zint', 'photoC_diat_zint',
                        'POC_FLUX_100m', 'CaCO3_FLUX_100m',
                        'diaz_Nfix', 'NOx_FLUX', 'NHy_FLUX', 'DENITRIF',
                        'SedDenitrif', 'DON_RIV_FLUX', 'DONr_RIV_FLUX',
                        'FG_CO2', 'O2' ,
                        'O2_under_thres' # add a thres dimension corresponding to limits
                       ]

    # experiments is a list of experiments to compute values for
    all_vars['experiments'] = dict()
    all_vars['experiments']['cesm1'] = [
                                        'cesm1_PI',
                                        'cesm1_PI_esm',
                                        'cesm1_hist',
                                        'cesm1_hist_esm',
                                        'cesm1_RCP85',
                                       ]
                   # CESM 2
    all_vars['experiments']['cesm2'] = [
                                        'cesm2_PI',
                                        'cesm2_hist',
                                        'cesm2_SSP1-2.6',
                                        'cesm2_SSP2-4.5',
                                        'cesm2_SSP3-7.0',
                                        'cesm2_SSP5-8.5',
                                       ]

    # experiment_longnames defines the table headers
    all_vars['experiment_longnames']={
                                      'cesm1_PI' : 'preindustrial (CESM1)',
                                      'cesm1_PI_esm' : 'preindustrial (CESM1, BPRP)',
                                      'cesm1_hist' : '1981-2005 (CESM1)',
                                      'cesm1_hist_esm' : '1990s (CESM1)',
                                      'cesm1_hist_RCP85' : '1990 - 2014 (CESM1)',
                                      'cesm1_RCP45' : 'RCP 4.5 2090s (CESM1)', # not available yet
                                      'cesm1_RCP85' : 'RCP 8.5 2090s (CESM1)',
                                      'cesm1_RCP85_esm' : 'RCP 8.5 2090s (CESM1)',
                                      'cesm2_PI' : 'preindustrial (CESM2)',
                                      'cesm2_hist' : '1990-2014 (CESM2)',
                                      'cesm2_SSP1-2.6' : 'RCP26 2090s (CESM2)',
                                      'cesm2_SSP2-4.5' : 'RCP45 2090s (CESM2)',
                                      'cesm2_SSP3-7.0' : 'RCP70 2090s (CESM2)',
                                      'cesm2_SSP5-8.5' : 'RCP85 2090s (CESM2)'
                                     }

    # experiment_dict determines which module version & intake data each experiment uses
    all_vars['experiment_dict'] = {
                                   'cesm1_PI' : ('cesm1', 'piControl'),
                                   'cesm1_PI_esm' : ('cesm1', 'esm-piControl'),
                                   'cesm1_hist' : ('cesm1', 'historical'),
                                   'cesm1_hist_esm' : ('cesm1', 'esm-hist'),
                                   'cesm1_RCP85' : ('cesm1', 'RCP-8.5'),
                                   'cesm1_RCP85_esm' : ('cesm1', 'esm-RCP-8.5'),
                                   'cesm2_PI' : ('cesm2', 'piControl'),
                                   'cesm2_hist' : ('cesm2', 'historical'),
                                   'cesm2_SSP1-2.6' : ('cesm2', 'SSP1-2.6'),
                                   'cesm2_SSP2-4.5' : ('cesm2', 'SSP2-4.5'),
                                   'cesm2_SSP3-7.0' : ('cesm2', 'SSP3-7.0'),
                                   'cesm2_SSP5-8.5' : ('cesm2', 'SSP5-8.5')
                                  }

    # NOTE: 2090-01-01 0:00:00 is the time stamp on the Dec 2089 monthly average
    #       So slice("2090", "2100") would actually return Dec 2090 - Nov 2099
    #       Specifying a day mid-month gets us to Jan 2090 - Dec 2099 (the 2090s)
    #       (this can be verified by looking at time bounds)
    time_slices_1990s = slice("1990-01-15", "2000-01-15")
    time_slices_2090s = slice("2090-01-15", "2100-01-15")
    time_slices_CMIP5 = slice("1981-01-15", "2006-01-15")
    time_slices_CMIP6 = slice("1990-01-15", "2015-01-15")

    all_vars['time_slices'] = dict()

    # 200 year averages for CESM1 PI runs, per Lindsay et al 2014
    # (He starts 30 years prior to branch point, so I will too)
    all_vars['time_slices']['cesm1_PI'] = slice(120, 320) # cfunits doesn't do years too far in past; this is 121-07-01 - 320-07-01
    all_vars['time_slices']['cesm1_PI_esm'] = slice(320, 520) # cfunits doesn't do years too far in past; this is 321-07-01 - 520-07-01
    # For CESM2, going from 50 years prior to first historical branch point
    #                  to 50 years after end of last historical member
    # TODO: These dates should be computed automatically based on intake metadata!
    all_vars['time_slices']['cesm2_PI'] = slice(550, 1070) # cfunits doesn't do years too far in past; this is 551-07-01 - 1070-07-01

    # Historical runs all use slightly different time periods
    # Note: that the annual mean data is actually running from July 1st to June 30th
    #       these slices were defined to work with monthly data, but pick up the correct years as well
    all_vars['time_slices']['cesm1_hist'] = time_slices_CMIP5 # per Lindsay et al 2014
    all_vars['time_slices']['cesm1_hist_esm'] = time_slices_1990s # per Moore et al 2013
    all_vars['time_slices']['cesm1_hist_RCP85'] = time_slices_CMIP6 # For our paper
    all_vars['time_slices']['cesm2_hist'] = time_slices_CMIP6 # For our paper

    # RCP runs use 2090s
    all_vars['time_slices']['cesm1_RCP45'] = time_slices_2090s
    all_vars['time_slices']['cesm1_RCP85'] = time_slices_2090s
    all_vars['time_slices']['cesm1_RCP85_esm'] = time_slices_2090s
    all_vars['time_slices']['cesm2_SSP1-2.6'] = time_slices_2090s
    all_vars['time_slices']['cesm2_SSP2-4.5'] = time_slices_2090s
    all_vars['time_slices']['cesm2_SSP3-7.0'] = time_slices_2090s
    all_vars['time_slices']['cesm2_SSP5-8.5'] = time_slices_2090s

    return all_vars

def get_table_specs(final_units, o2_levs=[]):
    table_specs = {
                  'POC' : {
                           'key' : 'Sinking POC at 100 m',
                           'units' : final_units['POC_FLUX_100m'],
                           'rounding' : 2
                          },
                  'CaCO3' : {
                             'key' : 'Sinking CaCO$_3$ at 100 m',
                             'units' : final_units['CaCO3_FLUX_100m'],
                             'rounding' : 3
                            },
                  'rain' : {
                            'key' : 'Rain ratio (CaCO$_3$/POC) at 100 m',
                            'units' : None,
                            'rounding' : 3
                           },
                  'NPP' : {
                           'key' : 'Net primary production, full depth',
                           'units' : final_units['photoC_TOT_zint'],
                           'rounding' : 1
                          },
                  'NPP_diat' : {
                                'key' : 'Diatom primary production, full depth',
                                'units' : '\%',
                                'rounding' : 0
                               },
                  'NPP_100m' : {
                                'key' : 'Net primary production, top 100m',
                                'units' : final_units['photoC_TOT_zint_100m'],
                                'rounding' : 1
                               },
                  'NPP_diat_100m' : {
                                     'key' : 'Diatom primary production, top 100m',
                                     'units' : '\%',
                                     'rounding' : 0
                                    },
                  'Nfix' : {
                            'key' : 'Nitrogen fixation',
                            'units' : final_units['diaz_Nfix'],
                            'rounding' : 0
                           },
                  'Ndep' : {
                            'key' : 'Nitrogen deposition',
                            'units' : final_units['NOx_FLUX'],
                            'rounding' : 1
                           },
                  'denitrif' : {
                                'key' : 'Water Column Denitrification',
                                'units' : final_units['DENITRIF'],
                                'rounding' : 0
                               },
                  'denitrif2' : {
                                 'key' : 'Sediment Denitrification',
                                 'units' : final_units['SedDenitrif'],
                                 'rounding' : 0
                                },
                  'rivflux' : {
                               'key' : 'Nitrogen River Flux',
                               'units' : final_units['DON_RIV_FLUX'],
                               'rounding' : 0
                              },
                  'Ncycle' : {
                              'key' : 'N cycle imbalance',
                              'units' : final_units['diaz_Nfix'],
                              'rounding' : 0
                             },
                  'CO2' : {
                           'key' : 'Air–sea CO2 flux',
                           'units' : final_units['FG_CO2'],
                           'rounding' : 2
                          },
                  'O2' : {
                          'key' : 'Mean ocean oxygen',
                          'rounding' : 0
                         }
                 }

    if final_units['O2'] == 'uM':
        table_specs['O2']['units'] = '$\mu$M'
    else:
        table_specs['O2']['units'] = final_units['O2']
    
    for o2_thres in o2_levs:
        table_specs[f'o2_under_{o2_thres}'] = {'key' : _O2_vol_keys(o2_thres), 'rounding' : 0}
        if final_units['O2_under_thres'] == 'Pm * m^2':
            table_specs[f'o2_under_{o2_thres}']['units'] = '10$^1$$^5$ m$^3$'
        else:
            table_specs[f'o2_under_{o2_thres}']['units'] = final_units['O2_under_thres']

    return table_specs

def compute_diagnostic_values(experiments, table_specs, ann_avg, time_slices, cesm_units, final_units, verbose=False):
    diagnostic_values = dict()
    kwargs = dict()
    kwargs['ann_avg'] = ann_avg
    kwargs['time_slices'] = time_slices
    kwargs['cesm_units'] = cesm_units
    kwargs['final_units'] = final_units
    for model_version in experiments:
        exp_loop = experiments[model_version]
        if model_version == 'cesm1' and 'cesm1_hist' in exp_loop and 'cesm1_RCP85' in exp_loop:
            exp_loop.append('cesm1_hist_RCP85')
        for exp in exp_loop:
            kwargs['exp'] = exp
            diagnostic_values[exp] = dict()
            # Compute each value by hand
            if verbose:
                print(f'Computing 100m POC flux for {exp}')
            diagnostic_values[exp][table_specs['POC']['key']] = _get_time_and_ensemble_mean('POC_FLUX_100m', **kwargs)

            if verbose:
                print(f'Computing 100m CaCO3 flux for {exp}')
            diagnostic_values[exp][table_specs['CaCO3']['key']] = _get_time_and_ensemble_mean('CaCO3_FLUX_100m', **kwargs)

            if verbose:
                print(f'Computing 100m rain rate for {exp}')
            try:
                diagnostic_values[exp][table_specs['rain']['key']] = (diagnostic_values[exp][table_specs['CaCO3']['key']] /
                                                     diagnostic_values[exp][table_specs['POC']['key']])
            except:
                print(f'   * Can not compute rain rate for {exp}')

            if verbose:
                print(f'Computing full depth net primary production for {exp}')
            diagnostic_values[exp][table_specs['NPP']['key']] = _get_time_and_ensemble_mean('photoC_TOT_zint', **kwargs)

            if verbose:
                print(f'Computing full depth primary production from diatoms for {exp}')
            try:
                diagnostic_values[exp][table_specs['NPP_diat']['key']] = 100*(_get_time_and_ensemble_mean('photoC_diat_zint', **kwargs) /
                                                            diagnostic_values[exp][table_specs['NPP']['key']])
            except:
                print(f'   * Can not compute primary production from diatoms for {exp}')

            if verbose:
                print(f'Computing top 100m net primary production for {exp}')
            diagnostic_values[exp][table_specs['NPP_100m']['key']] = _get_time_and_ensemble_mean('photoC_TOT_zint_100m', **kwargs)

            if verbose:
                print(f'Computing top 100m primary production from diatoms for {exp}')
            try:
                diagnostic_values[exp][table_specs['NPP_diat_100m']['key']] = 100*(_get_time_and_ensemble_mean('photoC_diat_zint_100m', **kwargs) /
                                                            diagnostic_values[exp][table_specs['NPP_100m']['key']])
            except:
                print(f'   * Can not compute primary production from diatoms for {exp}')

            if verbose:
                print(f'Computing Nfixation for {exp}')
            diagnostic_values[exp][table_specs['Nfix']['key']] = _get_time_and_ensemble_mean('diaz_Nfix', **kwargs)

            if verbose:
                print(f'Computing Ndep for {exp}')
            diagnostic_values[exp][table_specs['Ndep']['key']] = (_get_time_and_ensemble_mean('NOx_FLUX', **kwargs) +
                                                _get_time_and_ensemble_mean('NHy_FLUX', **kwargs))

            if verbose:
                print(f'Computing Water Column Denitrif for {exp}')
            diagnostic_values[exp][table_specs['denitrif']['key']] = _get_time_and_ensemble_mean('DENITRIF', **kwargs)

            if verbose:
                print(f'Computing Sediment Denitrif for {exp}')
            diagnostic_values[exp][table_specs['denitrif2']['key']] = _get_time_and_ensemble_mean('SedDenitrif', **kwargs)

            if verbose:
                print(f'Computing Nitrogen River Flux for {exp}')
            diagnostic_values[exp][table_specs['rivflux']['key']] = (_get_time_and_ensemble_mean('DON_RIV_FLUX', **kwargs) +
                                                   _get_time_and_ensemble_mean('DONr_RIV_FLUX', **kwargs))

            if verbose:
                print(f'Computing Nitrogen Cycle imbalance for {exp}')
            try:
                diagnostic_values[exp][table_specs['Ncycle']['key']] = (diagnostic_values[exp][table_specs['Ndep']['key']] +
                                                                        diagnostic_values[exp][table_specs['Nfix']['key']] -
                                                                        diagnostic_values[exp][table_specs['denitrif']['key']]
                                                                       )
                try:
                    diagnostic_values[exp][table_specs['Ncycle']['key']] = (diagnostic_values[exp][table_specs['Ncycle']['key']] -
                                                                            diagnostic_values[exp][table_specs['denitrif2']['key']] - 
                                                                            diagnostic_values[exp][table_specs['rivflux']['key']]
                                                                           )
                except:
                    print(f'   * No additional denitrification terms for {exp}')
                    pass
            except:
                print(f'   * Can not compute Ncycle imbalance for {exp}')

            if verbose:
                print(f'Computing air-sea CO2 Flux for {exp}')
            diagnostic_values[exp][table_specs['CO2']['key']] = _get_time_and_ensemble_mean('FG_CO2', **kwargs)

            # Update O2 units to account for fact that we are dividing my total volume
            if verbose:
                print(f'Computing O2 concentration for {exp}')
            try:
                diagnostic_values[exp][table_specs['O2']['key']] = _get_time_and_ensemble_mean('O2', **kwargs)
            except:
                print(f'   * Can not compute O2 concentration for {exp}')

            try:
                if exp in ann_avg['O2_under_thres']:
                    for n, o2_thres in enumerate(ann_avg['O2_under_thres'][exp]['o2_thres'].data):
                        if verbose:
                            print(f'Computing volume where O2 < {o2_thres} uM for {exp}')
                        diagnostic_values[exp][table_specs[f'o2_under_{o2_thres}']['key']] = _get_time_and_ensemble_mean('O2_under_thres', **kwargs)[n]
            except:
                print(f'   * Can not compute O2 volumes under thresholds for {exp}')

            if verbose and (model_version != list(experiments.keys())[-1] or exp != experiments[model_version][-1]):
                print('\n----\n')

    return diagnostic_values

def _get_time_and_ensemble_mean(variable, ann_avg, exp, time_slices, cesm_units, final_units):
    try:
        if exp in ['cesm1_PI', 'cesm1_PI_esm', 'cesm2_PI']:
            # Need isel instead of sel since PI slices are in index space rather than years
            ens_time_mean = (ann_avg[variable][exp][variable].isel(time=time_slices[exp]).mean('member_id')).mean('time').values
        else:
            ens_time_mean = (ann_avg[variable][exp][variable].sel(time=time_slices[exp]).mean('member_id')).mean('time').values
    except:
        print(f'   * Can not compute {variable} for {exp}')
        return('-')
    return((ens_time_mean * cesm_units[variable][exp]).to(final_units[variable]))

# Define keys that will go into table columns
def _O2_vol_keys(o2_thres):
    if o2_thres == 20:
        return f'OMZ volume (O$_2$ $<$20 $\mu$M)'
    return f'Volume where O$_2$ $<${o2_thres} $\mu$M'