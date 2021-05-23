import funnel as fn


# TODO: move this kind of thing to a yaml file with other units info
C_flux_vars = [
    'FG_CO2', 'photoC_TOT_zint_100m', 'photoC_diat_zint_100m', 'photoC_diaz_zint_100m',
    'photoC_sp_zint_100m', 'POC_FLUX_100m',
]


@fn.register_derived_var(
    varname='pCFC11',
    dependent_vars=['CFC11', 'TEMP', 'SALT'],
)
def derive_var_pCFC11(ds):
    """compute pCFC11"""
    from calc import calc_cfc11sol

    ds['pCFC11'] = ds['CFC11'] * 1e-9 / calc_cfc11sol(ds.SALT, ds.TEMP)
    ds.pCFC11.attrs['long_name'] = 'pCFC-11'
    ds.pCFC11.attrs['units'] = 'patm'
    if 'coordinates' in ds.TEMP.attrs:
        ds.pCFC11.attrs['coordinates'] = ds.TEMP.attrs['coordinates']
    ds.pCFC11.encoding = ds.TEMP.encoding
    return ds.drop(['CFC11', 'TEMP', 'SALT'])


@fn.register_derived_var(
    varname='pCFC12',
    dependent_vars=['CFC12', 'TEMP', 'SALT'],
)
def derive_var_pCFC12(ds):
    """compute pCFC12"""
    from calc import calc_cfc12sol

    ds['pCFC12'] = ds['CFC12'] * 1e-9 / calc_cfc12sol(ds['SALT'],ds['TEMP'])
    ds.pCFC12.attrs['long_name'] = 'pCFC-12'
    ds.pCFC12.attrs['units'] = 'patm'
    if 'coordinates' in ds.TEMP.attrs:
        ds.pCFC12.attrs['coordinates'] = ds.TEMP.attrs['coordinates']
    ds.pCFC12.encoding = ds.TEMP.encoding
    return ds.drop(['CFC12', 'TEMP', 'SALT'])      


@fn.register_derived_var(
    varname='Cant', 
    dependent_vars=['DIC', 'DIC_ALT_CO2'],
)
def derive_var_Cant(ds):
    """compute Cant"""
    ds['Cant'] = ds['DIC'] - ds['DIC_ALT_CO2']
    ds.Cant.attrs = ds.DIC.attrs
    ds.Cant.attrs['long_name'] = 'Anthropogenic CO$_2$'

    if 'coordinates' in ds.DIC.attrs:
        ds.Cant.attrs['coordinates'] = ds.DIC.attrs['coordinates']
    ds.Cant.encoding = ds.DIC.encoding
    return ds.drop(['DIC', 'DIC_ALT_CO2'])


@fn.register_derived_var(
    varname='Del14C', 
    dependent_vars=['ABIO_DIC14', 'ABIO_DIC'],
)
def derive_var_Del14C(ds):
    """compute Del14C"""
    ds['Del14C'] = 1000. * (ds['ABIO_DIC14'] / ds['ABIO_DIC'] - 1.)
    ds.Del14C.attrs = ds.ABIO_DIC14.attrs
    ds.Del14C.attrs['long_name'] = '$\Delta^{14}$C'
    ds.Del14C.attrs['units'] = 'permille'    

    if 'coordinates' in ds.ABIO_DIC14.attrs:
        ds.Del14C.attrs['coordinates'] = ds.ABIO_DIC14.attrs['coordinates']
    ds.Del14C.encoding = ds.ABIO_DIC14.encoding
    return ds.drop(['ABIO_DIC14', 'ABIO_DIC'])


@fn.register_derived_var(
    varname='SST', 
    dependent_vars=['TEMP'],
)
def derive_var_SST(ds):
    """compute SST"""
    ds['SST'] = ds['TEMP'].isel(z_t=0, drop=True)
    ds.SST.attrs = ds.TEMP.attrs
    ds.SST.attrs['long_name'] = 'SST'
    ds.SST.encoding = ds.TEMP.encoding
    if 'coordinates' in ds.TEMP.attrs:
        ds.SST.attrs['coordinates'] = ds.TEMP.attrs['coordinates'].replace('z_t', '')   
    return ds.drop('TEMP')   


