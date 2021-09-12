from functools import partial
import numpy as np
import xarray as xr

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


@fn.register_derived_var(
    varname='DOC_FLUX_IN_100m', 
    dependent_vars=[
        'DIA_IMPVF_DOC',
        'HDIFB_DOC',
        'WT_DOC',
    ],
)
def derive_var_DOC_FLUX_IN_100m(ds):
    """compute DOC flux across 100m (positive down)"""
    k_100m_top = np.where(ds.z_w_top == 100e2)[0][0]
    k_100m_bot = np.where(ds.z_w_bot == 100e2)[0][0]    
    DIA_IMPVF = ds.DIA_IMPVF_DOC.isel(z_w_bot=k_100m_bot)
    HDIFB = ds.HDIFB_DOC.isel(z_w_bot=k_100m_bot) * ds.dz[k_100m_bot]
    WT = (-1.0) * ds.WT_DOC.isel(z_w_top=k_100m_top) * ds.dz[k_100m_top]
    
    ds['DOC_FLUX_IN_100m'] = (DIA_IMPVF + HDIFB + WT)
    ds.DOC_FLUX_IN_100m.attrs = ds.WT_DOC.attrs
    ds.DOC_FLUX_IN_100m.attrs['long_name'] = 'DOC flux across 100 m (positive down)'
    ds.DOC_FLUX_IN_100m.attrs['units'] = 'nmol/s/cm^2'
    ds.DOC_FLUX_IN_100m.encoding = ds.WT_DOC.encoding
    return ds.drop(['DIA_IMPVF_DOC', 'HDIFB_DOC', 'WT_DOC'])   


@fn.register_derived_var(
    varname='DOCr_FLUX_IN_100m', 
    dependent_vars=[
        'DIA_IMPVF_DOCr',
        'HDIFB_DOCr',
        'WT_DOCr',
    ],
)
def derive_var_DOCr_FLUX_IN_100m(ds):
    """compute DOCr flux across 100m (positive down)"""
    k_100m_top = np.where(ds.z_w_top == 100e2)[0][0]
    k_100m_bot = np.where(ds.z_w_bot == 100e2)[0][0]    
    DIA_IMPVF = ds.DIA_IMPVF_DOCr.isel(z_w_bot=k_100m_bot)
    HDIFB = ds.HDIFB_DOCr.isel(z_w_bot=k_100m_bot) * ds.dz[k_100m_bot]
    WT = (-1.0) * ds.WT_DOCr.isel(z_w_top=k_100m_top) * ds.dz[k_100m_top]
    
    ds['DOCr_FLUX_IN_100m'] = (DIA_IMPVF + HDIFB + WT)
    ds.DOCr_FLUX_IN_100m.attrs = ds.WT_DOCr.attrs
    ds.DOCr_FLUX_IN_100m.attrs['long_name'] = 'DOCr flux across 100 m (positive down)'
    ds.DOCr_FLUX_IN_100m.attrs['units'] = 'nmol/s/cm^2'
    ds.DOCr_FLUX_IN_100m.encoding = ds.WT_DOCr.encoding
    return ds.drop(['DIA_IMPVF_DOCr', 'HDIFB_DOCr', 'WT_DOCr'])  


@fn.register_derived_var(
    varname='DOCt_FLUX_IN_100m', 
    dependent_vars=[
        'DIA_IMPVF_DOC',
        'HDIFB_DOC',
        'WT_DOC',        
        'DIA_IMPVF_DOCr',
        'HDIFB_DOCr',
        'WT_DOCr',
    ],
)
def derive_var_DOCt_FLUX_IN_100m(ds):
    """compute DOCt (DOC + DOCr) flux across 100m (positive down)"""
    k_100m_top = np.where(ds.z_w_top == 100e2)[0][0]
    k_100m_bot = np.where(ds.z_w_bot == 100e2)[0][0]    

    DIA_IMPVF = ds.DIA_IMPVF_DOC.isel(z_w_bot=k_100m_bot)
    HDIFB = ds.HDIFB_DOC.isel(z_w_bot=k_100m_bot) * ds.dz[k_100m_bot]
    WT = (-1.0) * ds.WT_DOC.isel(z_w_top=k_100m_top) * ds.dz[k_100m_top]
    
    DIA_IMPVF += ds.DIA_IMPVF_DOCr.isel(z_w_bot=k_100m_bot)
    HDIFB += ds.HDIFB_DOCr.isel(z_w_bot=k_100m_bot) * ds.dz[k_100m_bot]
    WT += (-1.0) * ds.WT_DOCr.isel(z_w_top=k_100m_top) * ds.dz[k_100m_top]
    
    ds['DOCt_FLUX_IN_100m'] = (DIA_IMPVF + HDIFB + WT)
    ds.DOCt_FLUX_IN_100m.attrs = ds.WT_DOC.attrs
    ds.DOCt_FLUX_IN_100m.attrs['long_name'] = 'Total DOC flux across 100 m (positive down)'
    ds.DOCt_FLUX_IN_100m.attrs['units'] = 'nmol/s/cm^2'
    ds.DOCt_FLUX_IN_100m.encoding = ds.WT_DOC.encoding   
    
    return ds.drop([
        'DIA_IMPVF_DOC', 'HDIFB_DOC', 'WT_DOC',
        'DIA_IMPVF_DOCr', 'HDIFB_DOCr', 'WT_DOCr'
    ])


@fn.register_derived_var(
    varname='DOCt', 
    dependent_vars=['DOC', 'DOCr'],
)
def derive_var_DOCt(ds):
    """compute DOCt"""
    ds['DOCt'] = ds['DOC'] + ds['DOCr']
    ds.DOCt.attrs = ds.DOC.attrs
    ds.DOCt.attrs['long_name'] = 'Dissolved Organic Carbon (total)'
    ds.DOCt.encoding = ds.DOC.encoding
    return ds.drop(['DOC', 'DOCr'])


@fn.register_derived_var(
    varname='DONt', 
    dependent_vars=['DON', 'DONr'],
)
def derive_var_DONt(ds):
    """compute DONt"""
    ds['DONt'] = ds['DON'] + ds['DONr']
    ds.DONt.attrs = ds.DON.attrs
    ds.DONt.attrs['long_name'] = 'Dissolved Organic Nitrogen (total)'
    ds.DONt.encoding = ds.DON.encoding
    return ds.drop(['DON', 'DONr'])


@fn.register_derived_var(
    varname='DOPt', 
    dependent_vars=['DOP', 'DOPr'],
)
def derive_var_DOPt(ds):
    """compute DOPt"""
    ds['DOPt'] = ds['DOP'] + ds['DOPr']
    ds.DOPt.attrs = ds.DOP.attrs
    ds.DOPt.attrs['long_name'] = 'Dissolved Organic Phosphorus (total)'
    ds.DOPt.encoding = ds.DOP.encoding
    return ds.drop(['DOP', 'DOPr'])


@fn.register_derived_var(
    varname='Omega_calc', 
    dependent_vars=['CO3', 'co3_sat_calc'],
)
def derive_var_Omega_calc(ds):
    """compute Omega calcite"""
    ds['Omega_calc'] = ds['CO3'] / ds['co3_sat_calc']
    ds.Omega_calc.attrs = ds.CO3.attrs
    ds.Omega_calc.attrs['long_name'] = '$\Omega_{calc}$'
    ds.Omega_calc.attrs['units'] = ''    
    ds.Omega_calc.encoding = ds.CO3.encoding
    return ds.drop(['CO3', 'co3_sat_calc'])


@fn.register_derived_var(
    varname='Omega_arag', 
    dependent_vars=['CO3', 'co3_sat_arag'],
)
def derive_var_Omega_arag(ds):
    """compute Omega aragonite"""
    ds['Omega_arag'] = ds['CO3'] / ds['co3_sat_arag']
    ds.Omega_arag.attrs = ds.CO3.attrs
    ds.Omega_arag.attrs['long_name'] = '$\Omega_{arag}$'
    ds.Omega_arag.attrs['units'] = ''        
    ds.Omega_arag.encoding = ds.CO3.encoding    
    return ds.drop(['CO3', 'co3_sat_arag'])



def snormalize(X, S, Sbar=35.):
    """compute salinity normalized tracer values"""
    return X * Sbar / S


@fn.register_derived_var(
    varname='sDIC', 
    dependent_vars=['DIC', 'SALT'],
)
def derive_var_sDIC(ds):
    """compute salinity normalized DIC"""
    ds['sDIC'] = snormalize(ds.DIC, ds.SALT, 35.)
    ds.sDIC.attrs = ds.DIC.attrs
    ds.sDIC.attrs['long_name'] = 'DIC (salinity normalized)'
    ds.sDIC.encoding = ds.DIC.encoding
    return ds.drop(['DIC', 'SALT'])


@fn.register_derived_var(
    varname='sALK', 
    dependent_vars=['ALK', 'SALT'],
)
def derive_var_sALK(ds):
    """compute salinity normalized ALK"""
    ds['sALK'] = snormalize(ds.ALK, ds.SALT, 35.)
    ds.sALK.attrs = ds.ALK.attrs
    ds.sALK.attrs['long_name'] = 'ALK (salinity normalized)'
    ds.sALK.encoding = ds.ALK.encoding
    return ds.drop(['ALK', 'SALT'])


@fn.register_derived_var(
    varname='zoo_prod_zint_100m', 
    dependent_vars=['graze_diat_zoo_zint_100m',
                    'graze_diaz_zoo_zint_100m',
                    'graze_sp_zoo_zint_100m',],
)
def derive_var_zoo_prod_zint_100m(ds):
    ds['zoo_prod_zint_100m'] = (
        ds.graze_diat_zoo_zint_100m + 
        ds.graze_diaz_zoo_zint_100m +
        ds.graze_sp_zoo_zint_100m 
    )
    ds.zoo_prod_zint_100m.attrs = ds.graze_diat_zoo_zint_100m.attrs
    ds.zoo_prod_zint_100m.attrs['long_name'] = 'Zooplankton production'
    ds.zoo_prod_zint_100m.encoding = ds.graze_diat_zoo_zint_100m.encoding
    return ds.drop([
        'graze_diat_zoo_zint_100m',
        'graze_diaz_zoo_zint_100m',
        'graze_sp_zoo_zint_100m',
    ])
        

def compute_zint(ds, varname, depth_slice=None):
    """compute vertical integral over `depth_slice`"""
    da = ds[varname]
    with xr.set_options(keep_attrs=True):
        if 'z_t' in da.dims:            
            dz = ds.dz
            if depth_slice is not None:
                dao = (dz * da).sel(z_t=depth_slice).sum('z_t')
            else:
                dao = (dz * da).sum('z_t')
        
        elif 'z_t_150m' in da.dims:
            dz = ds.dz.isel(z_t=slice(0, 15)).rename({'z_t': 'z_t_150m'})            
            if depth_slice is not None:            
                dao = (dz * da).sel(z_t_150m=depth_slice).sum('z_t_150m')
            else:
                dao = (dz * da).sum('z_t_150m')                
    dao.attrs['units'] = da.attrs['units'] + ' cm'
    ds[f'{varname}_zint'] = dao
    return ds.drop([varname])
    
    
def compute_zint_100m(ds, varname):
    """compute integral over top 100 m"""
    return compute_zint(
        ds, varname, depth_slice=slice(0, 100e2)
    ).rename({f'{varname}_zint': f'{varname}_zint_100m'})


# register vertical integrals for PFTs
for v in ['sp', 'diat', 'diaz', 'zoo']:
    func = partial(compute_zint_100m, varname=f'{v}C')
    fn.register_derived_var(func, varname=f'{v}C_zint_100m', dependent_vars=[f'{v}C'])
