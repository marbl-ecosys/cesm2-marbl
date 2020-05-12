
cdf_kwargs = dict(
    chunks=dict(time=180),
    decode_coords=False,
    decode_times=False,
)

unit_str_replacements = {
    'mmol/m^3 cm/s': 'nmol/cm^2/s',
    'centimeters': 'cm',
    'centimeter': 'cm',    
}
