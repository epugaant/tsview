import astropy.units as u
from astropy.units.core import UnitBase
from synphot import units
from astropy.time import Time

_u_vega = u.def_unit('VEGA')
VEGAMAG = u.mag(_u_vega)

PHOT_SYS = ['AB', 'ST', 'VEGA']
PHOT_SYS_MAG_UNIT = {'mag(AB)': u.ABmag, 'mag(ST)': u.STmag, 'mag(VEGA)': units.VEGAMAG}

SPECT_DENS = [
            'Jy', 'mJy', 'uJy',
            'W / (m2 Hz)', #'W / (Hz m2)',  # Order is different in astropy v5.3
            'eV / (s m2 Hz)', #'eV / (Hz s m2)',
            'erg / (s cm2)',
            'erg / (s cm2 Angstrom)', #'erg / (Angstrom s cm2)',
            'erg / (s cm2 Hz)', #'erg / (Hz s cm2)',
            'ph / (s cm2 Angstrom)', #'ph / (Angstrom s cm2)',
            'ph / (s cm2 Hz)', #'ph / (Hz s cm2)'
                ]

def equivalent_units(unit_str):
    '''Function to return spectral flux and spectral equivalences in BaseUnit and MagUnit as set of strings'''
    try:
        units = PHOT_SYS_MAG_UNIT[unit_str]
    except:
        units = u.Unit(unit_str)
        
    if units.physical_type in ['spectral flux density'] or (hasattr(units, 'physical_unit') and units.physical_unit != 'dimensionless'): # spectral flux
        eqv = u.spectral_density(1 * u.m)  # Value does not matter here.
        exclude_lower = {'flam', 'fnu', 'bol', 'photlam', 'photnu'}
        exclude_upper = {'FLAM', 'FNU', 'BOL', 'PHOTLAM', 'PHOTNU'}
        try:
            list_of_units = set(list(map(str, units.find_equivalent_units(
                include_prefix_units=False, equivalencies=eqv))) + PHOT_SYS + 
                                SPECT_DENS) - exclude_lower - exclude_upper
        except:
            list_of_units = set(PHOT_SYS + 
                                SPECT_DENS) - exclude_lower - exclude_upper
        
        # include magnitudes for three photometric systems
        list_of_units = list_of_units | {'mag({})'.format(val) for val in list_of_units if val in PHOT_SYS}
        # remove VEGA flux
        list_of_units.remove('VEGA')
        # remove self
        list_of_units.remove(unit_str)
    elif units.physical_type in ['length', 'frequency', 'wavenumber', 'energy']:  # spectral axis
        # prefer Hz over Bq and um over micron
        exclude = {'Bq', 'micron'}
        list_of_units = set(list(map(str, units.find_equivalent_units(
            include_prefix_units=False, equivalencies=u.spectral()))) + ['um']) - exclude
    else:
        list_of_units = []
    return list_of_units

def time_units(time):
    '''Function to return time equivalences in as set of strings'''
    exclude = {'cxcsec', 'datetime', 'gps', 'unix', 'unix_tai', 'ymdhms', 'datetime64'}
    list_of_units = set(list(map(str, time.FORMATS.keys()))) - exclude
    return list_of_units


l = equivalent_units('mJy') # 15 options
print(len(l))
l = equivalent_units('mag(AB)') # 15 options
l = equivalent_units('mag(VEGA)') # 15 options
l = equivalent_units('electron/s') #12 options (default)
l = equivalent_units('micron') #18 options
lt = time_units(Time.now())

#convert_unit uses 
unt = u.Unit(units.VEGAMAG)
if hasattr(unt, 'function_unit'):
    print(unt.function_unit, unt.physical_unit)
else:
    print(unt.physical_unit)
    
from astropy import units as u
from astropy.units import photometric
with u.add_enabled_units_context([photometric.AB, photometric.ST]):
    u.mJy.find_equivalent_units()  

u.mJy.find_equivalent_units()    
_u_vega = u.def_unit('VEGA')
VEGAMAG = u.mag(_u_vega)
_u_test = u.def_unit('test', 3*u.Jy)
u.add_enabled_units(_u_test)

# Register with astropy units
u.add_enabled_units([VEGAMAG])
u.mJy.find_equivalent_units()
list(map(str.lower,list(exclude)))





    