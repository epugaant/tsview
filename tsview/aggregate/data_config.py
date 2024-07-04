
import astropy.units as u
from synphot import units

DATA_DICT = {
'gaia': {
    'system': 'VEGAMAG',
    'zeropt': {'G': {'VEGAMAG': {'zp': 25.6873668671 * units.VEGAMAG + 0 * u.mag('s/electron'), 'e_zP': 0.0027553202* units.VEGAMAG + 0 * u.mag('s/electron'), 'lamb': 6217.51 * u.AA, 'f_zp_nu': 3228.75 * u.Jy}, 'AB': {'zp': 25.8010446445 * u.Unit("mag(AB s/electron)"), 'e_zP': 0.0027590522 * u.Unit("mag(AB s/electron)"), 'lamb': 6217.51 * u.AA, 'f_zp_nu': 3631 * u.Jy}}, 
        'BP': {'VEGAMAG': {'zp': 25.3385422158 * units.VEGAMAG + 0 * u.mag('s/electron'), 'e_zP': 0.0027901700* units.VEGAMAG + 0 * u.mag('s/electron'), 'lamb': 5109.71 * u.AA, 'f_zp_nu': 3552.01 * u.Jy}, 'AB': {'zp': 25.3539555559 * u.Unit("mag(AB s/electron)"), 'e_zP': 0.0023065687 * u.Unit("mag(AB s/electron)"), 'lamb': 5109.71 * u.AA, 'f_zp_nu': 3631 * u.Jy}}, 
        'RP': {'VEGAMAG': {'zp': 24.7478955012 * units.VEGAMAG + 0 * u.mag('s/electron'), 'e_zP': 0.0037793818* units.VEGAMAG + 0 * u.mag('s/electron'), 'lamb': 7769.02 * u.AA, 'f_zp_nu': 2554.95 * u.Jy}, 'AB': {'zp': 25.1039837393 * u.Unit("mag(AB s/electron)"), 'e_zP': 0.0015800349 * u.Unit("mag(AB s/electron)"), 'lamb': 7769.02 * u.AA, 'f_zp_nu': 3631 * u.Jy}}, 
        },
    'graphic': {'y': {'colname': 'flux', 'units': 'electron/s'},
              'y_err': {'colname': 'flux_error', 'units': 'electron/s'},
              'cid': 'band',
              'multi': None}    
        },
'jwst': {
    'system': 'VEGAMAG',
    'zeropt': {'P750L': {'VEGAMAG': { 'lamb': 8000 * u.AA}}
               },
    'graphic': {'y': {'colname': 'FLUX', 'units': 'Jy'},
              'y_err': {'colname': 'FLUX_ERROR', 'units': 'Jy'},
              'cid': None,
              'multi': 'FILTER',
              'cextra': 'WAVELENGTH'}    
        },
'xmm': {
    'system': 'VEGAMAG',
    'zeropt': {'EMOS1': {'VEGAMAG': { 'lamb': 200 * u.eV}},
               'EMOS2': {'VEGAMAG': { 'lamb': 200 * u.eV}},
               'EPN': {'VEGAMAG': { 'lamb': 200 * u.eV}},
               },
    'graphic': {'y': {'colname': 'RATE', 'units': 'cts/s'},
              'y_err': {'colname': 'ERROR', 'units': 'cts/s'},
              'cid': None,
              'multi': 'INSTRUME'}    
        }
}