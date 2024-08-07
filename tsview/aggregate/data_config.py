
import astropy.units as u
from synphot import units

# zeropt: cid/multi, PHOT_SYSTEM, data_dict_field
# '**.{{}}.**.{0}.**.{1}'.format('VEGAMAG', 'e_zP')


DATA_DICT = {
'gaia': {
    'system': 'VEGAMAG',
    'expr': '**.{{}}.**.{0}.**.{1}',
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
    'system': None,
    'expr': '**.{{}}.**.{1}',
    'zeropt': {'P750L': { 'lamb': 8000 * u.AA}
               },
    'graphic': {'y': {'colname': 'FLUX', 'units': 'Jy'},
              'y_err': {'colname': 'FLUX_ERROR', 'units': 'Jy'},
              'cid': None,
              'multi': 'FILTER',
              'cextra': 'WAVELENGTH'}    
        },

'xmm-epic': {
    'system': None,
    'expr': '**.{{}}.**.{1}',
    'zeropt': {'EMOS1': { 'lamb': 200 * u.eV},
               'EMOS2': { 'lamb': 200 * u.eV},
               'EPN': { 'lamb': 200 * u.eV}
               },
    'graphic': {'y': {'colname': 'RATE', 'units': 'cts/s'},
              'y_err': {'colname': 'ERROR', 'units': 'cts/s'},
              'cid': None,
              'multi': 'INSTRUME'}    
        }
}
# zeropt: cid/multi, PHOT_SYSTEM, data_dict_field
# '**.{{}}.**.{1}'.format(None, 'lamb')