import astropy.units as u
from astropy.units import Unit, Quantity
import copy
from functools import reduce
from pprint import pp as pp
import numpy as np
from astropy.table import Table
from synphot import units, SourceSpectrum

def my_key_exists(key, var):
    '''Check if any key exists in nested dictionary (var) and break upon the first occurence.
    It needs to keep track of the res status 
    Note: key in var.keys() only applies to first level'''
    res = False
    if isinstance(var, dict):
        for k, v in var.items():
            if k == key:
                return True
            if isinstance(v, dict):
                res = my_key_exists(key, v)
                if res:
                  break                
    return res

def dict_key_filter(obj, obj_filter):
    '''Filter dictionary by keys through a nested dictionary https://stackoverflow.com/questions/31710271/how-to-filter-by-keys-through-a-nested-dictionary-in-a-pythonic-way'''
    def inner_dict_key_filter(obj): return dict_key_filter(obj, obj_filter)
    def to_keep(subtree): return not isinstance(subtree, (dict, list)) or subtree

    def build_subtree(key, value):
        if key in obj_filter:
            return copy.deepcopy(value) # keep the branch
        elif isinstance(value, (dict, list)):
            return inner_dict_key_filter(value) # continue to search
        return [] # just an orphan value here

    if isinstance(obj, dict):
        key_subtree_pairs = ((key, build_subtree(key, value)) for key, value in obj.items())
        return {key:subtree for key, subtree in key_subtree_pairs if to_keep(subtree)}
    elif isinstance(obj, list):
        return list(filter(to_keep, map(inner_dict_key_filter, obj)))
    return []

def gen_dict_extract_slow(var, keys):
  '''Traverse to extract the values for all occurences of a list of keys https://stackoverflow.com/questions/9807634/find-all-occurrences-of-a-key-in-nested-dictionaries-and-lists'''
  for key in keys:
      if hasattr(var, 'items'):
         for k, v in var.items():
            if k == key:
               yield v
            if isinstance(v, dict):
               for result in gen_dict_extract_slow(v, [key]):
                  yield result
            elif isinstance(v, list):
               for d in v:
                  for result in gen_dict_extract_slow(d, [key]):
                     yield result    

def gen_dict_extract_fast(key, var):
    '''Traverse to extract the values for all occurences of a key https://stackoverflow.com/questions/9807634/find-all-occurrences-of-a-key-in-nested-dictionaries-and-lists
    Note: does not work to extract a subdictionary by key d = list(*gen_dict_extract_fast('G', gaia))'''
    if hasattr(var,'items'): # hasattr(var,'items') for python 3
        for k, v in var.items(): # var.items() for python 3
            if k == key:
                yield v
            if isinstance(v, dict):
                for result in gen_dict_extract_fast(key, v):
                    yield result
            elif isinstance(v, list):
                for d in v:
                    for result in gen_dict_extract_fast(key, d):
                        yield result

def gen_dict_extract_value(key, var):
    '''Traverse to extract the values for all occurences of a key https://stackoverflow.com/questions/9807634/find-all-occurrences-of-a-key-in-nested-dictionaries-and-lists
    Note: does not work to extract a subdictionary by key d = list(*gen_dict_extract_fast('G', gaia))'''
    if hasattr(var,'items'): # hasattr(var,'items') for python 3
        for k, v in var.items(): # var.items() for python 3
            if k == key:
                yield v.value
            if isinstance(v, dict):
                for result in gen_dict_extract_value(key, v):
                    yield result
            elif isinstance(v, list):
                for d in v:
                    for result in gen_dict_extract_value(key, d):
                        yield result
def gen_dict_extract_unit(key, var):
    '''Traverse to extract the values for all occurences of a key https://stackoverflow.com/questions/9807634/find-all-occurrences-of-a-key-in-nested-dictionaries-and-lists
    Note: does not work to extract a subdictionary by key d = list(*gen_dict_extract_fast('G', gaia))'''
    if hasattr(var,'items'): # hasattr(var,'items') for python 3
        for k, v in var.items(): # var.items() for python 3
            if k == key:
                yield v.unit
            if isinstance(v, dict):
                for result in gen_dict_extract_unit(key, v):
                    yield result
            elif isinstance(v, list):
                for d in v:
                    for result in gen_dict_extract_unit(key, d):
                        yield result


def main():
    gaia = {
    'band': ['G', 'BP', 'RP'],
    'zeropt': {'G': {'VEGAMAG': {'zp': 25.6873668671 * units.VEGAMAG, 'e_zP': 0.0027553202* units.VEGAMAG, 'lamb': 6251.50 * u.AA, 'f_zp_nu': 3599.12 * u.Jy}, 'AB': {'zp': 25.8010446445 * u.ABmag, 'e_zP': 0.0027590522 * u.ABmag, 'lamb': 6251.50 * u.AA, 'f_zp_nu': 3631 * u.Jy}}, 
          'BP': {'VEGAMAG': {'zp': 25.3385422158 * units.VEGAMAG, 'e_zP': 0.0027901700* units.VEGAMAG, 'lamb': 5124.20 * u.AA, 'f_zp_nu': 3256.01 * u.Jy}, 'AB': {'zp': 25.3539555559 * u.ABmag, 'e_zP': 0.0023065687 * u.ABmag, 'lamb': 5124.20 * u.AA, 'f_zp_nu': 3631 * u.Jy}}, 
          'RP': {'VEGAMAG': {'zp': 24.7478955012 * units.VEGAMAG, 'e_zP': 0.0037793818* units.VEGAMAG, 'lamb': 7829.65 * u.AA, 'f_zp_nu': 2530.35 * u.Jy}, 'AB': {'zp': 25.1039837393 * u.ABmag, 'e_zP': 0.0015800349 * u.ABmag, 'lamb': 7829.65 * u.AA, 'f_zp_nu': 3631 * u.Jy}}, 
           },
    'graphic': {'y': {'colname': 'flux', 'units': 'electron/s'},
            'y_err': {'colname': 'flux_error', 'units': 'electron/s'}}    
}

    tbl = Table.read('/Users/epuga/ESDC/TSViz/data/gaia/anonymous1690191210843_fits/EPOCH_PHOTOMETRY-Gaia DR3 4057091150787104896.fits', astropy_native=True)
    
    SYSTEM = 'AB'

    y_col = gaia['graphic']['y']['colname'] # has y_col.name and y_col.unit
    if my_key_exists('y_err', gaia):
        y_err_col = gaia['graphic']['y_err']['colname']
    if tbl[y_col].unit == "'electron'.s**-1":
        tbl[y_col].unit = "electron/s"
    if my_key_exists('y_err', gaia):
        if tbl[y_err_col].unit == "'electron'.s**-1":
            tbl[y_err_col].unit = "electron/s"
        
    print(tbl)

    
    original_unit = tbl[y_col].unit
    target_unit = u.mJy
    y_col_target = '{0}_{1}'.format(y_col, target_unit)
    if my_key_exists('y_err', gaia):
        y_err_col_target = '{0}_err_{1}'.format(y_col, target_unit)
    
    indices = ['band']

    tbl.add_column(1., name=y_col_target)
    #tbl.add_column(1., name=y_err_col_target)
    try:
        tbl[y_col_target] = tbl[y_col].to(target_unit)
    except:
        if my_key_exists('zp', gaia):
            gaia_system = dict_key_filter(gaia, [SYSTEM])
            zp_unit = Unit(*gen_dict_extract_unit('zp', gaia_system))
            if isinstance(zp_unit, u.MagUnit):
                if isinstance(tbl[y_col], u.Magnitude):
                    tbl[y_col_target] = tbl[y_col].quantity
                else:
                    tbl[y_col_target] = u.Magnitude(tbl[y_col].quantity)
                for index in indices:
                    # the index is tabulated in tbl
                    tbl_by_index = tbl.group_by(index)
                    #do we have the value of those indices as key in the data and calibration dictionary?
                    if my_key_exists(index, gaia):
                        for key, group in zip(tbl_by_index.groups.keys, tbl_by_index.groups):
                            id = key[index]
                            subdict = dict(*gen_dict_extract_slow(gaia_system, ([id])))
                            print(float(*gen_dict_extract_value('zp', subdict)))
                            print(u.Magnitude(*gen_dict_extract_fast('zp', subdict)) )
                            print(group[y_col_target].value)
                            group[y_col_target][:] = group[y_col_target].quantity + u.Magnitude(*gen_dict_extract_fast('zp', subdict)) 
                            print(group[y_col_target].value)
                            # if you do not have m_err, and need to derive from f_err and f in other units
                            #err_mag_system_filter = np.sqrt((-2.5 / np.log(10) * group[y_err_col].quantity / group[y_col].quantity)**2 + float(*gen_dict_extract_value('e_zp', subdict))**2)
                if target_unit.is_equivalent(u.mag): 
                    tbl[y_col_target] = tbl[y_col_target].value * zp_unit
                    #tbl[y_err_col_target] = tbl[y_err_col_target].value * zp_unit
                else:   
                    print(tbl) 
                    tbl[y_col_target] = (tbl[y_col_target].value * zp_unit).to(target_unit)
                    #tbl[y_err_col_target] = (tbl[y_err_col_target].value / 2.5 * np.log(10) * zp_unit).to(target_unit)
    
            else:
                #TODO: just do the
                pass 
        else:
            #There is no zeropoint and what is failing is the equivalency
            pass
    print(tbl)
            
    
# if not isinstance(wavelengths, u.Quantity):
#         wavelengths = wavelengths * u.AA

#     eqv = u.spectral_density(wavelengths)

#     # Use built-in astropy equivalencies
#     try:
#         out_flux = fluxes.to(out_flux_unit, eqv)


if __name__ == '__main__':
    main()
    pass
