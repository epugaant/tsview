import astropy.units as u
from astropy.units import Unit, Quantity
import copy
from functools import reduce
from pprint import pp as pp
import numpy as np
from astropy.table import Table
from synphot import units, SourceSpectrum

def keys_exists(element, *keys):
    '''
    Check if *keys (nested) exists in `element` (dict). https://stackoverflow.com/questions/43491287/elegant-way-to-check-if-a-nested-key-exists-in-a-dict
    '''
    if not isinstance(element, dict):
        raise AttributeError('keys_exists() expects dict as first argument.')
    if len(keys) == 0:
        raise AttributeError('keys_exists() expects at least two arguments, one given.')

    _element = element
    for key in keys:
        try:
            _element = _element[key]
        except KeyError:
            return False
    return True



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

#orig
def gen_dict_extract(var, key):
    if isinstance(var, dict):
        for k, v in var.items():
            if k == key:
                yield v
            if isinstance(v, (dict, list)):
                yield from gen_dict_extract(v, key)
    elif isinstance(var, list):
        for d in var:
            yield from gen_dict_extract(d, key)

            
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

def getpath(nested_dict, value, prepath=()):
    '''Search for a unique value and get the parent dictionary names (keys)  https://stackoverflow.com/questions/22162321/search-for-a-value-in-a-nested-dictionary-python'''
    for k, v in nested_dict.items():
        path = prepath + (k,)
        if v == value: # found value
            return path
        elif hasattr(v, 'items'): # v is a dict
            p = getpath(v, value, path) # recursive call
            if p is not None:
                return p
            
def find_paths(nested_dict, value, prepath=()):
    '''Search  for multiple occurences of a value and get the parent dictionary names (keys)  https://stackoverflow.com/questions/22162321/search-for-a-value-in-a-nested-dictionary-python'''
    for k, v in nested_dict.items():
        path = prepath + (k,)
        if v == value: # found value
            yield path
        elif hasattr(v, 'items'): # v is a dict
            yield from find_paths(v, value, path) 

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

def extract(var, key, context_keys=(), xpath=''):
    '''XPath to a given key bottom https://stackoverflow.com/questions/9807634/find-all-occurrences-of-a-key-in-nested-dictionaries-and-lists'''
    if isinstance(var, dict):
        if key in var:
            yield {f'{xpath}.{key}': var[key]} | {f'{xpath}.{key}': value for key, value in var.items() if key in context_keys}
        for subkey, value in var.items():
            yield from extract(value, key, context_keys, f'{xpath}.{subkey}')
    elif isinstance(var, list):
        for i, elem in enumerate(var):
            yield from extract(elem, key, context_keys, f'{xpath}[{i}]')
 
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
        # 'zeropt': {'G': {'VEGA': {'zp': 25.6873668671, 'e_zp': 0.0027553202, 'lamb_AA': 6251.50, 'f_zp_nu_Jy': 3599.12}, 'AB': {'zp': 25.8010446445, 'e_zp': 0.0027590522, 'lamb_AA': 6251.50, 'f_zp_nu_Jy': 3631}}, 
        #         'BP': {'VEGA': {'zp': 25.3385422158, 'e_zp': 0.0027901700, 'lamb_AA': 5124.20, 'f_zp_nu_Jy': 3256.01}, 'AB': {'zp': 25.3539555559, 'e_zp': 0.0023065687, 'lamb_AA': 5124.20, 'f_zp_nu_Jy': 3631}}, 
        #         'RP': {'VEGA': {'zp': 24.7478955012, 'e_zp': 0.0037793818, 'lamb_AA': 7829.65, 'f_zp_nu_Jy': 2530.35}, 'AB': {'zp': 25.1039837393, 'e_zp': 0.0015800349, 'lamb_AA': 7829.65, 'f_zp_nu_Jy': 3631}}, 
        #         },
        'zeropt': {'G': {'VEGAMAG': {'zp': 25.6873668671 * units.VEGAMAG, 'e_zP': 0.0027553202* units.VEGAMAG, 'lamb': 6251.50 * u.AA, 'f_zp_nu': 3599.12 * u.Jy}, 'AB': {'zp': 25.8010446445 * u.ABmag, 'e_zP': 0.0027590522 * u.ABmag, 'lamb': 6251.50 * u.AA, 'f_zp_nu': 3631 * u.Jy}}, 
          'BP': {'VEGAMAG': {'zp': 25.3385422158 * units.VEGAMAG, 'e_zP': 0.0027901700* units.VEGAMAG, 'lamb': 5124.20 * u.AA, 'f_zp_nu': 3256.01 * u.Jy}, 'AB': {'zp': 25.3539555559 * u.ABmag, 'e_zP': 0.0023065687 * u.ABmag, 'lamb': 5124.20 * u.AA, 'f_zp_nu': 3631 * u.Jy}}, 
          'RP': {'VEGAMAG': {'zp': 24.7478955012 * units.VEGAMAG, 'e_zP': 0.0037793818* units.VEGAMAG, 'lamb': 7829.65 * u.AA, 'f_zp_nu': 2530.35 * u.Jy}, 'AB': {'zp': 25.1039837393 * u.ABmag, 'e_zP': 0.0015800349 * u.ABmag, 'lamb': 7829.65 * u.AA, 'f_zp_nu': 3631 * u.Jy}}, 
           },
        'graphic': {'y': {'colname': 'flux', 'units': 'electron/s'},
                'y_err': {'colname': 'flux_error', 'units': 'electron/s'}}    
    }
    
    SYSTEM = 'AB'
    #TODO: implement workflow for VEGAMAG
    filter = 'G'
    # Does the field exist for this mission?
    # top-down levels
    print('zeropt (exists): {}'.format(keys_exists(gaia, "zeropt")))
    # any level
    print('{0} exists: {1}'.format('zp', my_key_exists('zp', gaia)))
    print('')
          
    #filter dictionary with keys
    gaia_sel = dict_key_filter(gaia, [filter])
    pp(gaia_sel)
    gaia_sel2 = dict_key_filter(gaia_sel, [SYSTEM])
    pp(gaia_sel2)
    pp(dict_key_filter(dict_key_filter(gaia, [filter]), [SYSTEM]))
    print('')
    
    # can we get the value of a key(s)?
    # What is the value from a (single/multiple) key, easier to get them in a final list
    # for x in gen_dict_extract_fast('band', gaia): #awful way to display, keep it to be informative
    #     print(x)
    # for x in gen_dict_extract_slow(gaia, (['band'])):
    #     print(x)
   
    print(*gen_dict_extract_fast('band', gaia)) #nifty way to display
    print(*gen_dict_extract_slow(gaia, (['zp', 'lamb'])))
    print(*gen_dict_extract_slow(gaia, (['G'])))
    print('Select System, filter and zp')
    pp(dict_key_filter(gaia, [SYSTEM]))
    d = dict(*gen_dict_extract_slow(dict_key_filter(gaia, [SYSTEM]), (['G'])))
    print(d)
    print(*gen_dict_extract_unit('zp', d))
    print(*gen_dict_extract_fast('zp', d))
    print('')
    
    
    # What is the dictionary value call for one occurence?
    print(getpath(gaia, 25.6873668671 *units.VEGAMAG))
    print(*find_paths(gaia, 'flux'))
    
    # What is the dictionary key xpath call? (interesting, but not so useful with dictionaries)
    pp(list(extract(gaia, 'zp')))
    # for various keys
    pp(reduce(lambda acc, elem: acc | elem, extract(gaia, 'zp', 'lamb_AA')))


    gaia = {
    'band': ['G', 'BP', 'RP'],
    # 'zeropt': {'G': {'VEGA': {'zp': {'value': 25.6873668671, 'unit': 'mag'}, 'e_zp': {'value': 0.0027553202, 'unit': 'mag'}, 'lamb': {'value': 6251.50, 'unit': 'Angstrom'}, 'f_zp_nu': {'value': 3599.12, 'unit': 'Jy'}}, 'AB': {'zp': {'value': 25.8010446445, 'unit': 'mag'}, 'e_zp': {'value': 0.0027590522, 'unit': 'mag'}, 'lamb': {'value': 6251.50, 'unit': 'Angstrom'}, 'f_zp_nu': {'value': 3631, 'unit': 'Jy'}}}, 
    #         'BP': {'VEGA': {'zp': {'value': 25.3385422158, 'unit': 'mag'}, 'e_zp': {'value': 0.0027901700, 'unit': 'mag'}, 'lamb': {'value': 5124.20, 'unit': 'Angstrom'}, 'f_zp_nu': {'value': 3256.01, 'unit': 'Jy'}}, 'AB': {'zp': {'value': 25.3539555559, 'unit': 'mag'}, 'e_zp': {'value': 0.0023065687, 'unit': 'mag'}, 'lamb': {'value': 5124.20, 'unit': 'Angstrom'}, 'f_zp_nu': {'value': 3631, 'unit': 'Jy'}}}, 
    #         'RP': {'VEGA': {'zp': {'value': 24.7478955012, 'unit': 'mag'}, 'e_zp': {'value': 0.0037793818, 'unit': 'mag'}, 'lamb': {'value': 7829.65, 'unit': 'Angstrom'}, 'f_zp_nu': {'value': 2530.35, 'unit': 'Jy'}}, 'AB': {'zp': {'value': 25.1039837393, 'unit': 'mag'}, 'e_zp': {'value': 0.0015800349, 'unit': 'mag'}, 'lamb': {'value': 7829.65, 'unit': 'Angstrom'}, 'f_zp_nu': {'value': 3631, 'unit': 'Jy'}}}, 
    #         },
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
            
    if not original_unit.is_equivalent(u.mag):
        tbl.add_column(1., name=y_col_target)
        tbl.add_column(1., name=y_err_col_target)
        # add zeropoints to transform to mag
        # zeropoints need to be added by band
        # indices are the columns that contain grouping indices and will correspond to different artists in the plot
        indices = ['band']
        for index in indices:
            # the index is tabulated in tbl
            tbl_by_index = tbl.group_by(index)
            #do we have the value of those indices as key in the data and calibration dictionary?
            if my_key_exists(index, gaia):
                filters = list(*gen_dict_extract_fast(index, gaia))
                
                if my_key_exists('zp', gaia):
                    # for filter in filters:
                    for key, group in zip(tbl_by_index.groups.keys, tbl_by_index.groups):
                        id = key[index]
                        subdict = dict(*gen_dict_extract_slow(dict_key_filter(gaia, [SYSTEM]), ([id])))
                        print(float(*gen_dict_extract_value('zp', subdict)))
                        print()
                        mag_system_filter = u.Magnitude(group[y_col].quantity) + float(*gen_dict_extract_value('zp', subdict)) * u.Unit('mag ({0})'.format(SYSTEM))
                        # if you do not have m_err, and need to derive from f_err and f in other units
                        err_mag_system_filter = np.sqrt((-2.5 / np.log(10) * group[y_err_col].quantity / group[y_col].quantity)**2 + float(*gen_dict_extract_value('e_zp', subdict))**2)
                        # parse to the new added column
                        if target_unit.is_equivalent(u.mag): 
                            group[y_col_target] = mag_system_filter.value * u.Unit('mag ({0})'.format(SYSTEM))
                            group[y_err_col_target] = err_mag_system_filter.value * u.Unit('mag ({0})'.format(SYSTEM))
                        else:    
                            group[y_col_target] = (mag_system_filter.value * u.Unit('mag ({0})'.format(SYSTEM))).to(target_unit)
                            group[y_err_col_target] = (err_mag_system_filter.value / 2.5 * np.log(10) * u.Unit('mag ({0})'.format(SYSTEM))).to(target_unit)
                    
                    # 
                    # add selected zp to entire subcolumn
    print()   
                    
        
    
    # exists?
    # how many bands there are? bands
    # extract value for every band
# if not isinstance(wavelengths, u.Quantity):
#         wavelengths = wavelengths * u.AA

#     eqv = u.spectral_density(wavelengths)

#     # Use built-in astropy equivalencies
#     try:
#         out_flux = fluxes.to(out_flux_unit, eqv)


if __name__ == '__main__':
    main()
    pass

def get(data, *keys):
    for key in keys:
        if key in data:
            data = data[key]
        else:
            return False
    return data

thing = {'outerkey':{'innerkey': 42, 'mostinnerkey': {'value': 5}}}
print(get(thing, 'outerkey', 'innerkey')) # prints 42
print(get(thing, 'outerkey', 'mostinnerkey', 'value')) # prints 42
print(get(thing, 'outerkey', 'invalid')) # prints False
print(get(thing, 'invalid', 'innerkey')) # prints False

print(get(thing, *('outerkey', 'mostinnerkey', 'value')))