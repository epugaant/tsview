from astropy.table import Table
import numpy as np
import astropy.units as u
from synphot import units
import glom
from synphot.spectrum import SourceSpectrum

gaia = {'zeropt': {'G': {'VEGAMAG': {'zp': 25.6873668671 * units.VEGAMAG, 'e_zP': 0.0027553202* units.VEGAMAG, 'lamb': 6251.50 * u.AA, 'f_zp_nu': 3599.12 * u.Jy}, 'AB': {'zp': 25.8010446445 * u.ABmag, 'e_zP': 0.0027590522 * u.ABmag, 'lamb': 6251.50 * u.AA, 'f_zp_nu': 3631 * u.Jy}}, 
        'BP': {'VEGAMAG': {'zp': 25.3385422158 * units.VEGAMAG, 'e_zP': 0.0027901700* units.VEGAMAG, 'lamb': 5124.20 * u.AA, 'f_zp_nu': 3256.01 * u.Jy}, 'AB': {'zp': 25.3539555559 * u.ABmag, 'e_zP': 0.0023065687 * u.ABmag, 'lamb': 5124.20 * u.AA, 'f_zp_nu': 3631 * u.Jy}}, 
        'RP': {'VEGAMAG': {'zp': 24.7478955012 * units.VEGAMAG, 'e_zP': 0.0037793818* units.VEGAMAG, 'lamb': 7829.65 * u.AA, 'f_zp_nu': 2530.35 * u.Jy}, 'AB': {'zp': 25.1039837393 * u.ABmag, 'e_zP': 0.0015800349 * u.ABmag, 'lamb': 7829.65 * u.AA, 'f_zp_nu': 3631 * u.Jy}}, 
        }
}

def add(x, val):
    return x + val
def multiply(x, val):
    #TODO: If x is a DexUnit, use the physical attribute to allow the conversion.
    return x * val

def create_mod_column_using_index(tbl, cid, colname, new_colname, expr, func):
    '''Create a new column and update the rows according to an index column'''
    tbl.add_index(cid)
    orig_unit = tbl[colname].unit
    key_indices = {}
    for key in tbl.group_by(cid).groups.keys:
        print(key)
        indx = tbl.loc_indices[key[cid]]
        #Astropy User Guide recommendation to modify Table using table[column][row] order
        target =  glom.glom(gaia, expr.format(key[cid]))
        #extract one element from a list
        [zp] = glom.flatten(target, levels=(expr.count('**') - 1))
        res = func(tbl[colname][indx] * orig_unit, zp)
        tbl[new_colname][indx] = res.value
        key_indices[key[cid]] = indx
    tbl[new_colname].unit = zp.unit
    return key_indices

def col_from_index(tbl, colname,  cid,  expr):
    '''Return a wavelength column from a column indexed'''
    #TODO: check that the column is a Quantity
    tbl.add_index(cid)
    tbl[colname] = np.nan
    for key in tbl.group_by(cid).groups.keys:
        print(key)
        indx = tbl.loc_indices[key[cid]]
        #Astropy User Guide recommendation to modify Table using table[column][row] order
        target =  glom.glom(gaia, expr.format(key[cid]))
        [value] = glom.flatten(target, levels=(expr.count('**') - 1))
        tbl[colname][indx] = value.value
    tbl[colname].unit = value.unit
    col = tbl[colname]
    tbl.remove_column(colname)
    return col

def main():        
    t = Table([['G', 'RP'], [3, 4], [5, 6]], names=('band', 'b', 'c'))
    t['c'].unit = "electron/s"
    
    #t_mod = create_mod_column_using_index(t, 'band', 'c', 'c_new')
    #print(t_mod)

    orig_unit = u.Unit( "electron/s") #u.Unit("count/s")
    target_unit = u.mJy
    sys = 'AB'
    zpt = 'zp'
    
    # No conversion necessary
    if orig_unit == target_unit:
        return t['c']
    
    # orig_unit (electron/s; count/s)------ zp_unit (mag)
    if orig_unit.physical_type == 'unknown' and all(isinstance(x.unit, u.MagUnit) for x in glom.glom(gaia, '**.{0}.{1}'.format(sys, zpt))):
        #Initialize column with given value
        t['c_new'] = u.Magnitude(t['c'].quantity)
        expr = '**.{{}}.**.{0}.**.{1}'.format(sys, zpt)
        key_indices = create_mod_column_using_index(t, 'band', 'c_new', 'c_new', expr, add)
    
    # orig_unit (&mag) --- zp_unit (Jy)
    if orig_unit.physical_type == 'spectral flux density' and isinstance(orig_unit, u.MagUnit) and not all(isinstance(x.unit, u.MagUnit) for x in glom.glom(gaia, '**.{0}.{1}'.format(sys, zpt))): 
        #Initialize column with given value
        t['c_new'] = u.Dex(-0.4*t['c'].quantity).physical #this is just a quantity
        expr = '**.{{}}.**.{0}.**.{1}'.format(sys, zpt)
        key_indices = create_mod_column_using_index(t, 'band', 'c_new', 'c_new', expr, multiply)
    
    wave = col_from_index(t, 'wave', 'band', '**.{{}}.**.{0}.**.{1}'.format('AB', 'lamb'))
    if not isinstance(wave, u.Quantity):
        wave = wave * u.AA
    eqv = u.spectral_density(wave)
        
    #convert unit
    try:
        t['c_new'] = t['c_new'].to(target_unit, eqv)
        pass
    except u.core.UnitConversionError:
        
        if t['c_new'].unit == units.VEGAMAG:
            vega = SourceSpectrum.from_vega()  # For unit conversion  
            t['c_new'].quantity = units.convert_flux(wave.quantity, t['c_new'].quantity, target_unit, vegaspec=vega)
        else:
            t['c_new'].quantity = units.convert_flux(wave.quantity, t['c_new'].quantity, target_unit) 
        return t['c_new']
        
if __name__ == '__main__':
    main()
    glomexpr='**.{0}.**.{1}.**.{2}'.format((key, system, field))
    target= glom.glom(gaia, glomexpr)
    [val] = glom.flatten(target, levels=2)
    #thanks to https://stackoverflow.com/questions/22219004/how-to-group-dataframe-rows-into-list-in-pandas-groupby
    df2 = df.groupby('band')[['time', 'flux', 'flux_error']].agg(list).rename(columns={'time': 'x', 'flux': 'y', 'flux_error': 'error_y'}).to_json(orient='index')

    pass


band = ['BP', 'G', 'RP']
mag = [12.751842803737022, 12.415125851786478, 11.902541753829738] * u.mag
flux = [108312.79547960439, 203655.62323629577, 137448.97291377262] * u.electron/u.s
flux_error = [181.4765160322774, 76.00599243480103, 183.6823434658511] * u.electron/u.s
t = Table([band, mag, flux, flux_error], names=('band', 'mag', 'flux', 'flux_error'))


