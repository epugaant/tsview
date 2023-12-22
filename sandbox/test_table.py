from astropy.table import Table, Column
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
        print(key[cid])
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
        print(key[cid])
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
    # t = Table([['G', 'RP'], [3, 4], [5, 6]], names=('band', 'b', 'flux'))
    # t['flux'].unit = "electron/s"
    band = ['BP', 'G', 'RP']
    mag = [12.751842803737022, 12.415125851786478, 11.902541753829738] * u.mag
    flux = [108312.79547960439, 203655.62323629577, 137448.97291377262] * u.electron/u.s
    flux_error = [181.4765160322774, 76.00599243480103, 183.6823434658511] * u.electron/u.s
    t = Table([band, mag, flux, flux_error], names=('band', 'mag', 'flux', 'flux_error'))

    
    #t_mod = create_mod_column_using_index(t, 'band', 'flux', 'flux_new')
    #print(t_mod)
    flux_col = 'flux'
    flux_col_new = flux_col+'_new'
    #orig_unit = u.Unit( "electron/s") #u.Unit("count/s")
    orig_unit = t[flux_col].unit
    target_unit = u.mJy
    sys = 'VEGAMAG'
    zpt = 'zp' # None if there is no zpt
    
    flux_col = 'mag'
    flux_col_new = flux_col+'_new'
    #orig_unit = u.Unit( "electron/s") #u.Unit("count/s")
    orig_unit = t[flux_col].unit
    target_unit = u.mJy
    sys = 'VEGAMAG'
    zpt = 'f_zp_nu' # None if there is no zpt

    if orig_unit.is_equivalent(u.mag):
        orig_unit = u.mag()
        
    # No conversion necessary
    if orig_unit == target_unit:
        return t[flux_col]
    
    # orig_unit (electron/s; count/s)------ zp_unit (mag)
    # orig_unit (&mag) --- zp_unit (Jy)
    if zpt:
        if orig_unit.physical_type in ('unknown', 'dimensionless'):
            if isinstance(orig_unit, u.MagUnit):
                #Initialize column with given value
                t[flux_col_new] = u.Dex(-0.4*t[flux_col].value).physical #this is just a quantity
            else:
                #Initialize column with given value
                t[flux_col_new] = u.Magnitude(t[flux_col].quantity)
            #check now zp_unit
            expr = '**.{{}}.**.{0}.**.{1}'.format(sys, zpt) # escaped for key=band
            if all(isinstance(x.unit, u.MagUnit) for x in glom.glom(gaia, '**.{0}.{1}'.format(sys, zpt))):
                key_indices = create_mod_column_using_index(t, 'band', flux_col_new, flux_col_new, expr, add)
            else:
                key_indices = create_mod_column_using_index(t, 'band', flux_col_new, flux_col_new, expr, multiply)
    else:
        print('No zeropoint step. Directly goint to conversion')
    
    wave = col_from_index(t, 'wave', 'band', '**.{{}}.**.{0}.**.{1}'.format('AB', 'lamb'))
    # if not isinstance(wave, (Column, u.Quantity)):
    #     wave = wave * u.AA
    eqv = u.spectral_density(wave)
        
    #convert unit from intermediate zp_unit to target_unit
    try:
        t[flux_col_new] = t[flux_col_new].to(target_unit, eqv)
        print(t)
    except u.core.UnitConversionError:
        #TODO: I could simplify the fist try wihouth equivalencies and move the wavelength in here
        if t[flux_col_new].unit == units.VEGAMAG or target_unit == units.VEGAMAG:
            vega = SourceSpectrum.from_vega()  # For unit conversion  
            t[flux_col_new] = units.convert_flux(wave.quantity, t[flux_col_new].quantity, target_unit, vegaspec=vega)
        else:
            t[flux_col_new] = units.convert_flux(wave.quantity, t[flux_col_new].quantity, target_unit) 
        print(t)
        return t[flux_col_new]
    
    
        
if __name__ == '__main__':
    main()
    # glomexpr='**.{0}.**.{1}.**.{2}'.format((key, system, field))
    # target= glom.glom(gaia, glomexpr)
    # [val] = glom.flatten(target, levels=2)
    # #thanks to https://stackoverflow.com/questions/22219004/how-to-group-dataframe-rows-into-list-in-pandas-groupby
    # df2 = df.groupby('band')[['time', 'flux', 'flux_error']].agg(list).rename(columns={'time': 'x', 'flux': 'y', 'flux_error': 'error_y'}).to_json(orient='index')

    pass




