from dataclasses import dataclass

import glom
import numpy as np

from astropy.time import Time
from astropy.table import Table
import astropy.units as u
from synphot import units
from synphot.spectrum import SourceSpectrum


gaia = {
    'zeropt': {'G': {'VEGAMAG': {'zp': 25.6873668671 * units.VEGAMAG + 0 * u.mag('s/electron'), 'e_zP': 0.0027553202* units.VEGAMAG + 0 * u.mag('s/electron'), 'lamb': 6217.51 * u.AA, 'f_zp_nu': 3228.75 * u.Jy}, 'AB': {'zp': 25.8010446445 * u.Unit("mag(AB s/electron)"), 'e_zP': 0.0027590522 * u.Unit("mag(AB s/electron)"), 'lamb': 6217.51 * u.AA, 'f_zp_nu': 3631 * u.Jy}}, 
        'BP': {'VEGAMAG': {'zp': 25.3385422158 * units.VEGAMAG + 0 * u.mag('s/electron'), 'e_zP': 0.0027901700* units.VEGAMAG + 0 * u.mag('s/electron'), 'lamb': 5109.71 * u.AA, 'f_zp_nu': 3552.01 * u.Jy}, 'AB': {'zp': 25.3539555559 * u.Unit("mag(AB s/electron)"), 'e_zP': 0.0023065687 * u.Unit("mag(AB s/electron)"), 'lamb': 5109.71 * u.AA, 'f_zp_nu': 3631 * u.Jy}}, 
        'RP': {'VEGAMAG': {'zp': 24.7478955012 * units.VEGAMAG + 0 * u.mag('s/electron'), 'e_zP': 0.0037793818* units.VEGAMAG + 0 * u.mag('s/electron'), 'lamb': 7769.02 * u.AA, 'f_zp_nu': 2554.95 * u.Jy}, 'AB': {'zp': 25.1039837393 * u.Unit("mag(AB s/electron)"), 'e_zP': 0.0015800349 * u.Unit("mag(AB s/electron)"), 'lamb': 7769.02 * u.AA, 'f_zp_nu': 3631 * u.Jy}}, 
        },
    'graphic': {'y': {'colname': 'flux', 'units': 'electron/s'},
              'y_err': {'colname': 'flux_error', 'units': 'electron/s'}}    
}

def column_from_dict_index(tbl, cid, data_dict, expr_dict):
    '''Function to return a new column and update the rows of a cid (e.g. band)'''
    tbl.add_index(cid)
    tbl['tmp_col'] = np.nan
    for key in tbl.group_by(cid).groups.keys:
        #print(key[cid])
        indx = tbl.loc_indices[key[cid]]
        #Astropy User Guide recommendation to modify Table using table[column][row] order
        target =  glom.glom(data_dict, expr_dict.format(key[cid]))
        #extract one element from a list
        [zp] = glom.flatten(target, levels=(expr_dict.count('**') - 1))
        tbl['tmp_col'][indx] = zp.value
    tbl['tmp_col'].unit = zp.unit
    new_col = tbl['tmp_col']
    tbl.remove_column('tmp_col')
    return new_col

def fluxToMag(flux):
    """ Return the magnitudes from flux quantities"""
    #return -2.5 * np.log10(flux.value)
    return u.Magnitude(flux)

def fluxErrToMag(flux, fluxerr):
    """ Return the magnitudes and associated errors from fluxes and flux error
    quantities. But doing the actual calculation in values"""
    mag = fluxToMag(flux)
    #return mag, -2.5 * np.log10( 1. - fluxerr / flux )
    return mag, u.Magnitude( 1. - fluxerr / flux ).value * mag.unit

def magToFlux(mag):
    """ Return the flux from magnitude quantities"""
    #return 10 ** (-0.4 * mag.value)
    return u.Dex(-0.4*mag.value).physical

def magErrToFlux(mag, err):
    """ Return the flux and associated error value from magnitude and mag error quantities"""
    flux = magToFlux(mag)
    return flux, flux * ( 1. - magToFlux(err) )

def data_convert_unit(t, flux_col, data_dict, sys, cid='band', zpt=None, target_unit=u.mJy, orig_unit=None, fluxe_col=None):
    '''Function to convert intrumental flux (or magnitude without physical type) to 
    calibrated physical type (using zeropoint)'''
    
    if t[flux_col].unit.is_equivalent(u.mag):
        t[flux_col].unit = u.mag()
        
    if fluxe_col in t.colnames:
        if t[fluxe_col].unit.is_equivalent(u.mag):
            t[fluxe_col].unit = u.mag()
    
    if not orig_unit:
        orig_unit = t[flux_col].unit
    
    if orig_unit.physical_type in ('unknown', 'dimensionless'):
        #check if field exists in dictionary
        if glom.glom(data_dict, '**.f_zp_nu'):
            f_zp_nu = column_from_dict_index(t, cid, data_dict, '**.{{}}.**.{0}.**.{1}'.format(sys, 'f_zp_nu'))
            if isinstance(orig_unit, u.MagUnit):
                if fluxe_col in t.colnames:
                    f, ef = f_zp_nu.quantity * magErrToFlux(t[flux_col].quantity, t[fluxe_col].quantity)# [u.Jy] * [dimensionless]
                else:
                    f = f_zp_nu.quantity * magToFlux(t[flux_col].quantity) 
            else:
                
                if glom.glom(data_dict, '**.zp'):
                    zp = column_from_dict_index(t, cid, data_dict, '**.{{}}.**.{0}.**.{1}'.format(sys, 'zp'))
                else:
                    print('No mag zeropoint exists in the data_dict')
                    return []
                if fluxe_col in t.colnames:
                    minst, eminst = fluxErrToMag(t[flux_col].quantity, t[fluxe_col].quantity)
                    if glom.glom(gaia, '**.e_zP'):
                        e_zp = column_from_dict_index(t, 'band', gaia, '**.{{}}.**.{0}.**.{1}'.format('VEGAMAG', 'e_zP'))
                        em = np.sqrt(eminst.value**2+ e_zp.value**2) * minst.unit # No astropy unit can do addition in quadrature
                    else:
                        em = eminst
                else:
                    minst = fluxToMag(t[flux_col].quantity)
                msys = minst + zp.quantity
                if fluxe_col in t.colnames:
                    f, ef = f_zp_nu.quantity * magErrToFlux(msys, em)
                else:
                    f = f_zp_nu.quantity * magToFlux(msys)
            if target_unit != f_zp_nu.unit:
                wave = column_from_dict_index(t, cid, data_dict, '**.{{}}.**.{0}.**.{1}'.format(sys, 'lamb'))
                vega = SourceSpectrum.from_vega()
                f = units.convert_flux(wave.quantity, f, target_unit, vegaspec=vega)
                if fluxe_col in t.colnames:
                    ef = units.convert_flux(wave.quantity, ef, target_unit, vegaspec=vega)               
        else:
            print('No flux zeropoint exists in the data_dict. We will do an approximation with convert_flux')
            #conversion
            wave = column_from_dict_index(t, cid, data_dict, '**.{{}}.**.{0}.**.{1}'.format(sys, 'lamb'))
            vega = SourceSpectrum.from_vega()  # For unit conversion  
            if isinstance(orig_unit, u.MagUnit):
                f = units.convert_flux(wave.quantity, t[flux_col].quantity, target_unit, vegaspec=vega)
                if fluxe_col in t.colnames:
                    msys_plus = t[flux_col].quantity + t[fluxe_col].quantity
                    f_plus = units.convert_flux(wave.quantity, msys_plus, target_unit, vegaspec=vega)
                    ef = abs(f_plus-f)
            else:
                if glom.glom(data_dict, '**.zp'):
                    zp = column_from_dict_index(t, cid, data_dict, '**.{{}}.**.{0}.**.{1}'.format(sys, 'zp'))
                else:
                    print('No mag zeropoint exists in the data_dict')
                    return []
                minst = fluxToMag(t[flux_col].quantity)
                msys = minst + zp.quantity
                f = units.convert_flux(wave.quantity, msys, target_unit, vegaspec=vega)
                if fluxe_col in t.colnames:
                    _, eminst = fluxErrToMag(t[flux_col].quantity, t[fluxe_col].quantity)
                    if glom.glom(gaia, '**.e_zP'):
                        e_zp = column_from_dict_index(t, 'band', gaia, '**.{{}}.**.{0}.**.{1}'.format('VEGAMAG', 'e_zP'))
                        em = np.sqrt(eminst.value**2+ e_zp.value**2) * minst.unit # No astropy unit can do addition in quadrature
                    else:
                        em = eminst
                    minst_plus = minst + em.value * u.mag()
                    msys_plus = minst_plus + zp.quantity
                    f_plus = units.convert_flux(wave.quantity, msys_plus, target_unit, vegaspec=vega)
                    ef = abs(f_plus-f)
                    
        if fluxe_col in t.colnames:
            return (f, ef)
        else:
            return (f)

@dataclass
class UnitHarmonizer:
    format_ref: str
    scale_ref: str
    ref_data: str
    ref_data_unit: str = 'FLUX'

#lists of times MaskedColumn and the sanitized astropy.table.Table
class DataProcess(object):
    '''
    Time series data aggregator and harmonizer class 
    list of `astropy.time.Time` and list of `astropy.table.QTable` per resource(votable) or extver(hdulist)
    '''
    def __init__(self, time_collection, table_collection, time_unit_ref=0, y_unit_ref=0):
        # Attribute definition for plotting data state
        self.time_collection = time_collection
        self.table_collection = table_collection
        self.time_unit_ref = time_unit_ref
        self.y_unit_ref = y_unit_ref
    def time_to_unit(time_list, colnames):
        '''Method to recursively convert to reference time format'''
    def y_column_list(table_collection, colnames):
        '''Method to select y_column and y_err_column in each element of a list'''
        y_data = [tbl.keep_columns(colnames) for tbl in table_collection]
        return y_data
    def to_unit(y_collection, err_y_collection, y_colname='flux', err_y_colname='flux_error'):
        '''Method to convert data column to reference unit equivalency
        Ref: jdaviz https://github.com/spacetelescope/jdaviz/blob/main/jdaviz/app.py UnitConverterWithSpectra class'''
        

    def group_by_index():
        '''Method to re-structure data by index ['band', 'instrument']'''

    
if __name__ == '__main__':

SYSTEM = 'VEGAMAG'
(f, ef) = data_convert_unit(table, 'flux', gaia, 'VEGAMAG', cid='band', target_unit=u.mJy, fluxe_col='flux_error')
