'''27.12.2023 Predecessor was test_table.py to get the errors right using Vizier sig_m calculation and test using the two use-cases'''

from astropy.table import Table, Column
import numpy as np
import astropy.units as u
from synphot import units
import glom
from synphot.spectrum import SourceSpectrum

gaia = {'zeropt': {'G': {'VEGAMAG': {'zp': 25.6873668671 * units.VEGAMAG + 0 * u.mag('s/electron'), 'e_zP': 0.0027553202* units.VEGAMAG + 0 * u.mag('s/electron'), 'lamb': 6251.50 * u.AA, 'f_zp_nu': 3599.12 * u.Jy}, 'AB': {'zp': 25.8010446445 * u.Unit("mag(AB s/electron)"), 'e_zP': 0.0027590522 * u.Unit("mag(AB s/electron)"), 'lamb': 6251.50 * u.AA, 'f_zp_nu': 3631 * u.Jy}}, 
        'BP': {'VEGAMAG': {'zp': 25.3385422158 * units.VEGAMAG + 0 * u.mag('s/electron'), 'e_zP': 0.0027901700* units.VEGAMAG + 0 * u.mag('s/electron'), 'lamb': 5124.20 * u.AA, 'f_zp_nu': 3256.01 * u.Jy}, 'AB': {'zp': 25.3539555559 * u.Unit("mag(AB s/electron)"), 'e_zP': 0.0023065687 * u.Unit("mag(AB s/electron)"), 'lamb': 5124.20 * u.AA, 'f_zp_nu': 3631 * u.Jy}}, 
        'RP': {'VEGAMAG': {'zp': 24.7478955012 * units.VEGAMAG + 0 * u.mag('s/electron'), 'e_zP': 0.0037793818* units.VEGAMAG + 0 * u.mag('s/electron'), 'lamb': 7829.65 * u.AA, 'f_zp_nu': 2530.35 * u.Jy}, 'AB': {'zp': 25.1039837393 * u.Unit("mag(AB s/electron)"), 'e_zP': 0.0015800349 * u.Unit("mag(AB s/electron)"), 'lamb': 7829.65 * u.AA, 'f_zp_nu': 3631 * u.Jy}}, 
        }
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
    """ Return the flux and associated errors from magnitude and mag error quantities"""
    flux = magToFlux(mag)
    return flux, flux * ( 1. - magToFlux(err) )

def data_convert_unit(t, flux_col, data_dict, sys, cid='band', zpt=None, target_unit=u.mJy, orig_unit=None, fluxe_col=None):
    '''Function to convert intrumental flux (or magnitude without physical type) to 
    calibrated physical type (using zeropoint)'''
    if not orig_unit:
        orig_unit = t[flux_col].unit
    
    if t[flux_col].unit.is_equivalent(u.mag):
        t[flux_col].unit = u.mag()
        
    if t[fluxe_col].unit.is_equivalent(u.mag):
        t[fluxe_col].unit = u.mag()
    
    if orig_unit.physical_type in ('unknown', 'dimensionless'):
        if zpt:
            #check if field exists in dictionary
            if glom.glom(data_dict, '**.'+zpt):
                zp = column_from_dict_index(t, cid, data_dict, '**.{{}}.**.{0}.**.{1}'.format(sys, zpt))
                
                if isinstance(orig_unit, u.MagUnit):
                    m, em = t[flux_col].quantity, t[fluxe_col].quantity
                else:
                    #instrumental flux
                    if fluxe_col in t.colnames:
                        m, em = fluxErrToMag(t[flux_col].quantity, t[fluxe_col].quantity) 
                    else:
                        m = fluxToMag(t[flux_col].quantity)
                    m = m + zp.quantity
                    if glom.glom(gaia, '**.e_zP'):
                        e_zp = column_from_dict_index(t, cid, data_dict, '**.{{}}.**.{0}.**.{1}'.format(sys, 'e_zP'))
                        em = np.sqrt(em.value**2+ e_zp.value**2) * m.unit # No astropy unit can do addition in quadrature
                    else:
                        em = em + zp.quantity
                # common to both cases
                if fluxe_col in t.colnames:
                    f, ef = magErrToFlux(m.quantity, em.quantity) * zp.quantity
                else:
                    f = magToFlux(m.quantity) * zp.quantity
            else:
                print('No zeropoint exists in the data_dict')
                return []
        else:
            print('No zeropoint given for a non physical unit')
            return []
        #conversion block
        wave = column_from_dict_index(t, cid, data_dict, '**.{{}}.**.{0}.**.{1}'.format(sys, 'lamb'))
        vega = SourceSpectrum.from_vega()  # For unit conversion  
        f_final = units.convert_flux(wave.quantity, f, target_unit, vegaspec=vega)
        if fluxe_col in t.colnames:
            ef_final = units.convert_flux(wave.quantity, ef, target_unit, vegaspec=vega)
            print(f_final, ef_final)
            return (f_final, ef_final)
        else:
            print(f_final)
            return (f_final)
    

def main():
    band = ['BP', 'G', 'RP']
    mag = [12.756858, 12.411411, 11.898971] * u.mag
    mag_error = [0.002881, 0.002765, 0.003798] * u.mag
    flux = [107813.60407285392, 204353.63903, 137901.80445866095]* u.electron/u.s
    flux_error = [71.20729, 44.28437, 48.130802]* u.electron/u.s
    t = Table([band, mag, mag_error, flux, flux_error], names=('band', 'mag', 'mag_error', 'flux', 'flux_error'))
    print(t)

    if t['mag'].unit.is_equivalent(u.mag):
        t['mag'].unit = u.mag()
        
    if t['mag_error'].unit.is_equivalent(u.mag):
        t['mag_error'].unit = u.mag()
        
    wave = column_from_dict_index(t, 'band', gaia, '**.{{}}.**.{0}.**.{1}'.format('VEGAMAG', 'lamb'))
    vega = SourceSpectrum.from_vega()  # For unit conversion  
    
    print('')
    print('m, sig_m (without physical type) to f, sig_f (mJy), going through f_nu zeropoint')
    f_zp_nu = column_from_dict_index(t, 'band', gaia, '**.{{}}.**.{0}.**.{1}'.format('VEGAMAG', 'f_zp_nu'))
    #f = magToFlux(t['mag'].quantity) * gaia['zeropt']['BP']['VEGAMAG']['f_zp_nu']
    f, ef = magErrToFlux(t['mag'].quantity, t['mag_error'].quantity) * f_zp_nu.quantity
    #already in Jy
    f_final = units.convert_flux(wave.quantity, f, u.mJy, vegaspec=vega)
    ef_final = units.convert_flux(wave.quantity, ef, u.mJy, vegaspec=vega)
    print(f_final, ef_final)
    
    print('')
    print('f, sig_f (instrumental flux) to m, sig_m (VEGAMAG), going through zp zeropoint, and then converting to Jy')
    # m = fluxToMag(t['flux'].quantity) + gaia['zeropt']['BP']['VEGAMAG']['zp']
    # m_inter = units.convert_flux(wave.quantity, m, u.mJy, vegaspec=vega)
    m, em = fluxErrToMag(t['flux'].quantity, t['flux_error'].quantity) 
    zp = column_from_dict_index(t, 'band', gaia, '**.{{}}.**.{0}.**.{1}'.format('VEGAMAG', 'zp'))
    e_zp = column_from_dict_index(t, 'band', gaia, '**.{{}}.**.{0}.**.{1}'.format('VEGAMAG', 'e_zP'))
    m = m + zp.quantity
    if glom.glom(gaia, '**.e_zP'):
        em = np.sqrt(em.value**2+ e_zp.value**2) * m.unit # No astropy unit can do addition in quadrature
    else:
        em = em + zp.quantity
    # Now from apparent magnitude mag(VEGA) to flux
    ff, eff = magErrToFlux(m, em) * f_zp_nu.quantity
    m_final =  units.convert_flux(wave.quantity, ff, u.mJy, vegaspec=vega)
    em_final = units.convert_flux(wave.quantity, eff, u.mJy, vegaspec=vega)
    print(m_final, em_final)
    pass

if __name__ == '__main__':
    main()
