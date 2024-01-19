'''27.12.2023 Predecessor was test_table.py to get the errors right using Vizier sig_m calculation and test using the two use-cases'''

from astropy.table import Table, Column
import numpy as np
import astropy.units as u
from synphot import units
import glom
from synphot.spectrum import SourceSpectrum

gaia = {'zeropt': {'G': {'VEGAMAG': {'zp': 25.6873668671 * units.VEGAMAG + 0 * u.mag('s/electron'), 'e_zP': 0.0027553202* units.VEGAMAG + 0 * u.mag('s/electron'), 'lamb': 6217.51 * u.AA, 'f_zp_nu': 3228.75 * u.Jy}, 'AB': {'zp': 25.8010446445 * u.Unit("mag(AB s/electron)"), 'e_zP': 0.0027590522 * u.Unit("mag(AB s/electron)"), 'lamb': 6217.51 * u.AA, 'f_zp_nu': 3631 * u.Jy}}, 
        'BP': {'VEGAMAG': {'zp': 25.3385422158 * units.VEGAMAG + 0 * u.mag('s/electron'), 'e_zP': 0.0027901700* units.VEGAMAG + 0 * u.mag('s/electron'), 'lamb': 5109.71 * u.AA, 'f_zp_nu': 3552.01 * u.Jy}, 'AB': {'zp': 25.3539555559 * u.Unit("mag(AB s/electron)"), 'e_zP': 0.0023065687 * u.Unit("mag(AB s/electron)"), 'lamb': 5109.71 * u.AA, 'f_zp_nu': 3631 * u.Jy}}, 
        'RP': {'VEGAMAG': {'zp': 24.7478955012 * units.VEGAMAG + 0 * u.mag('s/electron'), 'e_zP': 0.0037793818* units.VEGAMAG + 0 * u.mag('s/electron'), 'lamb': 7769.02 * u.AA, 'f_zp_nu': 2554.95 * u.Jy}, 'AB': {'zp': 25.1039837393 * u.Unit("mag(AB s/electron)"), 'e_zP': 0.0015800349 * u.Unit("mag(AB s/electron)"), 'lamb': 7769.02 * u.AA, 'f_zp_nu': 3631 * u.Jy}}, 
        }
}
gaio = {'zeropt': {'G': {'VEGAMAG': {'zp': 25.6873668671 * units.VEGAMAG + 0 * u.mag('s/electron'), 'e_zP': 0.0027553202* units.VEGAMAG + 0 * u.mag('s/electron'), 'lamb': 6217.51 * u.AA}, 'AB': {'zp': 25.8010446445 * u.Unit("mag(AB s/electron)"), 'e_zP': 0.0027590522 * u.Unit("mag(AB s/electron)"), 'lamb': 6217.51 * u.AA}}, 
        'BP': {'VEGAMAG': {'zp': 25.3385422158 * units.VEGAMAG + 0 * u.mag('s/electron'), 'e_zP': 0.0027901700* units.VEGAMAG + 0 * u.mag('s/electron'), 'lamb': 5109.71 * u.AA}, 'AB': {'zp': 25.3539555559 * u.Unit("mag(AB s/electron)"), 'e_zP': 0.0023065687 * u.Unit("mag(AB s/electron)"), 'lamb': 5109.71 * u.AA}}, 
        'RP': {'VEGAMAG': {'zp': 24.7478955012 * units.VEGAMAG + 0 * u.mag('s/electron'), 'e_zP': 0.0037793818* units.VEGAMAG + 0 * u.mag('s/electron'), 'lamb': 7769.02 * u.AA}, 'AB': {'zp': 25.1039837393 * u.Unit("mag(AB s/electron)"), 'e_zP': 0.0015800349 * u.Unit("mag(AB s/electron)"), 'lamb': 7769.02 * u.AA}}, 
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
                    #f, ef = magErrToFlux(t[flux_col].quantity, t[fluxe_col].quantity) * zp.quantity # [dimensionless] * [u.Jy]
                    f, ef = f_zp_nu.quantity * magErrToFlux(t[flux_col].quantity, t[fluxe_col].quantity)# [u.Jy] * [dimensionless]
                else:
                    #f = magToFlux(t[flux_col].quantity) * zp.quantity
                    f = f_zp_nu.quantity * magToFlux(t[flux_col].quantity) 
            else:
                
                if glom.glom(data_dict, '**.zp'):
                    zp = column_from_dict_index(t, cid, data_dict, '**.{{}}.**.{0}.**.{1}'.format(sys, 'zp'))
                else:
                    print('No mag zeropoint exists in the data_dict')
                    return []
                #instrumental flux
                # if fluxe_col in t.colnames:
                #     m, em = fluxErrToMag(t[flux_col].quantity, t[fluxe_col].quantity) 
                # else:
                #     m = fluxToMag(t[flux_col].quantity)
                # m = m + zp.quantity
                # if fluxe_col in t.colnames:
                #     if glom.glom(gaia, '**.e_zP'):
                #         e_zp = column_from_dict_index(t, cid, data_dict, '**.{{}}.**.{0}.**.{1}'.format(sys, 'e_zP'))
                #         em = np.sqrt(em.value**2+ e_zp.value**2)  # No astropy unit can do addition in quadrature
                # # slightly different input in both cases
                # if fluxe_col in t.colnames:
                #     f, ef = m, m + (em * u.mag)
                # else:
                #     f = m
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
        #     return []
            #conversion
            wave = column_from_dict_index(t, cid, data_dict, '**.{{}}.**.{0}.**.{1}'.format(sys, 'lamb'))
            vega = SourceSpectrum.from_vega()  # For unit conversion  
        # f_final = units.convert_flux(wave.quantity, f, target_unit, vegaspec=vega)
        # if fluxe_col in t.colnames:
        #     ef_final = units.convert_flux(wave.quantity, ef, target_unit, vegaspec=vega)
        #     print(f_final, ef_final)
        #     return (f_final, ef_final)
        # else:
        #     print(f_final)
        #     return (f_final)
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

def main():
    band = ['BP', 'G', 'RP']
    mag = [12.756858, 12.411411, 11.898971] * u.mag
    mag_error = [0.002881, 0.002765, 0.003798] * u.mag
    flux = [107813.60407285392, 204353.63903, 137901.80445866095]* u.electron/u.s
    flux_error = [71.20729, 44.28437, 48.130802]* u.electron/u.s
    t = Table([band, mag, mag_error, flux, flux_error], names=('band', 'mag', 'mag_error', 'flux', 'flux_error'))
    print(t)
    
    (f) = data_convert_unit(t, 'flux', gaia, 'VEGAMAG', cid='band', target_unit=u.mJy)
    print('data_convert_unit 1 : {}'.format(f))
    (f, ef) = data_convert_unit(t, 'mag', gaia, 'VEGAMAG', cid='band', target_unit=u.mJy, fluxe_col='mag_error')
    print('data_convert_unit 2 : {} {}'.format(f, ef))
    (f, ef) = data_convert_unit(t, 'flux', gaia, 'VEGAMAG', cid='band', target_unit=u.mJy, fluxe_col='flux_error')
    print('data_convert_unit 3 : {} {}'.format(f, ef))
    (f, ef) = data_convert_unit(t, 'mag', gaia, 'AB', cid='band', target_unit=u.mJy, fluxe_col='mag_error') # this one is not strictly correct.
    print('data_convert_unit 4 : {} {}'.format(f, ef))
    (f, ef) = data_convert_unit(t, 'flux', gaia, 'AB', cid='band', target_unit=u.mJy, fluxe_col='flux_error')
    print('data_convert_unit 5 : {} {}'.format(f, ef))

    (f) = data_convert_unit(t, 'flux', gaio, 'VEGAMAG', cid='band', target_unit=u.mJy)
    print('data_convert_unit 6 : {}'.format(f))
    t['mag'].unit = units.VEGAMAG
    #t['mag_error'].unit = units.VEGAMAG
    (f, ef) = data_convert_unit(t, 'mag', gaio, 'VEGAMAG', cid='band', target_unit=u.mJy, fluxe_col='mag_error')
    print('data_convert_unit 7 : {} {}'.format(f, ef))
    (f, ef) = data_convert_unit(t, 'flux', gaio, 'VEGAMAG', cid='band', target_unit=u.mJy, fluxe_col='flux_error')
    print('data_convert_unit 8 : {} {}'.format(f, ef))

    
    if t['mag'].unit.is_equivalent(u.mag):
        t['mag'].unit = u.mag()
        
    if t['mag_error'].unit.is_equivalent(u.mag):
        t['mag_error'].unit = u.mag()
        
    wave = column_from_dict_index(t, 'band', gaia, '**.{{}}.**.{0}.**.{1}'.format('VEGAMAG', 'lamb'))
    vega = SourceSpectrum.from_vega()  # For unit conversion  
    
    f_zp_nu = column_from_dict_index(t, 'band', gaia, '**.{{}}.**.{0}.**.{1}'.format('VEGAMAG', 'f_zp_nu'))
    zp = column_from_dict_index(t, 'band', gaia, '**.{{}}.**.{0}.**.{1}'.format('VEGAMAG', 'zp'))
    
    #Sannity check --[e/s] <-> mag
    print('flux [e/s] to mag')
    print(fluxToMag(t['flux'].quantity) + zp.quantity)
    print(t['mag'].quantity)

    print('mag to flux [e/s]')
    print(magToFlux(t['mag'].quantity - zp.quantity)) # minst
    print(t['flux'].quantity)
    
    #Conversion -- mag(VEGA)/f[e/s] -> Jy with fnu zeropoints 
    f_zp_nu.quantity * magToFlux(t['mag'].quantity)
    f_zp_nu.quantity * magToFlux(fluxToMag(t['flux'].quantity) + zp.quantity)
    
    # Do you have fnu zeropoints? --> block 1/2
    # do you start from magnitudes or flux? if flux, you need zp (check if you have e_zp)
    
    #1 Conversion with errors -- mag(VEGA)/f[e/s] -> Jy with fnu zeropoints
    f, ef = f_zp_nu.quantity * magErrToFlux(t['mag'].quantity, t['mag_error'].quantity) # already in Jy
    
    minst, eminst = fluxErrToMag(t['flux'].quantity, t['flux_error'].quantity)
    e_zp = column_from_dict_index(t, 'band', gaia, '**.{{}}.**.{0}.**.{1}'.format('VEGAMAG', 'e_zP'))
    if glom.glom(gaia, '**.e_zP'):
        em = np.sqrt(eminst.value**2+ e_zp.value**2) * minst.unit # No astropy unit can do addition in quadrature
    else:
        em = eminst
    msys = minst + zp.quantity
    f, ef = f_zp_nu.quantity * magErrToFlux(msys, em)
    
    #2 Conversion with errors -- mag(VEGA)/(f[e/s] + zp) -> Jy with convert_flux (dont have zeropoints) It is not equivalent
    f = units.convert_flux(wave.quantity, (t['mag'].quantity).value*units.VEGAMAG, u.Jy, vegaspec=vega)
    msys_plus = t['mag'].quantity + t['mag_error'].quantity
    f_plus = units.convert_flux(wave.quantity, msys_plus.value*units.VEGAMAG, u.Jy, vegaspec=vega)
    ef = f_plus-f
    print(f, ef)

    minst = fluxToMag(t['flux'].quantity)
    msys = minst + zp.quantity
    f = units.convert_flux(wave.quantity, msys, u.Jy, vegaspec=vega)
    minst_plus = fluxToMag(t['flux'].quantity+t['flux_error'].quantity)
    msys_plus = minst_plus + zp.quantity
    f_plus = units.convert_flux(wave.quantity, msys_plus, u.Jy, vegaspec=vega)
    ef = f_plus-f
    print(f, ef)
    
    ##########################
    print('')
    print('m, sig_m (without physical type) to f, sig_f (mJy), going through f_nu zeropoint')
    
    #f = magToFlux(t['mag'].quantity) * gaia['zeropt']['BP']['VEGAMAG']['f_zp_nu']
    minst = t['mag'].quantity - zp.quantity
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
    print('')
    
    t = Table([band, mag, mag_error, flux, flux_error], names=('band', 'mag', 'mag_error', 'flux', 'flux_error'))

    (f, ef) = data_convert_unit(t, 'mag', gaia, 'VEGAMAG', cid='band', zpt='f_zp_nu', target_unit=u.mJy, fluxe_col='mag_error')
    print('data_convert_unit 1 : {} {}'.format(f, ef))
    (f) = data_convert_unit(t, 'flux', gaia, 'VEGAMAG', cid='band', zpt='zp', target_unit=u.mJy)
    print('data_convert_unit 2 : {}'.format(f))
    (f, ef) = data_convert_unit(t, 'flux', gaia, 'VEGAMAG', cid='band', zpt='zp', target_unit=u.mJy, fluxe_col='flux_error')
    print('data_convert_unit 3 : {} {}'.format(f, f - ef))
    pass

if __name__ == '__main__':
    main()
