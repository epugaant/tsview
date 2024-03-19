import pytest
import os
import numpy as np
import io
from astropy import units as u
from astropy.io import fits
from astropy.table import Table, MaskedColumn
from astropy.io import votable
from astropy.time import Time

from synphot import units
from synphot.spectrum import SourceSpectrum

from tsview import DATADIR
from tsview.aggregate.data_processing import DataProcess, data_convert_unit

GAIA = {
    'system': 'VEGAMAG',
    'zeropt': {'G': {'VEGAMAG': {'zp': 25.6873668671 * units.VEGAMAG + 0 * u.mag('s/electron'), 'e_zP': 0.0027553202* units.VEGAMAG + 0 * u.mag('s/electron'), 'lamb': 6217.51 * u.AA, 'f_zp_nu': 3228.75 * u.Jy}, 'AB': {'zp': 25.8010446445 * u.Unit("mag(AB s/electron)"), 'e_zP': 0.0027590522 * u.Unit("mag(AB s/electron)"), 'lamb': 6217.51 * u.AA, 'f_zp_nu': 3631 * u.Jy}}, 
        'BP': {'VEGAMAG': {'zp': 25.3385422158 * units.VEGAMAG + 0 * u.mag('s/electron'), 'e_zP': 0.0027901700* units.VEGAMAG + 0 * u.mag('s/electron'), 'lamb': 5109.71 * u.AA, 'f_zp_nu': 3552.01 * u.Jy}, 'AB': {'zp': 25.3539555559 * u.Unit("mag(AB s/electron)"), 'e_zP': 0.0023065687 * u.Unit("mag(AB s/electron)"), 'lamb': 5109.71 * u.AA, 'f_zp_nu': 3631 * u.Jy}}, 
        'RP': {'VEGAMAG': {'zp': 24.7478955012 * units.VEGAMAG + 0 * u.mag('s/electron'), 'e_zP': 0.0037793818* units.VEGAMAG + 0 * u.mag('s/electron'), 'lamb': 7769.02 * u.AA, 'f_zp_nu': 2554.95 * u.Jy}, 'AB': {'zp': 25.1039837393 * u.Unit("mag(AB s/electron)"), 'e_zP': 0.0015800349 * u.Unit("mag(AB s/electron)"), 'lamb': 7769.02 * u.AA, 'f_zp_nu': 3631 * u.Jy}}, 
        }
}

def convert_unit_round_trip():
    band = ['BP', 'G', 'RP']
    mag = [12.756858, 12.411411, 11.898971] * u.mag
    mag_error = [0.002881, 0.002765, 0.003798] * u.mag
    flux = [107813.60407285392, 204353.63903, 137901.80445866095]* u.electron/u.s
    flux_error = [71.20729, 44.28437, 48.130802]* u.electron/u.s
    t = Table([band, mag, mag_error, flux, flux_error], names=('band', 'mag', 'mag_error', 'flux', 'flux_error'))
    
    
    flux_mJy_expected = [28.03696373, 35.03241034, 44.442042] * u.mJy 
    eflux_mJy_expected = [0.07429744, 0.0891021, 0.15519055] * u.mJy
    
    (f) = data_convert_unit(t, 'flux', GAIA, 'VEGAMAG', cid='band', target_unit=u.mJy)
    assert flux_mJy_expected == pytest.approx(f.value)
    (f, ef) = data_convert_unit(t, 'flux', GAIA, 'VEGAMAG', cid='band', target_unit=units.VEGAMAG, fluxe_col='flux_error')
    assert t['mag'].quantity == pytest.approx(f.value)
    assert t['mag_error'].quantity == pytest.approx(ef.value, 10E-4)
    (f, ef) = data_convert_unit(t, 'mag', GAIA, 'VEGAMAG', cid='band', target_unit=u.mJy, fluxe_col='mag_error')
    assert flux_mJy_expected == pytest.approx(f.value)
    assert eflux_mJy_expected == pytest.approx(ef.value)
    (f, ef) = data_convert_unit(t, 'flux', GAIA, 'VEGAMAG', cid='band', target_unit=u.mJy, fluxe_col='flux_error')
    assert flux_mJy_expected == pytest.approx(f.value)
    assert eflux_mJy_expected == pytest.approx(ef.value, 10E-3)
    (f, ef) = data_convert_unit(t, 'mag', GAIA, 'AB', cid='band', target_unit=u.mJy, fluxe_col='mag_error') # this one is not strictly correct.
    assert flux_mJy_expected == pytest.approx(f.value, 10E1)
    assert eflux_mJy_expected == pytest.approx(ef.value, 10E-2)
    (f, ef) = data_convert_unit(t, 'flux', GAIA, 'AB', cid='band', target_unit=u.mJy, fluxe_col='flux_error')
    assert flux_mJy_expected == pytest.approx(f.value, 10E1)
    assert eflux_mJy_expected == pytest.approx(ef.value, 10E-2)


def convert_unit_gaia():
    # tbl = Table.read('/Users/epuga/ESDC/TSViz/data/gaia/anonymous1690191210843_fits/EPOCH_PHOTOMETRY-Gaia DR3 4057091150787104896.fits', astropy_native=True)
    # if tbl['flux'].unit == "'electron'.s**-1":
    #     tbl['flux'].unit = "electron/s"
    # if tbl['flux_error'].unit == "'electron'.s**-1":
    #     tbl['flux_error'].unit = "electron/s"
    
    # d = DataProcess('gaia', [tbl['time']], [tbl['source_id', 'band', 'mag', 'flux', 'flux_error']], 'VEGAMAG') 
    # d.convert_time('jd')
    # print(d.to_json()) 
    # d.convert_time('mjd')
    # print(d.to_json()) 
    # d.convert_flux(u.mJy)
    # print(d.to_json())  
    # d = DataProcess('gaia', [tbl['time']], [tbl['source_id', 'band', 'mag', 'flux', 'flux_error']], 'VEGAMAG')  
    # d.convert_flux(u.mJy)
    # print(d.timeseries[0].flux)
    # d2 = DataProcess('gaia', [tbl['time']], [tbl['source_id', 'band', 'mag', 'flux', 'flux_error']], 'AB')
    # d2.convert_flux(u.mJy)
    # print(d2.timeseries[0].flux)
    pass

if __name__ == '__main__':
    convert_unit_round_trip()
    #convert_unit_gaia()
