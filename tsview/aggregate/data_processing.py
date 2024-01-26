from dataclasses import dataclass, field


import glom
import numpy as np

from astropy.time import Time
from astropy.table import Table, QTable
import astropy.units as u
from astropy.units import Quantity
from synphot import units
from synphot.spectrum import SourceSpectrum


DATA_DICT = {
'gaia': {
    'zeropt': {'G': {'VEGAMAG': {'zp': 25.6873668671 * units.VEGAMAG + 0 * u.mag('s/electron'), 'e_zP': 0.0027553202* units.VEGAMAG + 0 * u.mag('s/electron'), 'lamb': 6217.51 * u.AA, 'f_zp_nu': 3228.75 * u.Jy}, 'AB': {'zp': 25.8010446445 * u.Unit("mag(AB s/electron)"), 'e_zP': 0.0027590522 * u.Unit("mag(AB s/electron)"), 'lamb': 6217.51 * u.AA, 'f_zp_nu': 3631 * u.Jy}}, 
        'BP': {'VEGAMAG': {'zp': 25.3385422158 * units.VEGAMAG + 0 * u.mag('s/electron'), 'e_zP': 0.0027901700* units.VEGAMAG + 0 * u.mag('s/electron'), 'lamb': 5109.71 * u.AA, 'f_zp_nu': 3552.01 * u.Jy}, 'AB': {'zp': 25.3539555559 * u.Unit("mag(AB s/electron)"), 'e_zP': 0.0023065687 * u.Unit("mag(AB s/electron)"), 'lamb': 5109.71 * u.AA, 'f_zp_nu': 3631 * u.Jy}}, 
        'RP': {'VEGAMAG': {'zp': 24.7478955012 * units.VEGAMAG + 0 * u.mag('s/electron'), 'e_zP': 0.0037793818* units.VEGAMAG + 0 * u.mag('s/electron'), 'lamb': 7769.02 * u.AA, 'f_zp_nu': 2554.95 * u.Jy}, 'AB': {'zp': 25.1039837393 * u.Unit("mag(AB s/electron)"), 'e_zP': 0.0015800349 * u.Unit("mag(AB s/electron)"), 'lamb': 7769.02 * u.AA, 'f_zp_nu': 3631 * u.Jy}}, 
        },
    'graphic': {'y': {'colname': 'flux', 'units': 'electron/s'},
              'y_err': {'colname': 'flux_error', 'units': 'electron/s'},
              'cid': 'band' }    
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
class TimeSeries:
    time: Quantity
    flux: Quantity
    flux_error: Quantity|None = field(intit=False)
    time_format: str
    data_unit: u
    cid:  str|None = field(intit=False)
    
    def __repr__(self):
        newt = QTable([self.time, self.flux, self.flux_error], names=['time', 'flux', 'flux_error'])
        if self.cid:
            newt.add_column(t['band'])

#lists of Times and the sanitized astropy.table.Table
@dataclass
class DataProcess:
    mission: str
    time_collection: list[Time]
    table_collection: list[Table]
    #time_ref_format: str = field(init=False)
    #time_ref_scale: str = field(init=False)
    system: str
    time_target_format: str
    #data_orig_unit: u = field(init=False)
    data_target_unit: u = field(default= u.mJy)
    cid: str|None = field(init=False)
    y_colname: str = field(init=False)
    err_y_colname: str|None = field(init=False)
    timeseries: list[TimeSeries] = field(default_factory=list, init=False)
    
    def __post_init__(self):
        graphic_dict =  DATA_DICT[self.mission]['graphic']
        [self.cid] = glom.glom(graphic_dict, '**.cid')
        [self.y_colname] = glom.glom(graphic_dict, '**.y.colname')
        if glom.glom(graphic_dict, '**.y_err.colname'):
            [self.err_y_colname] = glom.glom(graphic_dict, '**.y_err.colname')
        if self.err_y_colname:
            self.timeseries = [TimeSeries(time, table[self.y_colname], table[self.err_y_colname], time.format, table[self.y_colname].unit, self.cid) for time, table in zip(self.time_collection, self.table_collection)]
        else:
            self.timeseries = [TimeSeries(time, table[self.y_colname], time.format, table[self.y_colname].unit, self.cid) for time, table in zip(self.time_collection, self.table_collection)]

     
    def convert_time(self, target_unit):
        for i in range(self.timeseries):
            if self.timeseries[i].time_format != target_unit:
                self.timeseries[i].time_format = target_unit

    def convert_flux(self, target_unit):
        zpt_dict = DATA_DICT[self.mission]['zeropt']
        for i in range(self.timeseries):
            if self.timeseries[i].data_unit != target_unit:
                if self.err_y_colname:
                    self.timeseries[i].flux, self.timeseries[i].flux_error = self.data_convert_unit(self.table_collection[i], self.y_colname, zpt_dict, self.system, cid=self.cid, target_unit=target_unit, fluxe_col=self.err_y_colname)
                else:
                    self.timeseries[i].flux = self.data_convert_unit(self.table_collection[i], self.y_colname, zpt_dict, self.system, cid=self.cid, target_unit=target_unit)

    def group_by_index(cid):
        '''Method to re-structure data by index ['band', 'instrument']'''

    
if __name__ == '__main__':

    tbl = Table.read('/Users/epuga/ESDC/TSViz/data/gaia/anonymous1690191210843_fits/EPOCH_PHOTOMETRY-Gaia DR3 4057091150787104896.fits', astropy_native=True)
    df = tbl.to_pandas()
    #thanks to https://stackoverflow.com/questions/22219004/how-to-group-dataframe-rows-into-list-in-pandas-groupby

    d = DataProcess(time, tbl, system, target_time_unit, target_flux_unit, cid)
    newt= QTable([time, f, ef], names=['time', 'flux', 'flux_error'])
    newt= QTable([time.jd, f, ef], names=['time', 'flux', 'flux_error'])
    newt.add_column(t['band'])
    df = newt.to_pandas()
    df2 = df.groupby('band')[['time', 'flux', 'flux_error']].agg(list).rename(columns={'time': 'x', 'flux': 'y', 'flux_error': 'error_y'}).to_json(orient='index')