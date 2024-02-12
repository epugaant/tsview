from dataclasses import dataclass, field


import glom
import numpy as np

from astropy.time import Time
from astropy.table import Table, QTable
import astropy.units as u
from astropy.units import Quantity
from synphot import units
from synphot.spectrum import SourceSpectrum

import pandas as pd
import json


DATA_DICT = {
'gaia': {
    'system': 'AB',
    'zeropt': {'G': {'VEGAMAG': {'zp': 25.6873668671 * units.VEGAMAG + 0 * u.mag('s/electron'), 'e_zP': 0.0027553202* units.VEGAMAG + 0 * u.mag('s/electron'), 'lamb': 6217.51 * u.AA, 'f_zp_nu': 3228.75 * u.Jy}, 'AB': {'zp': 25.8010446445 * u.Unit("mag(AB s/electron)"), 'e_zP': 0.0027590522 * u.Unit("mag(AB s/electron)"), 'lamb': 6217.51 * u.AA, 'f_zp_nu': 3631 * u.Jy}}, 
        'BP': {'VEGAMAG': {'zp': 25.3385422158 * units.VEGAMAG + 0 * u.mag('s/electron'), 'e_zP': 0.0027901700* units.VEGAMAG + 0 * u.mag('s/electron'), 'lamb': 5109.71 * u.AA, 'f_zp_nu': 3552.01 * u.Jy}, 'AB': {'zp': 25.3539555559 * u.Unit("mag(AB s/electron)"), 'e_zP': 0.0023065687 * u.Unit("mag(AB s/electron)"), 'lamb': 5109.71 * u.AA, 'f_zp_nu': 3631 * u.Jy}}, 
        'RP': {'VEGAMAG': {'zp': 24.7478955012 * units.VEGAMAG + 0 * u.mag('s/electron'), 'e_zP': 0.0037793818* units.VEGAMAG + 0 * u.mag('s/electron'), 'lamb': 7769.02 * u.AA, 'f_zp_nu': 2554.95 * u.Jy}, 'AB': {'zp': 25.1039837393 * u.Unit("mag(AB s/electron)"), 'e_zP': 0.0015800349 * u.Unit("mag(AB s/electron)"), 'lamb': 7769.02 * u.AA, 'f_zp_nu': 3631 * u.Jy}}, 
        },
    'graphic': {'y': {'colname': 'flux', 'units': 'electron/s'},
              'y_err': {'colname': 'flux_error', 'units': 'electron/s'},
              'cid': 'band' }    
        },
'jwst': {
    'system': 'AB',
    'zeropt': {'G': {'VEGAMAG': {'zp': 25.6873668671 * units.VEGAMAG + 0 * u.mag('s/electron'), 'e_zP': 0.0027553202* units.VEGAMAG + 0 * u.mag('s/electron'), 'lamb': 6217.51 * u.AA, 'f_zp_nu': 3228.75 * u.Jy}, 'AB': {'zp': 25.8010446445 * u.Unit("mag(AB s/electron)"), 'e_zP': 0.0027590522 * u.Unit("mag(AB s/electron)"), 'lamb': 6217.51 * u.AA, 'f_zp_nu': 3631 * u.Jy}}, 
        'BP': {'VEGAMAG': {'zp': 25.3385422158 * units.VEGAMAG + 0 * u.mag('s/electron'), 'e_zP': 0.0027901700* units.VEGAMAG + 0 * u.mag('s/electron'), 'lamb': 5109.71 * u.AA, 'f_zp_nu': 3552.01 * u.Jy}, 'AB': {'zp': 25.3539555559 * u.Unit("mag(AB s/electron)"), 'e_zP': 0.0023065687 * u.Unit("mag(AB s/electron)"), 'lamb': 5109.71 * u.AA, 'f_zp_nu': 3631 * u.Jy}}, 
        'RP': {'VEGAMAG': {'zp': 24.7478955012 * units.VEGAMAG + 0 * u.mag('s/electron'), 'e_zP': 0.0037793818* units.VEGAMAG + 0 * u.mag('s/electron'), 'lamb': 7769.02 * u.AA, 'f_zp_nu': 2554.95 * u.Jy}, 'AB': {'zp': 25.1039837393 * u.Unit("mag(AB s/electron)"), 'e_zP': 0.0015800349 * u.Unit("mag(AB s/electron)"), 'lamb': 7769.02 * u.AA, 'f_zp_nu': 3631 * u.Jy}}, 
        },
    'graphic': {'y': {'colname': 'FLUX', 'units': 'Jy'},
              'y_err': {'colname': 'FLUX_ERROR', 'units': 'Jy'},
              'cid': None,
              'cextra': 'WAVELENGTH',
              'multip': 'FILTER'}    
        }
}

def column_from_column_index(tbl, cid, data_dict, expr_dict):
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

def column_from_index(tbl, index, data_dict, expr_dict):
    target =  glom.glom(data_dict, expr_dict.format(index))
    [val] = glom.flatten(target, levels=(expr_dict.count('**') - 1))
    tbl['tmp_col'] = val
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

def data_convert_unit(t, flux_col, data_dict, sys, cid=None, id=None, target_unit=None, orig_unit=None, fluxe_col=None):
    '''Function to convert intrumental flux (or magnitude without physical type) to 
    calibrated physical type (using zeropoint)'''
    
    if t[flux_col].unit.is_equivalent(u.mag):
        t[flux_col].unit = u.mag()
        
    if fluxe_col in t.colnames:
        if t[fluxe_col].unit.is_equivalent(u.mag):
            t[fluxe_col].unit = u.mag()
            
    if not target_unit:
        targer_unit = u.mJy
        
    if not orig_unit:
        orig_unit = t[flux_col].unit
    
    if orig_unit.physical_type in ('unknown', 'dimensionless'):
        #check if field exists in dictionary
        if glom.glom(data_dict, '**.f_zp_nu'):
            if cid is not None and id is not None:
                raise ValueError('Should have cid, or id, but not both')
            elif cid:
                f_zp_nu = column_from_column_index(t, cid, data_dict, '**.{{}}.**.{0}.**.{1}'.format(sys, 'f_zp_nu'))
            elif id:
                f_zp_nu = column_from_index(t, id, data_dict, '**.{{}}.**.{0}.**.{1}'.format(sys, 'f_zp_nu'))
            if isinstance(orig_unit, u.MagUnit):
                if fluxe_col in t.colnames:
                    f, ef = f_zp_nu.quantity * magErrToFlux(t[flux_col].quantity, t[fluxe_col].quantity)# [u.Jy] * [dimensionless]
                else:
                    f = f_zp_nu.quantity * magToFlux(t[flux_col].quantity) 
            else:
                
                if glom.glom(data_dict, '**.zp'):
                    if cid is not None and id is not None:
                        raise ValueError('Should have cid, or id, but not both')
                    elif cid:
                        zp = column_from_column_index(t, cid, data_dict, '**.{{}}.**.{0}.**.{1}'.format(sys, 'zp'))
                    elif id:
                        zp = column_from_index(t, id, data_dict, '**.{{}}.**.{0}.**.{1}'.format(sys, 'zp'))
                else:
                    print('No mag zeropoint exists in the data_dict')
                    return []
                if fluxe_col in t.colnames:
                    minst, eminst = fluxErrToMag(t[flux_col].quantity, t[fluxe_col].quantity)
                    if glom.glom(data_dict, '**.e_zP'):
                        e_zp = column_from_column_index(t, 'band', data_dict, '**.{{}}.**.{0}.**.{1}'.format(sys, 'e_zP'))
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
                if cid is not None and id is not None:
                    raise ValueError('Should have cid, or id, but not both')
                elif cid:
                    wave = column_from_column_index(t, cid, data_dict, '**.{{}}.**.{0}.**.{1}'.format(sys, 'lamb'))
                elif id:
                    wave = column_from_index(t, id, data_dict, '**.{{}}.**.{0}.**.{1}'.format(sys, 'lamb')) 
                vega = SourceSpectrum.from_vega()
                f = units.convert_flux(wave.quantity, f, target_unit, vegaspec=vega)
                if fluxe_col in t.colnames:
                    ef = units.convert_flux(wave.quantity, ef, target_unit, vegaspec=vega)               
        else:
            print('No flux zeropoint exists in the data_dict. We will do an approximation with convert_flux')
            #conversion
            if cid is not None and id is not None:
                raise ValueError('Should have cid, or id, but not both')
            elif cid:
                wave = column_from_column_index(t, cid, data_dict, '**.{{}}.**.{0}.**.{1}'.format(sys, 'lamb'))
            elif id:
                wave = column_from_index(t, id, data_dict, '**.{{}}.**.{0}.**.{1}'.format(sys, 'lamb'))
            vega = SourceSpectrum.from_vega()  # For unit conversion  
            if isinstance(orig_unit, u.MagUnit):
                f = units.convert_flux(wave.quantity, t[flux_col].quantity, target_unit, vegaspec=vega)
                if fluxe_col in t.colnames:
                    msys_plus = t[flux_col].quantity + t[fluxe_col].quantity
                    f_plus = units.convert_flux(wave.quantity, msys_plus, target_unit, vegaspec=vega)
                    ef = abs(f_plus-f)
            else:
                if glom.glom(data_dict, '**.zp'):
                    if cid is not None and id is not None:
                        raise ValueError('Should have cid, or id, but not both')
                    elif cid:
                        zp = column_from_column_index(t, cid, data_dict, '**.{{}}.**.{0}.**.{1}'.format(sys, 'zp'))
                    elif id:
                        zp = column_from_index(t, id, data_dict, '**.{{}}.**.{0}.**.{1}'.format(sys, 'zp')) 
                else:
                    print('No mag zeropoint exists in the data_dict')
                    return []
                minst = fluxToMag(t[flux_col].quantity)
                msys = minst + zp.quantity
                f = units.convert_flux(wave.quantity, msys, target_unit, vegaspec=vega)
                if fluxe_col in t.colnames:
                    _, eminst = fluxErrToMag(t[flux_col].quantity, t[fluxe_col].quantity)
                    if glom.glom(data_dict, '**.e_zP'):
                        e_zp = column_from_column_index(t, 'band', data_dict, '**.{{}}.**.{0}.**.{1}'.format('VEGAMAG', 'e_zP'))
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
    flux_error: Quantity|None 
    cid_col : np.ndarray|None 
    time_format: str
    data_unit: u
    cid:  str|None 
    
    def to_json(self):
        '''This method composes the output pandas DataFrame to a json representation'''
        # if len(self.flux) == len(self.time):
        #     t = QTable([self.time, self.flux, self.flux_error], names=['time', 'flux', 'flux_error'])
        # else:
        #     t = QTable([self.flux, self.flux_error], names=['flux', 'flux_error'])
        if self.flux_error:
            t = QTable([self.flux, self.flux_error], names=['flux', 'flux_error'])
        else:
            t = QTable([self.flux], names=['flux'])
        t['time'] = self.time.to_value(self.time.format) # workaround as QTable cannot convert mixin columns (time) to pandas.
        df = t.to_pandas()
        #include cid column in pandas df
        if self.cid:
            df[self.cid] = self.cid_col.astype('U13')
            #thanks to https://stackoverflow.com/questions/22219004/how-to-group-dataframe-rows-into-list-in-pandas-groupby
            df_g = df.groupby('band')[['time', 'flux', 'flux_error']].agg(list).rename(columns={'time': 'x', 'flux': 'y', 'flux_error': 'error_y'}).to_json(orient='index')
            return df_g
        else:
            return df.rename(columns={'time': 'x', 'flux': 'y', 'flux_error': 'error_y'}).to_json()

#lists of Times and the sanitized astropy.table.Table
@dataclass
class DataProcess:
    mission: str
    time_collection: list[Time]
    table_collection: list[Table]
    system: str|None = field(default = None) #init=False)
    cid: str|None = field(init=False)
    time_format: str = field(init=False)
    data_unit: u = field(default= u.mJy, init=False)
    y_colname: str = field(init=False)
    err_y_colname: str|None = field(init=False)
    timeseries: list[TimeSeries] = field(default_factory=list, init=False)
    
    def __post_init__(self):
        graphic_dict =  DATA_DICT[self.mission]['graphic']
        
        if not self.system:
            try:
                self.system = DATA_DICT[self.mission]['system']
            except KeyError:
                raise Exception("Photometric system undefined in configuration")
        
        try:
            [self.cid] = glom.glom(graphic_dict, '**.cid')
        except:
            self.cid = None
        [self.y_colname] = glom.glom(graphic_dict, '**.y.colname')
        try:
            [self.err_y_colname] = glom.glom(graphic_dict, '**.y_err.colname')
        except:
            self.err_y_colname = None
            
        if self.err_y_colname:
            if self.cid:
                self.timeseries = [TimeSeries(time, table[self.y_colname].quantity, table[self.err_y_colname].quantity, table[self.cid].value, time.format, table[self.y_colname].unit, self.cid) for time, table in zip(self.time_collection, self.table_collection)]
            else:
                self.timeseries = [TimeSeries(time, table[self.y_colname].quantity, table[self.err_y_colname].quantity, time.format, table[self.y_colname].unit, self.cid) for time, table in zip(self.time_collection, self.table_collection)]
        else:
            if self.cid:
                self.timeseries = [TimeSeries(time, table[self.y_colname].quantity, table[self.cid].value, time.format, table[self.y_colname].unit, self.cid) for time, table in zip(self.time_collection, self.table_collection)]
            else:
                self.timeseries = [TimeSeries(time, table[self.y_colname].quantity, time.format, table[self.y_colname].unit) for time, table in zip(self.time_collection, self.table_collection)]
        self.time_format = self.time_collection[-1].format
        self.data_unit = self.table_collection[-1][self.y_colname].unit
     
    def convert_time(self, target_unit): #perhaps a dataclass_transform
        for i in range(len(self.timeseries)):
            if self.timeseries[i].time_format != target_unit:
                self.timeseries[i].time.format = self.timeseries[i].time_format = target_unit
        self.time_format = target_unit

    def convert_flux(self, target_unit):
        zpt_dict = DATA_DICT[self.mission]['zeropt']
        for i in range(len(self.timeseries)):
            if self.timeseries[i].data_unit != target_unit:
                if self.err_y_colname:
                    self.timeseries[i].flux, self.timeseries[i].flux_error = data_convert_unit(self.table_collection[i], self.y_colname, zpt_dict, self.system, cid=self.cid, target_unit=target_unit, fluxe_col=self.err_y_colname)
                else:
                    self.timeseries[i].flux = data_convert_unit(self.table_collection[i], self.y_colname, zpt_dict, self.system, cid=self.cid, target_unit=target_unit)
                self.timeseries[i].data_unit = target_unit
        self.data_unit = target_unit
    
    def to_json(self) -> str:
        return json.dumps(dict(zip(list(str(x) for x in range(len(self.timeseries))), [json.loads(timeseries.to_json()) for timeseries in self.timeseries]))) 
        '''Method to re-structure data by index ['band', 'instrument']'''

    
if __name__ == '__main__':

    tbl = Table.read('/Users/epuga/ESDC/TSViz/data/gaia/anonymous1690191210843_fits/EPOCH_PHOTOMETRY-Gaia DR3 4057091150787104896.fits', astropy_native=True)
    if tbl['flux'].unit == "'electron'.s**-1":
        tbl['flux'].unit = "electron/s"
    if tbl['flux_error'].unit == "'electron'.s**-1":
        tbl['flux_error'].unit = "electron/s"
    
    d = DataProcess('gaia', [tbl['time']], [tbl['source_id', 'band', 'mag', 'flux', 'flux_error']], 'VEGAMAG') 
    d.convert_time('jd')
    print(d.to_json()) 
    d.convert_time('mjd')
    print(d.to_json()) 
    d.convert_flux(u.mJy)
    print(d.to_json())  
    d = DataProcess('gaia', [tbl['time']], [tbl['source_id', 'band', 'mag', 'flux', 'flux_error']], 'VEGAMAG')  
    d.convert_flux(u.mJy)
    print(d.timeseries[0].flux)
    d2 = DataProcess('gaia', [tbl['time']], [tbl['source_id', 'band', 'mag', 'flux', 'flux_error']], 'AB')
    d2.convert_flux(u.mJy)
    print(d2.timeseries[0].flux)
    
    
    