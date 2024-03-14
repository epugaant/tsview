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
import plotly.graph_objects as go

DATA_DICT = {
'gaia': {
    'system': 'VEGAMAG',
    'zeropt': {'G': {'VEGAMAG': {'zp': 25.6873668671 * units.VEGAMAG + 0 * u.mag('s/electron'), 'e_zP': 0.0027553202* units.VEGAMAG + 0 * u.mag('s/electron'), 'lamb': 6217.51 * u.AA, 'f_zp_nu': 3228.75 * u.Jy}, 'AB': {'zp': 25.8010446445 * u.Unit("mag(AB s/electron)"), 'e_zP': 0.0027590522 * u.Unit("mag(AB s/electron)"), 'lamb': 6217.51 * u.AA, 'f_zp_nu': 3631 * u.Jy}}, 
        'BP': {'VEGAMAG': {'zp': 25.3385422158 * units.VEGAMAG + 0 * u.mag('s/electron'), 'e_zP': 0.0027901700* units.VEGAMAG + 0 * u.mag('s/electron'), 'lamb': 5109.71 * u.AA, 'f_zp_nu': 3552.01 * u.Jy}, 'AB': {'zp': 25.3539555559 * u.Unit("mag(AB s/electron)"), 'e_zP': 0.0023065687 * u.Unit("mag(AB s/electron)"), 'lamb': 5109.71 * u.AA, 'f_zp_nu': 3631 * u.Jy}}, 
        'RP': {'VEGAMAG': {'zp': 24.7478955012 * units.VEGAMAG + 0 * u.mag('s/electron'), 'e_zP': 0.0037793818* units.VEGAMAG + 0 * u.mag('s/electron'), 'lamb': 7769.02 * u.AA, 'f_zp_nu': 2554.95 * u.Jy}, 'AB': {'zp': 25.1039837393 * u.Unit("mag(AB s/electron)"), 'e_zP': 0.0015800349 * u.Unit("mag(AB s/electron)"), 'lamb': 7769.02 * u.AA, 'f_zp_nu': 3631 * u.Jy}}, 
        },
    'graphic': {'y': {'colname': 'flux', 'units': 'electron/s'},
              'y_err': {'colname': 'flux_error', 'units': 'electron/s'},
              'cid': 'band',
              'multi': None}    
        },
'jwst': {
    'system': 'VEGAMAG',
    'zeropt': {'P750L': {'VEGAMAG': { 'lamb': 8000 * u.AA}}
               },
    'graphic': {'y': {'colname': 'FLUX', 'units': 'Jy'},
              'y_err': {'colname': 'FLUX_ERROR', 'units': 'Jy'},
              'cid': None,
              'multi': 'FILTER',
              'cextra': 'WAVELENGTH'}    
        }
}

def column_from_column_index(tbl, cid, data_dict, expr_dict):
    '''Function to return a new column and update by rows of a cid (e.g. band), and have the entire table sorted according to the cid values'''
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
    tbl.sort('tmp_col')
    new_col = tbl['tmp_col']
    tbl.remove_column('tmp_col')
    return new_col

def column_from_index(tbl, index, data_dict, expr_dict):
    '''Function to return a new column and update entire table sorted according to the values of the index'''
    target =  glom.glom(data_dict, expr_dict.format(index))
    [val] = glom.flatten(target, levels=(expr_dict.count('**') - 1))
    tbl['tmp_col'] = val
    tbl.sort('tmp_col')
    new_col = tbl['tmp_col']
    tbl.remove_column('tmp_col')
    return new_col

def convert_flux_by_wave(wave, flux, target_unit, vegaspec):
    if not isinstance(wave, u.Quantity):
            raise TypeError('Wave is Not a Quantity.')
    if not isinstance(flux, u.Quantity):
            raise TypeError('Flux is Not a Quantity.')
    try:
        flux_final = units.convert_flux(wave, flux, target_unit, vegaspec=vegaspec)
    #If VEGAMAG, vegaspec will need non-duplicated wavelengths; workaround to group by wavelengths
    except:
        wave_uniq = np.unique(wave.value)
        value = np.empty(0)
        for w in wave_uniq:
            flux_new = units.convert_flux(w, flux[np.where(np.isclose(wave.value, w))], target_unit, vegaspec=vegaspec)
            value = np.append(value, flux_new.value)
        flux_final = u.Quantity(value, flux_new.unit)
    return flux_final

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

def column_factory(t, colname_basis,  data_dict, sys, expr, cid, id):
    '''Wrapper function to handle two ways to create a new column'''
    if cid is not None and id is not None:
        raise ValueError('Should have cid, or id, but not both')
    elif cid:
        col = column_from_column_index(t, cid, data_dict, expr.format(sys, colname_basis))
    elif id:
        col = column_from_index(t, id, data_dict, expr.format(sys, colname_basis))
    return col



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
        
    if orig_unit == target_unit:
        print('Orig unit {0} is equal to target unit {1}'.format(orig_unit, target_unit))
        return 
    
    if orig_unit.physical_type in ('unknown', 'dimensionless'):
        #check if f_zp_mi field exists in dictionary
        if glom.glom(data_dict, '**.f_zp_nu'):
            # if cid is not None and id is not None:
            #     raise ValueError('Should have cid, or id, but not both')
            # elif cid:
            #     f_zp_nu = column_from_column_index(t, cid, data_dict, '**.{{}}.**.{0}.**.{1}'.format(sys, 'f_zp_nu'))
            # elif id:
            #     f_zp_nu = column_from_index(t, id, data_dict, '**.{{}}.**.{0}.**.{1}'.format(sys, 'f_zp_nu'))
            f_zp_nu = column_factory(t, 'f_zp_nu', data_dict, sys, '**.{{}}.**.{0}.**.{1}', cid, id)
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
            # from intermediate units, to target untis
            if target_unit != f_zp_nu.unit:
                if cid is not None and id is not None:
                    raise ValueError('Should have cid, or id, but not both')
                elif cid:
                    wave = column_from_column_index(t, cid, data_dict, '**.{{}}.**.{0}.**.{1}'.format(sys, 'lamb'))
                elif id:
                    wave = column_from_index(t, id, data_dict, '**.{{}}.**.{0}.**.{1}'.format(sys, 'lamb')) 
                vega = SourceSpectrum.from_vega()
                f_int = f
                f = convert_flux_by_wave(wave.quantity, f, target_unit, vegaspec=vega)
                if fluxe_col in t.colnames: 
                    msys_plus = f_int + ef
                    f_plus = convert_flux_by_wave(wave.quantity, msys_plus, target_unit, vegaspec=vega)
                    #ef = abs(f_plus-f)
                    # For the VEGAMAG case, where the unit is not recognized
                    ef = u.Quantity(abs(f_plus.value - f.value), f_plus.unit)          
        else:
            print('No flux zeropoint exists in the data_dict. We will do an approximation with convert_flux')
            #conversion using synphot conver_flux with model when VEGA photometric system
            if cid is not None and id is not None:
                raise ValueError('Should have cid, or id, but not both')
            elif cid:
                wave = column_from_column_index(t, cid, data_dict, '**.{{}}.**.{0}.**.{1}'.format(sys, 'lamb'))
            elif id:
                wave = column_from_index(t, id, data_dict, '**.{{}}.**.{0}.**.{1}'.format(sys, 'lamb'))
            vega = SourceSpectrum.from_vega()  # For unit conversion  
            if isinstance(orig_unit, u.MagUnit):
                f = convert_flux_by_wave(wave.quantity, t[flux_col].quantity, target_unit, vegaspec=vega)
                if fluxe_col in t.colnames:
                    msys_plus = t[flux_col].quantity + t[fluxe_col].quantity
                    f_plus = convert_flux_by_wave(wave.quantity, msys_plus, target_unit, vegaspec=vega)
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
                f = convert_flux_by_wave(wave.quantity, msys, target_unit, vegaspec=vega)
                if fluxe_col in t.colnames:
                    _, eminst = fluxErrToMag(t[flux_col].quantity, t[fluxe_col].quantity)
                    if glom.glom(data_dict, '**.e_zP'):
                        e_zp = column_from_column_index(t, 'band', data_dict, '**.{{}}.**.{0}.**.{1}'.format('VEGAMAG', 'e_zP'))
                        em = np.sqrt(eminst.value**2+ e_zp.value**2) * minst.unit # No astropy unit can do addition in quadrature
                    else:
                        em = eminst
                    minst_plus = minst + em.value * u.mag()
                    msys_plus = minst_plus + zp.quantity
                    f_plus = convert_flux_by_wave(wave.quantity, msys_plus, target_unit, vegaspec=vega)
                    ef = abs(f_plus-f)
    else:
        if cid is not None and id is not None:
            raise ValueError('Should have cid, or id, but not both')
        elif cid:
            wave = column_from_column_index(t, cid, data_dict, '**.{{}}.**.{0}.**.{1}'.format(sys, 'lamb'))
        elif id:
            wave = column_from_index(t, id, data_dict, '**.{{}}.**.{0}.**.{1}'.format(sys, 'lamb'))
        vega = SourceSpectrum.from_vega() 
        f = convert_flux_by_wave(wave.quantity, t[flux_col].quantity, target_unit, vegaspec=vega)
        if fluxe_col in t.colnames:
            msys_plus = t[flux_col].quantity + t[fluxe_col].quantity
            f_plus = convert_flux_by_wave(wave.quantity, msys_plus, target_unit, vegaspec=vega)
            ef = abs(f_plus-f) 
                    
    if fluxe_col in t.colnames:
        return (f, ef)
    else:
        return (f)


@dataclass
class TimeSeries:
    time: u.Quantity
    time_format: str
    flux: u.Quantity
    data_unit: u
    flux_error: u.Quantity|None = field(default = None)
    id_col : np.ndarray|None = field(default = None)
    id:  str|None = field(default = None)
    
    def to_pandas(self):
        '''This method converts a QTable into a pandas DataFrame for organization'''
        if self.flux_error is not None:
            t = QTable([self.flux, self.flux_error], names=['flux', 'flux_error'])
        else:
            t = QTable([self.flux], names=['flux'])
        t['time'] = self.time.to_value(self.time.format) # workaround as QTable cannot convert mixin columns (time) to pandas Dataframe
        df = t.to_pandas()
        return df
    
    def to_json(self):
        '''This method composes the output pandas DataFrame to a json representation using groupby on cid or directly converts to string'''
        df = self.to_pandas()
        #include cid column in pandas df
        if len(self.id_col) == len(self.flux):
            df[self.id] = self.id_col.astype('U13')
            #thanks to https://stackoverflow.com/questions/22219004/how-to-group-dataframe-rows-into-list-in-pandas-groupby
            df_g = df.groupby(self.id)[['time', 'flux', 'flux_error']].agg(list).rename(columns={'time': 'x', 'flux': 'y', 'flux_error': 'error_y'}).to_json(orient='index')
            return df_g
        elif len(self.id_col) == 1:
            return json.dumps({self.id_col[0]: df.rename(columns={'time': 'x', 'flux': 'y', 'flux_error': 'error_y'}).to_dict(orient='list')})
        else:
            return json.dumps(df.rename(columns={'time': 'x', 'flux': 'y', 'flux_error': 'error_y'}).to_dict(orient='list'))
        


#lists of Times and the sanitized astropy.table.Table
@dataclass
class DataProcess:
    mission: str
    time_collection: list[Time]
    table_collection: list[Table]
    system: str|None = field(default = None) #init=False)
    cid: str|None = field(init=False)
    multi: str|None = field(init=False)
    time_format: str = field(init=False)
    time_scale: str = field(init=False)
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
        try:
            [self.multi] = glom.glom(graphic_dict, '**.multi')
        except:
            self.multi = None

        [self.y_colname] = glom.glom(graphic_dict, '**.y.colname')
        try:
            [self.err_y_colname] = glom.glom(graphic_dict, '**.y_err.colname')
        except:
            self.err_y_colname = None
            
        if self.err_y_colname:
            if self.cid:
                self.timeseries = [TimeSeries(time, time.format, table[self.y_colname].quantity, table[self.y_colname].unit, table[self.err_y_colname].quantity, table[self.cid].value, self.cid) for time, table in zip(self.time_collection, self.table_collection)]
            elif self.multi:
                self.timeseries = [TimeSeries(time, time.format, table[self.y_colname].quantity, table[self.y_colname].unit, table[self.err_y_colname].quantity, np.array([table.meta[self.multi]], dtype='object'), self.multi) for time, table in zip(self.time_collection, self.table_collection)]
            else:
                self.timeseries = [TimeSeries(time, time.format, table[self.y_colname].quantity, table[self.y_colname].unit, table[self.err_y_colname].quantity) for time, table in zip(self.time_collection, self.table_collection)]
        else:
            if self.cid:
                self.timeseries = [TimeSeries(time, time.format, table[self.y_colname].quantity, table[self.y_colname].unit, table[self.cid].value, self.cid) for time, table in zip(self.time_collection, self.table_collection)]
            elif self.multi:
                self.timeseries = [TimeSeries(time, time.format, table[self.y_colname].quantity, table[self.y_colname].unit, np.array([table.meta[self.multi]], dtype='object'), self.multi) for time, table in zip(self.time_collection, self.table_collection)]                
            else:
                self.timeseries = [TimeSeries(time, time.format, table[self.y_colname].quantity, table[self.y_colname].unit) for time, table in zip(self.time_collection, self.table_collection)]
        self.time_format = self.time_collection[-1].format
        self.time_scale = self.time_collection[-1].scale
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
                    if self.cid:
                        self.timeseries[i].flux, self.timeseries[i].flux_error = data_convert_unit(self.table_collection[i], self.y_colname, zpt_dict, self.system, cid=self.cid, target_unit=target_unit, fluxe_col=self.err_y_colname)
                    elif self.multi:
                        self.timeseries[i].flux, self.timeseries[i].flux_error = data_convert_unit(self.table_collection[i], self.y_colname, zpt_dict, self.system, id=self.table_collection[i].meta[self.multi], target_unit=target_unit, fluxe_col=self.err_y_colname)
                else:
                    if self.cid:
                        self.timeseries[i].flux = data_convert_unit(self.table_collection[i], self.y_colname, zpt_dict, self.system, cid=self.cid, target_unit=target_unit)
                    elif self.multi:
                        self.timeseries[i].flux = data_convert_unit(self.table_collection[i], self.y_colname, zpt_dict, self.system, id=self.table_collection[i].meta[self.multi], target_unit=target_unit)
                self.timeseries[i].data_unit = target_unit
        self.data_unit = target_unit
    
    def to_json(self) -> str:
        return json.dumps(dict(zip(list(str(x) for x in range(len(self.timeseries))), [json.loads(timeseries.to_json()) for timeseries in self.timeseries]))) 
        '''Method to re-structure data by index ['band', 'instrument']'''
        
    def to_plotly(self) -> str:
        fig = go.Figure()
        for timeseries in self.timeseries:
            df = timeseries.to_pandas().rename(columns={'time': 'x', 'flux': 'y', 'flux_error': 'error_y'})
            if len(timeseries.id_col) == len(timeseries.flux):
                df[timeseries.id] = timeseries.id_col.astype('U13')
                #thanks to https://stackoverflow.com/questions/22219004/how-to-group-dataframe-rows-into-list-in-pandas-groupby
                for index, df_group in df.groupby(timeseries.id):
                    x = df_group.x
                    y = df_group.y
                    if timeseries.flux_error is None:
                        error_y = None
                    else:
                        error_y = dict(
                                type='data', # value of error bar given in data coordinates
                                array=df_group.error_y,
                                visible=True)
                    fig.add_trace(go.Scatter(x=x, y=y, error_y=error_y, name=index)) 
            elif len(timeseries.id_col) == 1: 
                x = df.x
                y = df.y
                if timeseries.flux_error is None:
                    error_y = None
                else:
                    error_y = dict(
                            type='data', # value of error bar given in data coordinates
                            array=df.error_y,
                            visible=True)
                fig.add_trace(go.Scatter(x=x, y=y, error_y=error_y, name=timeseries.id_col[0]))
            else:
                x = df.x
                y = df.y
                if timeseries.flux_error is None:
                    error_y = None
                else:
                    error_y = dict(
                            type='data', # value of error bar given in data coordinates
                            array=df.error_y,
                            visible=True)
                fig.add_trace(go.Scatter(x=x, y=y, error_y=error_y))
        fig.update_layout(legend_title_text = self.mission)
        fig.update_xaxes(title_text='Time [{0} in {1}]'.format(self.time_format, self.time_scale.upper()))
        fig.update_yaxes(title_text='{0} [{1}]'.format(self.y_colname.capitalize(), self.data_unit.to_string()))
        return fig.to_json(validate=True)

    
if __name__ == '__main__':

    tbl = Table.read('/Users/epuga/ESDC/TSViz/data/gaia/anonymous1690191210843_fits/EPOCH_PHOTOMETRY-Gaia DR3 4057091150787104896.fits', astropy_native=True)
    if tbl['flux'].unit == "'electron'.s**-1":
        tbl['flux'].unit = "electron/s"
    if tbl['flux_error'].unit == "'electron'.s**-1":
        tbl['flux_error'].unit = "electron/s"
    
    d = DataProcess('gaia', [tbl['time']], [tbl['source_id', 'band', 'mag', 'flux', 'flux_error']], 'VEGAMAG') 
    #d.convert_time('jd')
    #print(d.to_json()) 
    #d.convert_flux(u.mJy)
    #print(d.to_json()) 
    d.convert_flux(units.VEGAMAG)
    print(d.to_json()) 
    print(d.to_plotly())
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
    
    tbl2 =  Table.read('/Users/epuga/ESDC/TSViz/data/gaia/swo/Gaia DR3 4057091150787104896_ALL_INDIVIDUAL_fits/EPOCH_PHOTOMETRY-Gaia DR3 4057091150787104896.fits', astropy_native=True)
    
    import os
    from tsview import DATADIR
    from tsview.io.fits_parser import ts_fits_reader
    
    filename = 'jw02783-o002_t001_miri_p750l-slitlessprism_x1dints.fits'
    time, data = ts_fits_reader(os.path.join(DATADIR, filename))
    new_data = [tbl['WAVELENGTH', 'FLUX', 'FLUX_ERROR'] for tbl in data]
    d = DataProcess('jwst', time, new_data, 'VEGAMAG')
    d.convert_time('jd')
    print(d.to_json())
    d.convert_flux(u.mJy)
    print(d.to_plotly())
    pass
    
    