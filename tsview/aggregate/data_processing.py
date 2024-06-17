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

from tsview.aggregate.data_config import DATA_DICT 


COLS = ['time', 'flux', 'flux_error']
NEW_COLS = ['x', 'y', 'error_y']
    

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
    #tbl.sort('tmp_col')
    new_col = tbl['tmp_col']
    tbl.remove_column('tmp_col')
    return new_col

def column_from_index(tbl, index, data_dict, expr_dict):
    '''Function to return a new column and update entire table sorted according to the values of the index'''
    target =  glom.glom(data_dict, expr_dict.format(index))
    [val] = glom.flatten(target, levels=(expr_dict.count('**') - 1))
    tbl['tmp_col'] = val
    #tbl.sort('tmp_col')
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
    
    #if data product units are instrumental (uncalibrated)
    if orig_unit.physical_type in ('unknown', 'dimensionless'):
        #check if f_zp_nu field exists in dictionary
        if glom.glom(data_dict, '**.f_zp_nu'):
            f_zp_nu = column_factory(t, 'f_zp_nu', data_dict, sys, '**.{{}}.**.{0}.**.{1}', cid, id)
            if isinstance(orig_unit, u.MagUnit):
                if fluxe_col in t.colnames:
                    f, ef = f_zp_nu.quantity * magErrToFlux(t[flux_col].quantity, t[fluxe_col].quantity)# [u.Jy] * [dimensionless]
                else:
                    f = f_zp_nu.quantity * magToFlux(t[flux_col].quantity) 
            #convert to logaritmic scales for compatibility with zp
            else:
                if glom.glom(data_dict, '**.zp'):
                    zp = column_factory(t, 'zp', data_dict, sys, '**.{{}}.**.{0}.**.{1}', cid, id)
                else:
                    print('No mag zeropoint exists in the data_dict')
                    return []
                if fluxe_col in t.colnames:
                    minst, eminst = fluxErrToMag(t[flux_col].quantity, t[fluxe_col].quantity)
                    #include zp uncertainty in propagation, if you have it
                    if glom.glom(data_dict, '**.e_zP'):
                        e_zp = column_from_column_index(t, 'band', data_dict, '**.{{}}.**.{0}.**.{1}'.format(sys, 'e_zP'))
                        em = np.sqrt(eminst.value**2+ e_zp.value**2) * minst.unit # No astropy unit can do addition in quadrature
                    else:
                        em = eminst
                    em = em + (0 * zp.unit)
                else:
                    minst = fluxToMag(t[flux_col].quantity)
                msys = minst + zp.quantity
                if target_unit != units.VEGAMAG:
                    # if target_unit != (zp.unit - u.mag(1/orig_unit)):
                    #     msys = minst + zp.quantity
                    # else:
                    #     msys = minst
                    
                    #now back to flux
                    if fluxe_col in t.colnames:
                        f, ef = f_zp_nu.quantity * magErrToFlux(msys, em)
                    else:
                        f = f_zp_nu.quantity * magToFlux(msys)
                else:
                    if fluxe_col in t.colnames:
                        f, ef = msys, em
                    else:
                        f = msys
                # from intermediate units, to target untis
            if target_unit != f_zp_nu.unit and target_unit != units.VEGAMAG:
                wave = column_factory(t, 'lamb', data_dict, sys, '**.{{}}.**.{0}.**.{1}', cid, id)
                vega = SourceSpectrum.from_vega()
                f_int = f
                f = convert_flux_by_wave(wave.quantity, f, target_unit, vegaspec=vega)
                if fluxe_col in t.colnames: 
                    f_plus = f_int + ef
                    f_plus = convert_flux_by_wave(wave.quantity, f_plus, target_unit, vegaspec=vega)
                    #ef = abs(f_plus-f)
                    # For the VEGAMAG case, where the unit is not recognized
                    ef = u.Quantity(abs(f_plus.value - f.value), f_plus.unit)          
        else:
            print('No flux zeropoint exists in the data_dict. We will do an approximation with convert_flux')
            #conversion using synphot conver_flux with model when VEGA photometric system
            wave = column_factory(t, 'lamb', data_dict, sys, '**.{{}}.**.{0}.**.{1}', cid, id)
            vega = SourceSpectrum.from_vega()  # For unit conversion  
            if isinstance(orig_unit, u.MagUnit):
                f = convert_flux_by_wave(wave.quantity, t[flux_col].quantity, target_unit, vegaspec=vega)
                if fluxe_col in t.colnames:
                    msys_plus = t[flux_col].quantity + t[fluxe_col].quantity
                    f_plus = convert_flux_by_wave(wave.quantity, msys_plus, target_unit, vegaspec=vega)
                    ef = abs(f_plus-f)
            else:
                if glom.glom(data_dict, '**.zp'):
                     zp = column_factory(t, 'zp', data_dict, sys, '**.{{}}.**.{0}.**.{1}', cid, id)
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
    #if physical units are calibrated (e.g. mag(AB))
    else:
        wave = column_factory(t, 'lamb', data_dict, sys, '**.{{}}.**.{0}.**.{1}', cid, id)
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

PHOT_SYS = ['AB', 'ST', 'VEGA']
PHOT_SYS_MAG_UNIT = {'mag(AB)': u.ABmag, 'mag(ST)': u.STmag, 'mag(VEGA)': units.VEGAMAG}

SPECT_DENS = [
            'Jy', 'mJy', 'uJy',
            'W / (m2 Hz)', #'W / (Hz m2)',  # Order is different in astropy v5.3
            'eV / (s m2 Hz)', #'eV / (Hz s m2)',
            'erg / (s cm2)',
            'erg / (s cm2 Angstrom)', #'erg / (Angstrom s cm2)',
            'erg / (s cm2 Hz)', #'erg / (Hz s cm2)',
            'ph / (s cm2 Angstrom)', #'ph / (Angstrom s cm2)',
            'ph / (s cm2 Hz)', #'ph / (Hz s cm2)'
                ]

def equivalent_units(unit_str):
    '''Function to return spectral flux and spectral equivalences in BaseUnit and MagUnit as set of strings'''
    try:
        units = PHOT_SYS_MAG_UNIT[unit_str]
    except:
        units = u.Unit(unit_str)
        
    if units.physical_type in ['spectral flux density'] or (hasattr(units, 'physical_unit') and units.physical_unit != 'dimensionless'): # spectral flux
        eqv = u.spectral_density(1 * u.m)  # Value does not matter here.
        exclude_lower = {'flam', 'fnu', 'bol', 'photlam', 'photnu'}
        exclude_upper = {'FLAM', 'FNU', 'BOL', 'PHOTLAM', 'PHOTNU'}
        try:
            list_of_units = set(list(map(str, units.find_equivalent_units(
                include_prefix_units=False, equivalencies=eqv))) + PHOT_SYS + 
                                SPECT_DENS) - exclude_lower - exclude_upper
        except:
            list_of_units = set(PHOT_SYS + 
                                SPECT_DENS) - exclude_lower - exclude_upper
        
        # include magnitudes for three photometric systems
        list_of_units = list_of_units | {'mag({})'.format(val) for val in list_of_units if val in PHOT_SYS}
        # remove VEGA flux
        list_of_units.remove('VEGA')
        # remove self
        list_of_units.remove(unit_str)
    elif units.physical_type in ['length', 'frequency', 'wavenumber', 'energy']:  # spectral axis
        # prefer Hz over Bq and um over micron
        exclude = {'Bq', 'micron'}
        list_of_units = set(list(map(str, units.find_equivalent_units(
            include_prefix_units=False, equivalencies=u.spectral()))) + ['um']) - exclude
    else:
        list_of_units = set([])
    return list(list_of_units)

def time_units(time):
    '''Function to return time equivalences in as set of strings'''
    exclude = {'cxcsec', 'datetime', 'gps', 'unix', 'unix_tai', 'ymdhms', 'datetime64'}
    list_of_units = set(list(map(str, time.FORMATS.keys()))) - exclude
    # remove self
    list_of_units.remove(time.format)
    return list(list_of_units)


@dataclass
class TimeSeries:
    time: u.Quantity
    time_format: str
    flux: u.Quantity
    data_unit: u
    flux_error: u.Quantity|None = field(default = None)
    id_col : np.ndarray|None = field(default = None)
    id:  str|None = field(default = None)
    extra_col : u.Quantity|None = field(default = None)
    extra:  str|None = field(default = None)
    
    def to_pandas(self):
        '''This method converts a QTable into a pandas DataFrame for organization'''
        if self.flux_error is not None:
            t = QTable([self.flux, self.flux_error], names=['flux', 'flux_error'])
        else:
            t = QTable([self.flux], names=['flux'])
        t['time'] = self.time.to_value(self.time.format) # workaround as QTable cannot convert mixin columns (time) to pandas Dataframe
        if self.extra is not None:
            t['extra'] = self.extra_col
        df = t.to_pandas()
        return df
    
    def to_json(self):
        '''This method composes the output pandas DataFrame to a json representation using groupby on cid or directly converts to string'''
        if self.extra is not None:
            COLS.append('extra')
            NEW_COLS.append('z')
            
        df = self.to_pandas()
        #include cid column in pandas df
        if len(self.id_col) == len(self.flux):
            df[self.id] = self.id_col.astype('U13')
            #thanks to https://stackoverflow.com/questions/22219004/how-to-group-dataframe-rows-into-list-in-pandas-groupby
            df_g = df.groupby(self.id)[COLS].agg(list).rename(columns=dict(zip(COLS, NEW_COLS))).to_json(orient='index')
            return df_g
        elif len(self.id_col) == 1:
            return json.dumps({self.id_col[0]: df.rename(columns=dict(zip(COLS, NEW_COLS))).to_dict(orient='list')})
        else:
            return json.dumps(df.rename(columns=dict(zip(COLS, NEW_COLS))).to_dict(orient='list'))
        


#lists of Times and the sanitized astropy.table.Table
@dataclass
class DataProcess:
    mission: str
    time_collection: list[Time]
    table_collection: list[Table]
    system: str|None = field(default = None) #init=False)
    cextra: str|None = field(init=False)
    cid: str|None = field(init=False)
    multi: str|None = field(init=False)
    time_format: str = field(init=False)
    alt_time_format: list = field(init=False)
    time_scale: str = field(init=False)
    data_unit: u = field(default= u.mJy, init=False)
    alt_data_unit: list = field(init=False)
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
        try:
            [self.cextra] = glom.glom(graphic_dict, '**.cextra')
        except:
            self.cextra = None

        [self.y_colname] = glom.glom(graphic_dict, '**.y.colname')
        try:
            [self.err_y_colname] = glom.glom(graphic_dict, '**.y_err.colname')
        except:
            self.err_y_colname = None
            
        if self.err_y_colname:
            if self.cid:
                if self.cextra:
                    self.timeseries = [TimeSeries(time, time.format, table[self.y_colname].quantity, table[self.y_colname].unit, table[self.err_y_colname].quantity, table[self.cid].value, self.cid, table[self.cextra].quantity, self.cextra) for time, table in zip(self.time_collection, self.table_collection)]
                else:
                    self.timeseries = [TimeSeries(time, time.format, table[self.y_colname].quantity, table[self.y_colname].unit, table[self.err_y_colname].quantity, table[self.cid].value, self.cid) for time, table in zip(self.time_collection, self.table_collection)]
            elif self.multi:
                if self.cextra:
                    self.timeseries = [TimeSeries(time, time.format, table[self.y_colname].quantity, table[self.y_colname].unit, table[self.err_y_colname].quantity, np.array([table.meta[self.multi]], dtype='object'), self.multi, table[self.cextra].quantity, self.cextra) for time, table in zip(self.time_collection, self.table_collection)]
                else:
                    self.timeseries = [TimeSeries(time, time.format, table[self.y_colname].quantity, table[self.y_colname].unit, table[self.err_y_colname].quantity, np.array([table.meta[self.multi]], dtype='object'), self.multi) for time, table in zip(self.time_collection, self.table_collection)]
            else:
                if self.cextra:
                    self.timeseries = [TimeSeries(time, time.format, table[self.y_colname].quantity, table[self.y_colname].unit, table[self.err_y_colname].quantity, table[self.cextra].quantity, self.cextra) for time, table in zip(self.time_collection, self.table_collection)]
                else:
                    self.timeseries = [TimeSeries(time, time.format, table[self.y_colname].quantity, table[self.y_colname].unit, table[self.err_y_colname].quantity) for time, table in zip(self.time_collection, self.table_collection)]
        else:
            if self.cid:
                if self.cextra:
                    self.timeseries = [TimeSeries(time, time.format, table[self.y_colname].quantity, table[self.y_colname].unit, table[self.cid].value, self.cid, table[self.cextra].quantity, self.cextra) for time, table in zip(self.time_collection, self.table_collection)]
                else:
                    self.timeseries = [TimeSeries(time, time.format, table[self.y_colname].quantity, table[self.y_colname].unit, table[self.cid].value, self.cid) for time, table in zip(self.time_collection, self.table_collection)]
            elif self.multi:
                if self.cextra:
                    self.timeseries = [TimeSeries(time, time.format, table[self.y_colname].quantity, table[self.y_colname].unit, np.array([table.meta[self.multi]], dtype='object'), self.multi, table[self.cextra].quantity, self.cextra) for time, table in zip(self.time_collection, self.table_collection)]
                else: 
                    self.timeseries = [TimeSeries(time, time.format, table[self.y_colname].quantity, table[self.y_colname].unit, np.array([table.meta[self.multi]], dtype='object'), self.multi) for time, table in zip(self.time_collection, self.table_collection)]                
            else:
                if self.cextra:
                    self.timeseries = [TimeSeries(time, time.format, table[self.y_colname].quantity, table[self.y_colname].unit, table[self.cextra].quantity, self.cextra ) for time, table in zip(self.time_collection, self.table_collection)]
                else:
                    self.timeseries = [TimeSeries(time, time.format, table[self.y_colname].quantity, table[self.y_colname].unit) for time, table in zip(self.time_collection, self.table_collection)]
                
        self.time_format = self.time_collection[-1].format
        self.alt_time_format = time_units(self.timeseries[-1].time)
        self.time_scale = self.time_collection[-1].scale
        self.data_unit = self.table_collection[-1][self.y_colname].unit
        self.alt_data_unit = equivalent_units(self.data_unit.to_string())
     
    def convert_time(self, target_unit): #perhaps a dataclass_transform
        for i in range(len(self.timeseries)):
            if self.timeseries[i].time_format != target_unit:
                self.timeseries[i].time.format = self.timeseries[i].time_format = target_unit
        self.time_format = target_unit
        self.alt_time_format = time_units(self.timeseries[-1].time)

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
        self.alt_data_unit = equivalent_units(self.data_unit.to_string())
    
    def to_json(self) -> str:
        return json.dumps(dict(zip(list(str(x) for x in range(len(self.timeseries))), [json.loads(timeseries.to_json()) for timeseries in self.timeseries]))) 
        '''Method to re-structure data by index ['band', 'instrument']'''

    def create_scatter(self, x, y, error_y, index=None, z=None) -> go.Scatter:
        scatter = go.Scatter(x=x, y=y, error_y=error_y, name=index)
        scatter.mode = 'markers'
        scatter.hovertemplate = r'%{yaxis.title.text}: %{y}<br>%{xaxis.title.text}: %{x}'
        if z:
            scatter.marker=dict(
                cmax=max(z),
                cmin=min(z),
                color=z,
                colorbar=dict(
                    title='cextra'
                ),
                colorscale="Turbo"
            )
        return scatter
    

    def to_plotly(self, time=True) -> str:
        '''Function to generate the plotly.Figure object using graph_object'''
        if self.cextra is not None:
            COLS.append('extra')
            NEW_COLS.append('z')
        fig = go.Figure()
        for timeseries in self.timeseries:
            df = timeseries.to_pandas().rename(columns=dict(zip(COLS, NEW_COLS)))
            if len(timeseries.id_col) == len(timeseries.flux):
                df[timeseries.id] = timeseries.id_col.astype('U13')
                #thanks to https://stackoverflow.com/questions/22219004/how-to-group-dataframe-rows-into-list-in-pandas-groupby
                for index, df_group in df.groupby(timeseries.id):
                    x = df_group.x 
                    y = df_group.y
                    if self.cextra:
                        z = df_group.z
                    else:
                        tbl = Table([timeseries.id_col.data], names=[timeseries.id])
                        z = column_factory(tbl, 'lamb', DATA_DICT, self.system, '**.{{}}.**.{0}.**.{1}', self.cid, self.multi).quantity.to(u.AA, equivalencies=u.spectral()).value
                    if timeseries.flux_error is None:
                        error_y = None
                    else:
                        error_y = dict(
                                type='data', # value of error bar given in data coordinates
                                array=df_group.error_y,
                                visible=True)
                    fig.add_trace(self.create_scatter(x, y, error_y, index)) if time else fig.add_trace(self.create_scatter(z, y, error_y, index))
            elif len(timeseries.id_col) == 1: 
                x = df.x 
                y = df.y
                if self.cextra:
                    z = df.z
                else:
                    tbl = Table([timeseries.id_col], names=[timeseries.id])
                    z = column_factory(tbl, 'lamb', DATA_DICT, self.system, '**.{{}}.**.{0}.**.{1}', self.cid, self.multi).quantity.to(u.AA, equivalencies=u.spectral()).value
                if timeseries.flux_error is None:
                    error_y = None
                else:
                    error_y = dict(
                            type='data', # value of error bar given in data coordinates
                            array=df.error_y,
                            visible=True)
                fig.add_trace(self.create_scatter(x, y, error_y, timeseries.id_col[0]), z=self.cextra) if time else fig.add_trace(self.create_scatter(z, y, error_y, timeseries.id_col[0]))
            else:
                x = df.x 
                y = df.y
                if self.cextra:
                    z = df.z
                else:
                    tbl = Table([timeseries.id_col], names=[timeseries.id])
                    z = column_factory(tbl, 'lamb', DATA_DICT, self.system, '**.{{}}.**.{0}.**.{1}', self.cid, self.multi).quantity.to(u.AA, equivalencies=u.spectral()).value
                if timeseries.flux_error is None:
                    error_y = None
                else:
                    error_y = dict(
                            type='data', # value of error bar given in data coordinates
                            array=df.error_y,
                            visible=True)
                fig.add_trace(self.create_scatter(x, y, error_y)) if time else fig.add_trace(self.create_scatter(z, y, error_y))
        fig.update_layout(legend_title_text = self.mission)
        fig.update_xaxes(title_text='Time [{0} in {1}]'.format(self.time_format, self.time_scale.upper())) if time else fig.update_xaxes(title_text='Wave [{0}]'.format(u.AA.to_string().upper()))
        fig.update_yaxes(title_text='{0} [{1}]'.format(self.y_colname.capitalize(), self.data_unit.to_string()))

        plotly_json = fig.to_json(validate=True)
        
        json_dictionary = json.loads(plotly_json)
        json_dictionary['alt_time_units'] = self.alt_time_format
        json_dictionary['alt_flux_units'] = self.alt_data_unit
        return json.dumps(json_dictionary)

    
if __name__ == '__main__':


    tbl = Table.read('/Users/epuga/ESDC/TSViz/data/gaia/anonymous1690191210843_fits/EPOCH_PHOTOMETRY-Gaia DR3 4057091150787104896.fits', astropy_native=True)
    if tbl['flux'].unit == "'electron'.s**-1":
        tbl['flux'].unit = "electron/s"
    if tbl['flux_error'].unit == "'electron'.s**-1":
        tbl['flux_error'].unit = "electron/s"
    
    d = DataProcess('gaia', [tbl['time']], [tbl['source_id', 'band', 'mag', 'flux', 'flux_error']], 'VEGAMAG') 
    #d.convert_time('jd')
    #print(d.to_json()) 
    d.convert_flux(u.mJy)
    #print(d.to_json()) 
    print(d.timeseries[0].id_col)
    print(d.timeseries[0].flux)
    d.convert_flux(units.VEGAMAG)
    print(d.timeseries[0].id_col)
    print(d.timeseries[0].flux)
    #print(d.to_json()) 
    print(d.to_plotly(time=True))
    print(d.to_plotly(time=False)) #plotly wave vs flux

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
    
    tbl2 =  Table.read('/Users/epuga/ESDC/TSViz/data/gaia/swo/Gaia DR3 4057091150787104896_ALL_INDIVIDUAL_fits/EPOCH_PHOTOMETRY-Gaia DR3 4057091150787104896.fits', astropy_native=True)
    
    import os
    from tsview import DATADIR
    from tsview.io.fits_parser import ts_fits_reader
    
    filename = 'jw02783-o002_t001_miri_p750l-slitlessprism_x1dints.fits'
    time, data = ts_fits_reader(os.path.join(DATADIR, filename))
    new_data = [tbl['WAVELENGTH', 'FLUX', 'FLUX_ERROR'] for tbl in data]
    d = DataProcess('jwst', time, new_data, 'VEGAMAG')
    print(d.timeseries[0].extra_col)
    #d.convert_time('jd')
    print(d.to_json())
    d.convert_flux(u.mJy)
    print(d.to_plotly(time=False))
    print(d.timeseries[0].extra)
    pass
    
    