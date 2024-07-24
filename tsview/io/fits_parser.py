# Licensed under a 3-clause BSD style license - see LICENSE.rst
import os
import warnings
import numpy as np

from astropy.io import registry, fits
from astropy.table import Table
from astropy.time import Time, TimeDelta

from tsview import DATADIR
from tsview.utils.timer import timer_func

__all__ = ["ts_fits_reader"]

def sanitize_weakrefs(times):
    if any(time.scale in ['tt'] for time in times):
        times_new = [Time(time.value, scale='tt', format=time.format) for time in times if time.scale == 'tt']
        return times_new
    else:
        return times

def locate_extension_with_keyword(tinfo, filename, keyword):
    for i in range(len(tinfo)):
        hdr = fits.getheader(filename, i, ignore_missing_simple=True)
        if keyword in hdr:
            break
    return i, hdr

def int_times_ext_to_time_arr(filename, scale=None):
    '''Function to read FITS extension INT_TIMES as table and generate a single time object from int_mid column in mjd'''
    # Time
    t = Table.read(filename, hdu="INT_TIMES", astropy_native=True) 
    if scale is None:
        fmt = 'BJD'.lower() 
    if any(fmt in col.lower() for col in t.colnames):   
        # subset of BJD columns
        sub_cols = [col for col in t.colnames if fmt in col.lower()]
        # sub-subset of int_mid colum
        int_mid = [col for col in sub_cols if 'mid' in col]
        # if there is only one column
        if len(int_mid) == 1:
            # get the scale from the name of the column
            scale = int_mid[0].split('_')[-1].lower() # 'tdb'
            format = 'mjd'
    
    time = Time(t[int_mid[0]].data, format=format, scale=scale)

    return time

def extract1d_to_timeseries(filename, times, data, hdr, tinfo, ts_keys, scale=None, array_int_times = False):
    if scale is None:
        scale = 'TDB'.lower()
    #times, data = [], []
    hdr_time = dict(filter(lambda i:i[0] in ts_keys, hdr.items())) 
    mask = tinfo['Name'] == 'EXTRACT1D'
    extver = tinfo['Ver'][mask].data.tolist()
    if array_int_times:
        times.append(int_times_ext_to_time_arr(filename))
    for ver in extver:
        tab = Table.read(filename, hdu=("EXTRACT1D", ver), astropy_native=True)
        # propagate time information from the primary to all extract1d
        tab.meta.update(hdr_time)
        data.append(tab)
        if not array_int_times:
            times.append(Time(tab.meta[scale.upper()+'-MID'], scale=scale, format='mjd'))
            
    return times, data

def TSrebin(tstable, DTin, DTout):
    """
    Re-binning function
    
    Input: TimeSeries Table, Input Time bin and Output Time Bin
    
    #*****************************************************************************
    #
    #    WARNING: DTin should be previously read from the FITS Keyword: 'TIMEDEL'
    #
    #*****************************************************************************    
    """
    import sys
    import warnings
 
    import numpy as np
    from astropy.table import Table
    #
    # Check input Table Columns
    #
    colNames         = tstable.keys()                           # Input ColumnNames
    mandatoryColumns = ['TIME','RATE','ERROR','FRACEXP']        # Mandatory ColumnNames
    diff = set(mandatoryColumns) - set(colNames)
    fracexp_chk = 0
    if diff == {'FRACEXP'}:
      #print('FRACEXP is missing. fracexp = 1')
      fracexp_chk = 1
    elif len(diff) != 0:
      print('    Missing Columns in input Table : ', diff, '\n Exiting.')
      sys.exit(1)
 
    #
    # Re-Binning of th Table 'tstable'
    #
    tstable.sort('TIME')
    ntime=tstable['TIME']
    nrate=np.nan_to_num(tstable['RATE'])
    nerror=np.nan_to_num(tstable['ERROR'])

    nrate1=np.nan_to_num(tstable['RATE1'])
    nerror1=np.nan_to_num(tstable['RATE1_ERR'])

    nrate2=np.nan_to_num(tstable['RATE2'])
    nerror2=np.nan_to_num(tstable['RATE2_ERR'])

    nrate3=np.nan_to_num(tstable['RATE3'])
    nerror3=np.nan_to_num(tstable['RATE3_ERR'])

    nrate4=np.nan_to_num(tstable['RATE4'])
    nerror4=np.nan_to_num(tstable['RATE4_ERR'])

    nrate5=np.nan_to_num(tstable['RATE5'])
    nerror5=np.nan_to_num(tstable['RATE5_ERR'])

    
    if fracexp_chk == 1:
      ### fracexp = 1
      fracexp=np.linspace(1, 1, len(ntime))
    else:
      fracexp=np.nan_to_num(tstable['FRACEXP'])
 
    tcum=np.nancumsum(tstable['TIME'].value)
    rcum=np.nancumsum(nrate*fracexp)                            # Errors propagation are weighted by the FRACEXP from the TS
    ecum=np.nancumsum(nerror**2*fracexp)

    rcum1=np.nancumsum(nrate1*fracexp)                            # Errors propagation are weighted by the FRACEXP from the TS
    ecum1=np.nancumsum(nerror1**2*fracexp)

    rcum2=np.nancumsum(nrate2*fracexp)                            # Errors propagation are weighted by the FRACEXP from the TS
    ecum2=np.nancumsum(nerror2**2*fracexp)

    rcum3=np.nancumsum(nrate3*fracexp)                            # Errors propagation are weighted by the FRACEXP from the TS
    ecum3=np.nancumsum(nerror3**2*fracexp)

    rcum4=np.nancumsum(nrate4*fracexp)                            # Errors propagation are weighted by the FRACEXP from the TS
    ecum4=np.nancumsum(nerror4**2*fracexp)

    rcum5=np.nancumsum(nrate5*fracexp)                            # Errors propagation are weighted by the FRACEXP from the TS
    ecum5=np.nancumsum(nerror5**2*fracexp)

    fcum=np.nancumsum(fracexp)
    nin=len(rcum)                                               # Number of In points
 
    t_bin_frac   = int(np.round( DTout / DTin ))
 
    if t_bin_frac <= 1:                                         # No re-binning
        ntstable = tstable
        #new_time_bin = time_bin is there any iterative process?
        print('    No re-binning needed for this Light Curve. Input / Output time bins ratio : ', t_bin_frac)
    else:
        nout=int(nin/t_bin_frac)                     # Number of Out points. t_bin_frac == DTout / DTin ('TIMEDEL'). Re-binning factor
        #print('    Input / Output time bins / ratio : ', DTin, 'sec', ' / ',DTout, 'sec', ' / ', t_bin_frac)
        #print('    Input / Output dots number       : ', nin, ' / ',nout)
 
        i1 = np.arange(nout)*t_bin_frac
        i2 = np.concatenate(    (  i1[1:]-1 , np.array([nin-1])  )    )
 
        with warnings.catch_warnings():
          warnings.simplefilter("ignore", category=RuntimeWarning)              #
          if fracexp_chk == 0:          # TS has FRACEXP column => rate and error are weighted by FRACEXP
            nrate  = (    (rcum[i2]-rcum[i1]) / (fcum[i2]-fcum[i1])    )
            nerror = (np.sqrt(np.nan_to_num(ecum[i2]-ecum[i1]))/(fcum[i2]-fcum[i1]))

            nrate1  = (    (rcum1[i2]-rcum1[i1]) / (fcum[i2]-fcum[i1])    )
            nerror1 = (np.sqrt(np.nan_to_num(ecum1[i2]-ecum1[i1]))/(fcum[i2]-fcum[i1]))

            nrate2  = (    (rcum2[i2]-rcum2[i1]) / (fcum[i2]-fcum[i1])    )
            nerror2 = (np.sqrt(np.nan_to_num(ecum2[i2]-ecum2[i1]))/(fcum[i2]-fcum[i1]))

            nrate3  = (    (rcum3[i2]-rcum3[i1]) / (fcum[i2]-fcum[i1])    )
            nerror3 = (np.sqrt(np.nan_to_num(ecum3[i2]-ecum3[i1]))/(fcum[i2]-fcum[i1]))

            nrate4  = (    (rcum4[i2]-rcum4[i1]) / (fcum[i2]-fcum[i1])    )
            nerror4 = (np.sqrt(np.nan_to_num(ecum4[i2]-ecum4[i1]))/(fcum[i2]-fcum[i1]))

            nrate5  = (    (rcum5[i2]-rcum5[i1]) / (fcum[i2]-fcum[i1])    )
            nerror5 = (np.sqrt(np.nan_to_num(ecum5[i2]-ecum5[i1]))/(fcum[i2]-fcum[i1]))
 
            ntime= Time((tstable['TIME'][i2].value+tstable['TIME'][i1].value)/2., scale='tt', format=tstable['TIME'].format)
            nfrac=(fcum[i2]-fcum[i1])/(i2-i1)
 
            ntstable = Table([ntime, nrate, nerror, nrate1, nerror1,nrate2, nerror2,nrate3, nerror3,nrate4, nerror4,nrate5, nerror5, nfrac],\
                   names=['TIME','RATE','ERROR','RATE1','RATE1_ERR','RATE2','RATE2_ERR','RATE3','RATE3_ERR','RATE4','RATE4_ERR','RATE5','RATE5_ERR','FRACEXP'])
          else:                         # TS has NOT FRACEXP column => fracexp = 1
            nrate  = (    (rcum[i2]-rcum[i1]) / (fcum[i2]-fcum[i1])    )
            nerror = (np.sqrt(np.nan_to_num(ecum[i2]-ecum[i1]))/(fcum[i2]-fcum[i1]))

            nrate1  = (    (rcum1[i2]-rcum1[i1]) / (fcum[i2]-fcum[i1])    )
            nerror1 = (np.sqrt(np.nan_to_num(ecum1[i2]-ecum1[i1]))/(fcum[i2]-fcum[i1]))

            nrate2  = (    (rcum2[i2]-rcum2[i1]) / (fcum[i2]-fcum[i1])    )
            nerror2 = (np.sqrt(np.nan_to_num(ecum2[i2]-ecum2[i1]))/(fcum[i2]-fcum[i1]))

            nrate3  = (    (rcum3[i2]-rcum3[i1]) / (fcum[i2]-fcum[i1])    )
            nerror3 = (np.sqrt(np.nan_to_num(ecum3[i2]-ecum3[i1]))/(fcum[i2]-fcum[i1]))

            nrate4  = (    (rcum4[i2]-rcum4[i1]) / (fcum[i2]-fcum[i1])    )
            nerror4 = (np.sqrt(np.nan_to_num(ecum4[i2]-ecum4[i1]))/(fcum[i2]-fcum[i1]))

            nrate5  = (    (rcum5[i2]-rcum5[i1]) / (fcum[i2]-fcum[i1])    )
            nerror5 = (np.sqrt(np.nan_to_num(ecum5[i2]-ecum5[i1]))/(fcum[i2]-fcum[i1]))
 
            ntime= Time((tstable['TIME'][i2].value+tstable['TIME'][i1].value)/2., scale='tt', format=tstable['TIME'].format)
 
            ntstable = Table([ntime, nrate, nerror, nrate1, nerror1],\
                   names=['TIME','RATE','ERROR','RATE1','RATE1_ERR'])
 
 
    return ntstable


def table_to_timeseries(tbl, times, data):
    if any(isinstance(col, Time) for col in tbl.itercols()):
        for i, col in enumerate(tbl.itercols()): 
            if isinstance(col, Time):
                break
        col.format = 'jd'
        times.append(col)
        tbl.remove_column(tbl.colnames[i])
        
        for colname in tbl.colnames:
            unitstr = str(tbl[colname].unit)
            if unitstr == "e-/s":
                tbl[colname].unit = "electron/s"
            elif unitstr ==  "'electron'.s**-1":
                tbl[colname].unit = "electron/s"
            elif unitstr == "ct / s":
                tbl[colname].unit = "ct/s"
        
    data.append(tbl)
            
    return times, data


def ts_single_fits_reader(filename, times, data):
    """
    This serves as the FITS reader time series files within tsview. Based on kepler_fits_reader from astropy-timeseries. 

    This allows reading a supported FITS file using syntax such as::
    >>> from tsview.io import ts_fits_reader
    >>> time, data = ts_fits_reader('<name of fits file>')  # doctest: +SKIP

    Parameters
    ----------
    filename: `str`, `pathlib.Path`
        File to load.

    Returns
    -------
    list of `astropy.time.Time` and list of `astropy.table.QTable`
        time object and data for time series

    """
    
    # Get fits structure info in a table
    finfo = fits.info(filename, output=False) # list of tuples
    tinfo = Table(rows=finfo, names=('No.', 'Name', 'Ver', 'Type', 'Cards', 'Dimensions', 'Format', '')) # Convert into table
    tinfo.pprint()
    
    # Primary header without reading file
    hdr = fits.getheader(filename, extname="PRIMARY", ignore_missing_simple=True)
    
    i, hdr = locate_extension_with_keyword(tinfo, filename, 'telescop')
    # Get the telescope
    telescop = hdr['telescop'].lower()

    
    match telescop:
        case 'jwst':
            MULTI_KEYWORDS = ['FILTER']
            TIMESERIES_KEYWORDS = ['TIMESYS', 'TIMEUNIT', 
                    'EXPSTART', 'EXPMID', 'EXPEND', 'EFFEXPTM', #exposures related
                    'BARTDELT', 'BSTRTIME', 'BENDTIME', 'BMIDTIME', 'HELIDELT', 'HSTRTIME', 'HENDTIME', 'HMIDTIME', #exposures related on reference positions
                    'NINTS', 'EFFINTTM', 'INTSTART', 'INTEND' #integrations related
                    ]
            if hdr['TIMESYS'] != 'TDB':
                warnings.warn('Initial header Timesys is \'{0}\''.format(hdr['TIMESYS']))
            
            # if 'INT_TIMES' in tinfo['Name']: 
            #     times = []
            #     time = int_times_ext_to_time_arr(filename)
            #     times.append(time)
            # else:
            #     # TODO: Create the Time object from the exposure time only
            #     pass
            
            if 'EXTRACT1D' in tinfo['Name']:
                times, data = extract1d_to_timeseries(filename, times, data, hdr, tinfo, MULTI_KEYWORDS+TIMESERIES_KEYWORDS, array_int_times = False)
            elif 'SCI'in tinfo['Name']:
                # It is a rateints so the Type is ImageHDU and that is a cube to inspect in slices dimension?
                # TODO: Consult with Javier as this may be the handshake to cubeviewer.
                pass
        case 'xmm':
            
            MULTI_KEYWORDS = ['INSTRUME', 'FILTER', 'DSVAL4' ]
            hdr_multi = dict(filter(lambda i:i[0] in MULTI_KEYWORDS, hdr.items())) 
            
            header = fits.getheader(filename, extname="RATE", ignore_missing_simple=True)
            filter_keyword = header['FILTER']
            instrume_keyword = header['INSTRUME']
            tstart = header['TSTART']
            tstop = header['TSTOP']

            # Save the content in a variable
            header_keywords = {'FILTER': filter_keyword, 'INSTRUME': instrume_keyword, 'TSTART':tstart,'TSTOP':tstop}
            
            tbl = Table.read(filename, format='fits', astropy_native=True)
            
            # Bin data: tsr    = timeseries inicial # Tabla  ntsr = new timeseries, ya rebinned
            points_density = 512
            time_bin = tbl['TIMEDEL'][0] #binsize, same for all values in tseries
            tstart=header_keywords['TSTART']
            tstop=header_keywords['TSTOP']
            new_time_bin = (tstop-tstart) / ( points_density-1 )
            t_bin_frac=int(np.round(new_time_bin/time_bin))
        
            if t_bin_frac <= 1:
                rebin = False
            else:
                rebin = True
        
            if rebin:
                tbl = TSrebin(tbl, time_bin, new_time_bin)
            else:
                tbl = tbl
            # End bin data
            tbl.meta.update(hdr_multi)
            times, data = table_to_timeseries(tbl, times, data)
                
        case other:
            #times, data = [], []
            tbl = Table.read(filename, format='fits', astropy_native=True)
            times, data = table_to_timeseries(tbl, times, data)
        
            # raise NotImplementedError("{} is not implemented, only JWST is "
            #                       "supported through this reader".format(telescop))
    
    #TODO: Deal with masked data 
    
    #created for the xmm timeseries, as the scale in Terrestrial time gives error due to weakrefs
    times = sanitize_weakrefs(times)

    return times, data

@timer_func
def ts_fits_reader(name):
    
    times = []
    data = []

    if isinstance(name, list):
        for filename in name:
            times, data = ts_single_fits_reader(filename, times, data)
    else:
        times, data = ts_single_fits_reader(name, times, data)
            
    return times, data


# registry.register_reader('jwst.fits', TimeSeries, ts_fits_reader)

if __name__ == '__main__':
    
    filename = 'jw02783-o002_t001_miri_p750l-slitlessprism_x1dints.fits'
    time, data = ts_fits_reader(os.path.join(DATADIR, filename))
    
    filename = 'P0505720401PNX000SRCTSR800C.fits'
    time, data = ts_fits_reader(os.path.join(DATADIR, filename))
    
    
    import requests,io, gzip
    import tempfile
    
    # tap_server = 'https://jwst.esac.esa.int/server/tap/sync'
    # q = 'SELECT a.artifactid FROM  jwst.artifact AS a  WHERE a.uri LIKE \'%{0}%.fits\' AND a.uri LIKE \'%{1}%\' ORDER BY a.filename DESC'.format('x1dints','jw02783-o002_t001_miri_p750l-slitlessprism')
    # r = requests.post(tap_server, data = {
    #     'REQUEST':'doQuery',
    #     'LANG':'ADQL',
    #     'FORMAT':'json',
    #     'PHASE':'RUN',
    #     'QUERY':q
    #     })
    # artifact_id = r.json()['data'][0][0]
    # url_dp = 'https://jwst.esac.esa.int/server/data?ARTIFACTID={}&RETRIEVAL_TYPE=PRODUCT'.format(artifact_id)
    # # function dp_request
    # resp = requests.get(url_dp, timeout=1, verify=True)
    # #fits_content = io.BytesIO(gzip.decompress(resp.content))
    # # function get_data
    # fits_content = gzip.decompress(resp.content)
    # with tempfile.NamedTemporaryFile(delete=True) as fp:
    #     fp.write(fits_content)
    #     # Determine content-type in response (VOTable, FITS or csv)
    #     if resp.headers['content-type'] == 'application/fits':
    #         if os.path.exists(fp.name):
    #             time, data = ts_fits_reader(fp.name)
    
    url_dp = 'https://nxsa.esac.esa.int/nxsa-sl/servlet/data-action-aio?obsno={0}&sourceno={1}&extension=FTZ&level=PPS&instname=PN&name=SRCTSR&expflag=X'.format('0505720401', '00C')
    resp = requests.get(url_dp, timeout=1, verify=True)
    fits_content = gzip.decompress(resp.content)
    with tempfile.NamedTemporaryFile(delete=True) as fp:
        fp.write(fits_content)      
        # Determine content-type in response (VOTable, FITS or csv)
        if resp.headers['content-type'] in ['image/fits','application/fits', 'application/x-tar']:
            if os.path.exists(fp.name):
                time, data = ts_fits_reader(fp.name)
            pass
