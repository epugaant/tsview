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
    times_new = [Time(time.value, scale='tt', format=time.format) for time in times if time.scale == 'tt']
    return times_new

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
        case other:
            #times, data = [], []
            tbl = Table.read(filename, format='fits', astropy_native=True)
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
    
    filename = 'P0505720401PNS001SRCTSR800C.FTZ'
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
