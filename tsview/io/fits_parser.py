# Licensed under a 3-clause BSD style license - see LICENSE.rst
import os
import warnings
import numpy as np

from astropy.io import registry, fits
from astropy.table import Table
from astropy.time import Time, TimeDelta

from tsview import DATADIR

__all__ = ["ts_fits_reader"]


def ts_fits_reader(filename):
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
    `astropy.time.Time` and list of `astropy.table.QTable`
        time object and data for time series

    """
    # Get fits structure info in a table
    finfo = fits.info(filename, output=False) # list of tuples
    tinfo = Table(rows=finfo, names=('No.', 'Name', 'Ver', 'Type', 'Cards', 'Dimensions', 'Format', '')) # Convert into table
    
    # Primary header without reading file
    hdr = fits.getheader(filename, extname="PRIMARY", ignore_missing_simple=True)
    
    # Get the telescope
    telescop = hdr['telescop'].lower()

    if telescop == 'jwst':
        ts_keys = ['TIMESYS', 'TIMEUNIT', 
                    'EXPSTART', 'EXPMID', 'EXPEND', 'EFFEXPTM', #exposures related
                    'BARTDELT', 'BSTRTIME', 'BENDTIME', 'BMIDTIME', 'HELIDELT', 'HSTRTIME', 'HENDTIME', 'HMIDTIME', #exposures related on reference positions
                    'NINTS', 'EFFINTTM', 'INTSTART', 'INTEND' #integrations related
                    ]
        if hdr['TIMESYS'] != 'TDB':
            warnings.warn('Initial header Timesys is \'{0}\''.format(hdr['TIMESYS']))
        #hdu = hdulist['LIGHTCURVE']
    # elif telescop == 'kepler':
    #     hdu = hdulist[1]
    else:
        raise NotImplementedError("{} is not implemented, only JWST is "
                                  "supported through this reader".format(telescop))

    # Time
    if 'INT_TIMES' in tinfo['Name']: 
        t = Table.read(filename, hdu="INT_TIMES") 
        
        format = 'BJD'.lower()
        if any(format in col.lower() for col in t.colnames):   
            sub_cols = [col for col in t.colnames if format in col.lower()]
            int_mid = [col for col in sub_cols if 'mid' in col]
            if len(int_mid) == 1:
                
                scale = int_mid[0].split('_')[-1].lower()
                format = 'mjd'
        
        time = Time(t[int_mid[0]].data, format=format, scale=scale)

    else:
        # TODO: Create the Time object from the exposure time only
        pass
    
    #Data 
    # x1dints   
    if 'EXTRACT1D' in tinfo['Name']:
        # subdictionary related to time that we want to have in meta for every table
        hdr_time = dict(filter(lambda i:i[0] in ts_keys, hdr.items())) 
        mask = tinfo['Name'] == 'EXTRACT1D'
        extver = tinfo['Ver'][mask].data.tolist()
        data = []
        for ver in extver:
            tab = Table.read(filename, hdu=("EXTRACT1D", ver))
            # propagate time information from the primary to all extract1d
            tab.meta.update(hdr_time)
            data.append(tab)
    # rateints
    elif 'SCI'in tinfo['Name']:
        # It is a rateints so the Type is ImageHDU and that is a cube to inspect in slices dimension?
        # TODO: Consult with Javier
        pass
    '''    if hdu.header['EXTVER'] > 1:
        raise NotImplementedError("Support for {0} v{1} files not yet "
                                  "implemented".format(hdu.header['TELESCOP'], hdu.header['EXTVER']))

    # Check time scale
    if hdu.header['TIMESYS'] != 'TDB':
        raise NotImplementedError("Support for {0} time scale not yet "
                                  "implemented in {1} reader".format(hdu.header['TIMESYS'], hdu.header['TELESCOP']))

    tab = Table.read(hdu, format='fits')

    # Some KEPLER files have a T column instead of TIME.
    if "T" in tab.colnames:
        tab.rename_column("T", "TIME")

    for colname in tab.colnames:
        # Fix units
        if tab[colname].unit == 'e-/s':
            tab[colname].unit = 'electron/s'
        if tab[colname].unit == 'pixels':
            tab[colname].unit = 'pixel'

        # Rename columns to lowercase
        tab.rename_column(colname, colname.lower())

    # Filter out NaN rows
    nans = np.isnan(tab['time'].data)
    if np.any(nans):
        warnings.warn('Ignoring {0} rows with NaN times'.format(np.sum(nans)))
    tab = tab[~nans]

    # Time column is dependent on source and we correct it here
    reference_date = Time(hdu.header['BJDREFI'], hdu.header['BJDREFF'],
                          scale=hdu.header['TIMESYS'].lower(), format='jd')
    time = reference_date + TimeDelta(tab['time'].data)
    time.format = 'isot'

    # Remove original time column
    tab.remove_column('time')'''

    return time, data


# registry.register_reader('jwst.fits', TimeSeries, ts_fits_reader)

if __name__ == '__main__':
    
    from datetime import datetime
    start_time = datetime.now()
    
    filename = 'jw02783-o002_t001_miri_p750l-slitlessprism_x1dints.fits'
    time, data = ts_fits_reader(os.path.join(DATADIR, filename))
    
    end_time = datetime.now()
    delta = end_time - start_time
    print('Time to execute app is {} seconds'.format(delta.total_seconds()))
    
    
    import requests,io, gzip
    import tempfile
    
    start_time = datetime.now()
    
    tap_server = 'https://jwst.esac.esa.int/server/tap/sync'
    q = 'SELECT a.artifactid FROM  jwst.artifact AS a  WHERE a.uri LIKE \'%{0}%.fits\' AND a.uri LIKE \'%{1}%\' ORDER BY a.filename DESC'.format('x1dints','jw02783-o002_t001_miri_p750l-slitlessprism')
    r = requests.post(tap_server, data = {
        'REQUEST':'doQuery',
        'LANG':'ADQL',
        'FORMAT':'json',
        'PHASE':'RUN',
        'QUERY':q
        })
    artifact_id = r.json()['data'][0][0]
    url_dp = 'https://jwst.esac.esa.int/server/data?ARTIFACTID={}&RETRIEVAL_TYPE=PRODUCT'.format(artifact_id)
    resp = requests.get(url_dp, timeout=1, verify=True)
    #fits_content = io.BytesIO(gzip.decompress(resp.content))
    fits_content = gzip.decompress(resp.content)
    with tempfile.NamedTemporaryFile(delete=True) as fp:
        fp.write(fits_content)
        # Determine content-type in response (VOTable, FITS or csv)
        if resp.headers['content-type'] == 'application/fits':
            if os.path.exists(fp.name):
                time, data = ts_fits_reader(fp.name)
    
    end_time = datetime.now()
    delta = end_time - start_time
    print('Time to execute app is {} seconds'.format(delta.total_seconds()))