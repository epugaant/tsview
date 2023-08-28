# Licensed under a 3-clause BSD style license - see LICENSE.rst
import warnings

import numpy as np

from astropy.io import registry, fits
from astropy.table import Table
from astropy.time import Time, TimeDelta

__all__ = ["ts_fits_reader"]


def ts_fits_reader(filename):
    """
    This serves as the FITS reader time series files within tsview. Based on kepler_fits_reader from astropy-timeseries. 

    This allows reading a supported FITS file using syntax such as::
    >>> from tsview.io import ts_fits_reader
    >>> timeseries = TimeSeries.read('<name of fits file>', format='kepler.fits')  # doctest: +SKIP

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
    tinfo = Table(rows=finfo, names=('No.', 'Name', 'Ver', 'Type', 'Cards', 'Dimensions', 'Format')) # Convert into table
    
    # Primary header without reading file
    hdr = fits.getheader(filename, extname="PRIMARY")

    
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
        raise NotImplementedError("{} is not implemented, only KEPLER or TESS are "
                                  "supported through this reader".format(hdulist[0].header['telescop']))

    # Time
    if 'INT_TIMES' in tinfo['Name']: 
        t = Table.read(filename, hdu="INT_TIMES") 
        
        format = 'BJD'.lower()
        if any(format in col.lower() for col in t.colnames):   
            sub_cols = [col for col in t.colnames if format in col.lower()]
            int_mid = [col for col in sub_cols if 'mid' in col]
            if len(int_mid) == 1:
                scale = int_mid[0].split('_')[-1].lower()
        
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
    if 'SCI'in tinfo['Name']:
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
    tab.remove_column('time')
'''

    return time, data

registry.register_reader('kepler.fits', TimeSeries, ts_fits_reader)
registry.register_reader('kepler.fits', TimeSeries, kepler_fits_reader)
registry.register_reader('tess.fits', TimeSeries, kepler_fits_reader)