import pytest
import os
import numpy as np
import json
from astropy import units as u
from astropy.io import fits
from astropy.table import Table, MaskedColumn


from tsview import DATADIR
from tsview.io.fits_parser import ts_fits_reader

def test_fits_parser():
    # test #1 loading x1dints
    '''    No.     Name    Ver      Type    Cards Dimensions                         Format                         col7
    int64    str9   int64    str11    int64   object                           str54                          str1
    ----- --------- ----- ----------- ----- ---------- ------------------------------------------------------ ----
        0   PRIMARY     1  PrimaryHDU   332         ()                                                            
        1 INT_TIMES     1 BinTableHDU    24   10R x 7C                                  [J, D, D, D, D, D, D]     
        2 EXTRACT1D     1 BinTableHDU    75 388R x 18C [D, D, D, D, D, D, D, D, D, D, D, J, D, D, D, D, D, D]     
        3 EXTRACT1D     2 BinTableHDU    75 388R x 18C [D, D, D, D, D, D, D, D, D, D, D, J, D, D, D, D, D, D]     
        4 EXTRACT1D     3 BinTableHDU    75 388R x 18C [D, D, D, D, D, D, D, D, D, D, D, J, D, D, D, D, D, D]     
        5 EXTRACT1D     4 BinTableHDU    75 388R x 18C [D, D, D, D, D, D, D, D, D, D, D, J, D, D, D, D, D, D]     
        6 EXTRACT1D     5 BinTableHDU    75 388R x 18C [D, D, D, D, D, D, D, D, D, D, D, J, D, D, D, D, D, D]     
        7 EXTRACT1D     6 BinTableHDU    75 388R x 18C [D, D, D, D, D, D, D, D, D, D, D, J, D, D, D, D, D, D]     
        8 EXTRACT1D     7 BinTableHDU    75 388R x 18C [D, D, D, D, D, D, D, D, D, D, D, J, D, D, D, D, D, D]     
        9 EXTRACT1D     8 BinTableHDU    75 388R x 18C [D, D, D, D, D, D, D, D, D, D, D, J, D, D, D, D, D, D]     
        10 EXTRACT1D     9 BinTableHDU    75 388R x 18C [D, D, D, D, D, D, D, D, D, D, D, J, D, D, D, D, D, D]     
        11 EXTRACT1D    10 BinTableHDU    75 388R x 18C [D, D, D, D, D, D, D, D, D, D, D, J, D, D, D, D, D, D]     
        12      ASDF     1 BinTableHDU    11    1R x 1C                                              [651792B]    
        ''' 
    rows_expected = 13
    id_expected = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
    extname_expected = ['PRIMARY', 'INT_TIMES', 'EXTRACT1D', 'EXTRACT1D', 'EXTRACT1D', 'EXTRACT1D', 'EXTRACT1D', 'EXTRACT1D', 'EXTRACT1D', 'EXTRACT1D', 'EXTRACT1D', 'EXTRACT1D','ASDF']
    extver_expected = [1, 1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 1]
    cards_expected = [332, 24, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 11]
    format_int_expected = ['J', 'D', 'D', 'D', 'D', 'D', 'D']
    format_extract1d_expected = ['D', 'D', 'D', 'D', 'D', 'D', 'D', 'D', 'D', 'D', 'D', 'J', 'D', 'D', 'D', 'D', 'D', 'D']
    format_asdf_expected = ['651792B']
    
    filename = 'jw02783-o002_t001_miri_p750l-slitlessprism_x1dints.fits'
    finfo = fits.info(os.path.join(DATADIR, filename), output=False) # list of tuples
    tinfo = Table(rows=finfo, names=('No.', 'Name', 'Ver', 'Type', 'Cards', 'Dimensions', 'Format', '')) # Convert into table
    
    assert len(tinfo) == rows_expected
    # integers and strings
    for meas, truth in zip(tinfo["No."], id_expected):
        assert truth == pytest.approx(meas)
    for meas, truth in zip(tinfo["Name"], extname_expected):
        assert truth == pytest.approx(meas)
    for meas, truth in zip(tinfo["Ver"], extver_expected):
        assert truth == pytest.approx(meas)
    for meas, truth in zip(tinfo["Cards"], cards_expected):
        assert truth == pytest.approx(meas)
    for meas, truth in zip(tinfo['Format'][1].tolist().strip('][').split(', '), format_int_expected):
        assert truth == pytest.approx(meas)
    for meas, truth in zip(tinfo['Format'][2].tolist().strip('][').split(', '), format_extract1d_expected):
        assert truth == pytest.approx(meas)
    for meas, truth in zip(tinfo['Format'][12].tolist().strip('][').split(', '), format_asdf_expected):
        assert truth == pytest.approx(meas)
    

    # test time values
    time_expected = ([59989.97705226, 59989.9772382 , 59989.97742414, 59989.97761009,
       59989.97779603, 59989.97798198, 59989.97816792, 59989.97835386,
       59989.97853981, 59989.97872575])
    rows_expected = 10
    time, data = ts_fits_reader(os.path.join(DATADIR, filename))
    assert len(time) == rows_expected
    assert time.format == 'mjd'
    assert time.scale == 'tdb'
    assert np.array(time_expected) == pytest.approx(time.value)
    # data values test
    assert isinstance(data, list)
    assert len(data) == rows_expected 
    assert isinstance(data[0], Table)
    assert isinstance(data[0]['FLUX'], MaskedColumn)
    # only data spot checks as the table has 388 values
    assert data[0]['FLUX'].value[0] == pytest.approx(0.15152622916204114)
    assert data[0]['FLUX'].mask[0] == False
    assert data[0]['FLUX'].unit == u.Unit("Jy")
 
 
     # test #3 loading rateints

if __name__ == '__main__':
    test_fits_parser()
    
