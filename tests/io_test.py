import pytest
import os
import numpy as np
import io
from astropy import units as u
from astropy.io import fits
from astropy.table import Table, MaskedColumn
from astropy.io import votable
from astropy.time import Time


from tsview import DATADIR
from tsview.io.fits_parser import ts_fits_reader
from tsview.io.vo_parser import locate_time_column, col_to_jd, ts_votable_reader
from astropy.io.votable.table import is_votable

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
    

    # test time values in mjd
    time_expected = ([59989.97705226, 59989.9772382 , 59989.97742414, 59989.97761009,
       59989.97779603, 59989.97798198, 59989.97816792, 59989.97835386,
       59989.97853981, 59989.97872575])
    rows_expected = 10
    time, data = ts_fits_reader(os.path.join(DATADIR, filename))
    assert len(time) == rows_expected
    assert time[0].format == 'mjd'
    assert time[0].scale == 'tdb'
    assert np.array(time_expected) == pytest.approx([t.value for t in time])
    # data values test
    assert isinstance(data, list)
    assert len(data) == rows_expected 
    assert isinstance(data[0], Table)
    assert isinstance(data[0]['FLUX'], MaskedColumn)
    # only data spot checks as the table has 388 values
    assert data[0]['FLUX'].value[0] == pytest.approx(0.15152622916204114)
    assert data[0]['FLUX'].mask[0] == False
    assert data[0]['FLUX'].unit == u.Unit("Jy")
 
    #test chandra time values in jd
    filename = 'chandra_time.fits'
    time_expected = ([2457414.26033393, 2457414.2603339287])
    rows_expected = 1
    time, data = ts_fits_reader(os.path.join(DATADIR, filename))
    assert len(time) == rows_expected
    assert time[0].format == 'mjd'
    assert time[0].scale == 'tt'
    assert np.array(time_expected) == pytest.approx(time[0].jd)
    
     # test #3 loading rateints

def test_vo_parser():
    # time system with literal xml example
    TIME_SERIES_LITERALS = b'''<?xml version="1.0" encoding="UTF-8"?>
    <VOTABLE version="1.4" xmlns="http://www.ivoa.net/xml/VOTable/v1.3">
    <RESOURCE>
        <COOSYS ID="system" epoch="J2015.5" system="ICRS"/>
        <TIMESYS ID="time_frame" refposition="BARYCENTER" timeorigin="2455197.5" timescale="TCB"/>
        <TABLE name="ts_data">
        <FIELD datatype="double" name="obs_time" ucd="time.epoch" unit="d" ref="time_frame"/>
        <FIELD datatype="float" name="flux" ucd="phot.flux;em.opt.V" unit="s**-1"/>
        <FIELD datatype="float" name="mag" ucd="phot.mag;em.opt.V" unit="mag"/>
        <FIELD datatype="float" name="flux_error" ucd="stat.error;phot.flux;em.opt.V" 
            unit="s**-1"/>
        <PARAM datatype="double" name="ra" ucd="pos.eq.ra" value="45.7164887146879" ref="system"/>
        <PARAM datatype="double" name="dec" ucd="pos.eq.dec" value="1.18583048057467" ref="system"/>
        <DATA>
            <TABLEDATA>
            <TR>
                <TD>1821.2846388435</TD>
                <TD>168.358</TD>
                <TD>20.12281560517953</TD>
                <TD>8.71437</TD>
            </TR>
            </TABLEDATA>
        </DATA>
        </TABLE>
    </RESOURCE>
    </VOTABLE>'''
    rows_expected = 1
    time_expected = ([1821.2846388435])
    time_fto_jd_expected = ([2457018.7846388435]) # time from time origin
    time_fto_mjd_expected = ([57018.2846388435])
    
    vot = votable.parse(io.BytesIO(TIME_SERIES_LITERALS))
    timesystems = vot.resources[0].time_systems[0] # only element in a HomogeneousList, that contains dictionaries
    time_series = vot.resources[0].tables[0]

    # Time system comes in an `astropy.utils.collections.HomogeneousList` of TimeSys objects with attributes
    # ID, refposition, timeorigin and timescale, accessible via timesystems.timescale
    
    # Convert TimeSys attributes to dictionary
    times_meta = {key: getattr(timesystems, key) for key in timesystems._attr_list} 
    # just extract the times MaskedColumn
    t = locate_time_column(time_series, times_meta)

    time = Time(t, scale=timesystems.timescale.lower(),format='jd')
    print(dir(times_meta))
    # this tests the Column and MaskedColumn functionality of to_jd function
    time_fto_jd = col_to_jd(t, times_meta)
    
    assert len(time) == rows_expected
    assert time.format == 'jd'
    assert time.scale == 'tcb'
    assert np.array(time_expected) == pytest.approx(time.value)
    
    assert len(time_fto_jd) == rows_expected
    assert time_fto_jd.format == 'jd'
    assert time_fto_jd.scale == 'tcb'
    assert np.array(time_fto_jd_expected) == pytest.approx(time_fto_jd.value) 
    assert np.array(time_fto_mjd_expected) == pytest.approx(time_fto_jd.mjd) 
    
    # test with multi in xml and global ts_votable_reader
    filename = 'votable_multi_example.xml'
    vot = votable.parse(os.path.join(DATADIR, filename))
    rows_expected = 5
    time_fto0_expected = [2450000.00032, 2450000.00048, 2450000.00064, 2450000.0008,  2450000.00048,
    2450000.00048, 2450000.00048, 2450000.00048, 2450000.00048, 2450000.00112,
    2450000.00112, 2450000.00112, 2450000.00112, 2450000.00112, 2450000.00112, 2450000.00144]
    
    time_fto1_expected = [2450000.0014401, 2450000.0014401, 2450000.0014401, 2450000.0014401,
    2450000.0014401,  2450000.00160008, 2450000.00176006, 2450000.00192004,
    2450000.00208002, 2450000.00208002, 2450000.00208002, 2450000.00208002,
    2450000.00208002, 2450000.00208002, 2450000.00239998, 2450000.00239998]
    
    time_fto2_expected = [2450000.00764061, 2450000.00764061, 2450000.00764061, 2450000.00764061,
    2450000.00780061, 2450000.00796062, 2450000.00812062, 2450000.00844064,
    2450000.00844064, 2450000.00860064, 2450000.00876064, 2450000.00860064,
    2450000.00876064, 2450000.00892065, 2450000.00892065, 2450000.00572057]
    
    time_fto3_expected = [2450000.00029799, 2450000.00061799, 2450000.00125799, 2450000.00157799,
    2450000.00173799, 2450000.00189799, 2450000.00205799, 2450000.00173799,
    2450000.00189799, 2450000.00285799, 2450000.00269799, 2450000.00253799,
    2450000.00237799, 2450000.00317799, 2450000.00333799, 2450000.00349799]
    
    time_fto4_expected = [2450000.0098859 , 2450000.00972519, 2450000.00956448,
    2450000.00956448, 2450000.00956448, 2450000.00956448,
    2450000.01020367, 2450000.01020367, 2450000.01036438,
    2450000.01052509, 2450000.0106858 , 2450000.01036438,
    2450000.01052509, 2450000.0106858 , 2450000.0106858 ,
    2450000.01052509, 2450000.01036438, 2450000.01020367]

    times, data = ts_votable_reader(vot)
    
    assert len(times) == rows_expected
    assert np.array(time_fto0_expected) == pytest.approx(times[0].value)
    assert times[0].scale == 'tcb'
    assert np.array(time_fto1_expected) == pytest.approx(times[1].value)
    assert times[1].scale == 'tcb'
    assert np.array(time_fto2_expected) == pytest.approx(times[2].value)
    assert times[2].scale == 'tt'
    assert np.array(time_fto3_expected) == pytest.approx(times[3].value)
    assert times[3].scale == 'utc'
    assert np.array(time_fto4_expected) == pytest.approx(times[4].value)
    assert times[4].scale == 'tcg'
   
    # test ts_fits_reader with Gaia file
    pass

if __name__ == '__main__':
    test_fits_parser()
    test_vo_parser()
    
