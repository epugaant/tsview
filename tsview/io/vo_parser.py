from astropy.io.votable.table import parse as vo_parse
from astropy.io.votable.table import is_votable
from flask import jsonify

from astropy.table import Table
from astropy.table import Column, MaskedColumn
from astropy.time import Time, TimeDelta

import datetime
import io
import numpy as np

from astropy.io import votable
from astropy import coordinates
from astropy import table
from astropy import units
from astropy.time import Time

MAGIC_TIMEORIGINS = {
  "MJD-origin": 2400000.5,
  "JD-origin": 0,
  }

TIMESYS_SCALES_TO_ASTROPY = {
    'TAI': 'tai',
    'TT': 'tt',
    'UT': 'ut1',
    'UTC': 'utc',
#    'GPS': is handled separately, needs to be researched
    'TCG': 'tcg',
    'TCB': 'tcb',
    'TDB': 'tdb',
}
     
def locate_time_column(table, times_meta):
    """returns a numpy array of times in the VOTable table together with
    its TIMESYS metadata.

    The time is the first column referencing something in timesystems.
    """
    for col_index, field in enumerate(table.fields):
        if field.ref in times_meta.values(): # can also be getattr(column, 'ref')
            break
    else:
        raise Exception("No column referencing a TIMESYS found")

    return table.to_table().columns[col_index]

def to_jd(times, times_meta):
    """returns Time given in JD, for Time instances and table columns, considering the timeorigin.

    This is where we use FIELD metadata, TIMESYS/@timeorigin and 
    TIMESYS/@timescale
    """
    try:
        astropy_scale =TIMESYS_SCALES_TO_ASTROPY[times_meta["timescale"]]
    except KeyError:
        raise Exception("Unsupported timescale: %s"%times_meta["timescale"])

    if "timeorigin" not in times_meta:
        
        if isinstance(times, (Column, MaskedColumn)):
            # no timeorigin given: it's some sort of civic year
            if times.dtype=='O': 
                # These are strings and hence presumably timestamps
                # (we should rather inspect xtype here)
                # This should really work like this:
                # times = Time(times, format="isot", scale=astropy_scale)
                # but because of breakage in Debian stretch we do it half-manually:
                # Note: the following was modified to work on files, not on literals
                # times = Time(
                #     [datetime.datetime.strptime(
                #             v.split(".")[0], "%Y-%m-%dT%H:%M:%S")
                #         for v in times],
                #     scale=astropy_scale)
                times = Time(times.astype(str), format="isot", scale=astropy_scale)
                times.format = 'jd'

            else:
                # in VOTable, these must be julian or besselian years
                if times.unit not in ["yr", "byr"]:
                    raise Exception("Floats without timeorigin only allowed when"
                        " unit is year.")
                else:    
                    times = Time(times, 
                        format={"yr": "jyear", "byr": "byear"}[str(times.unit)],
                        scale=astropy_scale)
                    times.format = 'jd'
        elif isinstance(times, Time):
            times.format = 'jd'

    else:
        if isinstance(times, (Column, MaskedColumn)):
            times = Time(times.to(units.d)+times_meta["timeorigin"]*units.d,
            format="jd", scale=astropy_scale)
        elif isinstance(times, Time):
            if not isinstance(times_meta["timeorigin"], str):
                times = times + TimeDelta(times_meta["timeorigin"], 
                                             format="jd", scale=astropy_scale) 
            else:
                times = times + TimeDelta(MAGIC_TIMEORIGINS.get(times_meta["timeorigin"]), 
                                             format="jd", scale=astropy_scale)
        else:
            pass
        times.format = 'jd'

    return times



def ts_votable_reader(vot):
    '''Function to parse votable filename or VOTable into Time and data Table objects. The function deals with multiple resources returning two lists '''
    
    try:
        # If vot is a filename, it has a read attribute
        if is_votable(vot): # With Gaia we have 'astropy.io.votable.tree.VOTableFile'
            vot = vo_parse(vot)
    except:
        pass
        
    #timesystem = vot.resources[0].time_systems
    # We assume that all tables in resource are the same
    tbl = vot.get_first_table().to_table(use_names_over_ids=True)
    # tbl = vot.resources[0].tables[0]

    # Print the table column information
    tinfo = tbl.info(out=None)
    tinfo.pprint()

    # Create an empty list for our results
    times = []
    data = []
    
    # Flatten down all tables
    for resource in vot.resources:
        for i, table in enumerate(resource.tables):
            timesystems = resource.time_systems[i] # only element in a HomogeneousList, that contains dictionaries
            print(timesystems)
            times_meta = {key: getattr(timesystems, key) for key in timesystems._attr_list if getattr(timesystems, key) is not None} 
            
            astropy_scale = TIMESYS_SCALES_TO_ASTROPY[times_meta["timescale"]]

            # just extract the times MaskedColumn
            t = locate_time_column(table, times_meta)
            print(t.unit)

            if t.dtype=='O': 
                tt = Time(t.astype(str), format="isot", scale=astropy_scale)
            else:
                #{"yr": "jyear", "byr": "byear", "d": "jd"}[str(t.unit)]
                tt = Time(t, format='jd', scale=astropy_scale)
            tt = to_jd(tt, times_meta)

            tbl = table.to_table()                 
            if times_meta['ID'] in tbl.colnames:
                tbl.remove_column(times_meta['ID'])
            
            times.append(tt)
            data.append(tbl)

    return times, data

if __name__ == '__main__':
    
    import requests, os, io
    import tempfile
    
    import json
    from astropy.utils.misc import JsonCustomEncoder
    import numpy as np
    
    url = 'https://gea.esac.esa.int/data-server/data?RETRIEVAL_TYPE=EPOCH_PHOTOMETRY&ID={0}&FORMAT=votable&RELEASE=Gaia+DR3&DATA_STRUCTURE=INDIVIDUAL'.format('Gaia+DR3+4111834567779557376')
    resp = requests.get(url, timeout=1, verify=True)
    res = resp.content
    if isinstance(resp.content, bytes):
        f = io.BytesIO(res)
    else:
        f = res
    print(type(f)) # '_io.BytesIO'
    vot = vo_parse(f)
    # t = Table.read('aj285677t3_votable.xml')
    # vot.to_xml('test.xml')
    # if is_votable(res):
    with tempfile.NamedTemporaryFile(delete=True) as fp:
        fp.write(res)
        # Determine content-type in response (VOTable, FITS or csv)
        if resp.headers['content-type'] == 'application/x-votable+xml':
            if os.path.exists(fp.name):
                times, data = ts_votable_reader(fp.name)
                if not isinstance(times, list): 
                    t = [times]
                else:
                    pass
            result = json.dumps([obj.info._represent_as_dict() for obj in times], cls=JsonCustomEncoder)
            #this part only concerns the API
            #print(json.dumps(times, cls=JsonCustomEncoder)) #does not work
            # if times[0].masked:
            #     times[0]._time.jd
            #     times[0].value.data #ndarray witn nans or datetime if it only one point
            # if data[0].masked:
            #     result = data[0].to_pandas().to_json() #immune to mask

            pass
                
