import io
import datetime

from astropy.io import votable
from astropy import table
from astropy.table import Column, MaskedColumn
from astropy import units
from astropy.time import Time, TimeDelta

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
    """returns Time instances given in JD.

    This is where we use FIELD metadata, TIMESYS/@timeorigin and 
    TIMESYS/@timescale
    """
    try:
        astropy_scale =TIMESYS_SCALES_TO_ASTROPY[times_meta["timescale"]]
    except KeyError:
        raise Exception("Unsupported timescale: %s"%times_meta["timescale"])

    if "timeorigin" not in times_meta:
        # no timeorigin given: it's some sort of civic year
        if times.dtype=='O': 
            # These are strings and hence presumably timestamps
            # (we should rather inspect xtype here)
            # This should really work like this:
            # times = Time(times, format="isot", scale=astropy_scale)
            # but because of breakage in Debian stretch we do it half-manually:
            # times = Time(
            #     [datetime.datetime.strptime(
            #             v.split(b".")[0].decode("ascii"), "%Y-%m-%dT%H:%M:%S")
            #         for v in times],
            #     scale=astropy_scale)
            times.format = 'jd'

        else:
            # in VOTable, these must be julian or besselian years
            if times.unit not in ["yr", "byr"]:
                raise Exception("Floats without timeorigin only allowed when"
                    " unit is year.")
            times = Time(times, 
                format={"yr": "jyear", "byr": "byear"}[str(times.unit)],
                scale=astropy_scale)
            times.format = 'jd'

    else:
      if isinstance(times, (Column, MaskedColumn)):
        times = Time(
          times.to(units.d)+times_meta["timeorigin"]*units.d,
          format="jd", scale=astropy_scale)
      elif isinstance(times, Time):
        times = times + TimeDelta(times_meta["timeorigin"], format="jd",
                     scale=astropy_scale) 
      else:
        pass
    
    return times

if __name__ == '__main__':

    vot = votable.parse(io.BytesIO(TIME_SERIES_LITERALS))
    print(len(vot.resources))
    timesystems = vot.resources[0].time_systems[0] # only element in a HomogeneousList, that contains dictionaries
    time_series = vot.resources[0].tables[0]

    # Time system comes in an `astropy.utils.collections.HomogeneousList` of TimeSys objects
    # timesystems.ID
    # timesystems.refposition
    # timesystems.timeorigin
    # timesystems.timescale
    
    # Convert TimeSys attributes to dictionary
    times_meta = {key: getattr(timesystems, key) for key in timesystems._attr_list} 
    # just extract the times MaskedColumn
    t = locate_time_column(time_series, times_meta)

    times = Time(t, scale=timesystems.timescale.lower(),format='jd')
    print(dir(times_meta))
    print(times.value)
    times2 = to_jd(times, times_meta)
    print(times2.value)
    #times = to_tcb_barycenter(times, times_meta, obs_location)

    apy_table = time_series.to_table()

