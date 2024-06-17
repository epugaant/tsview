from astropy.io import fits
from astropy.io.votable import parse as orig_parse
from astropy.io.votable.table import parse as vo_parse
from astropy.io.votable.tree import VOTableFile, Resource, Table as tree_tbl
from astropy.table import Table 
from astropy.io.fits.connect import read_table_fits
from astropy.io.votable.table import from_table
from astropy import units as u

import astropy.utils.data as aud
import io

# The JWST models:
from stdatamodels.jwst import datamodels

with aud.get_readable_fileobj('https://jwst.esac.esa.int/server/data?ARTIFACTID=f2c33ddd-9971-4598-833f-3b851dcb99c9&RETRIEVAL_TYPE=PRODUCT') as f:
    t = orig_parse(f)
    #character_as_bytes=False
    #convert_bytestring_to_unicode()[source]

filename = '/Users/epuga/ESDC/TSViz/data/jwst/test/jw02783-o002_t001_miri_p750l-slitlessprism_x1dints.fits'

with fits.open(filename) as hdulist:
    hdr = hdulist[0].header
    
# When you read as table, you have the attributes colnames, dtype, meta and info. 
# Units normally has to come from 
vot = Table.read(filename, format='fits')#needs an iterator and config
print(type(vot))
print(vot.info(out=None))

# Trying to cast it into a votable
cards = fits.info(filename, output=False)
nr_green_cards = len([s for s in cards if not any(x in s for x in ['PRIMARY', 'ASDF'])])
#I want to directly read a fits as VOTableFile, not as an astropy table
votable = VOTableFile()
# ...with one resource...
resource = Resource()
votable.resources.append(resource)
for i in range(1,nr_green_cards):
    #if isinstance(hdu, fits.BinTableHDU):
        # ... with one table
        t = Table.read(filename, format='fits', hdu=i)
        resource.tables.append(tree_tbl.from_table(votable, t.convert_bytestring_to_unicode()))
    # else:
    #     pass
        

    vot_0 = vo_parse(hdulist)
    print(type(vot_0))

    vot_1 = VOTableFile.parse(filename)
    print(type(vot_1))


#JWST x1dints.fits
# How to operate with datamodels
x1d = datamodels.open(filename)
x1d.info()
'''There is also a utility method for finding elements in the metadata schema. search_schema will search the schema for the given substring in metadata names as well as their documentation. The search is case-insensitive'''
x1d.search_schema('meta') #searches for anything that contains that string in the final name!!!!
x1d.search(value='UTC') #exact searches, outputs 
x1d.search(key='time_sys')
x1d.meta.instance
x1d.spec[0].instance

# Following recommendation of James Davies
# reading in as hdulist --> heavier
with fits.open(filename) as hdulist:
    x1d_ver = [hdu.ver for hdu in hdulist if hdu.name == "EXTRACT1D"]
for ver in x1d_ver:
    tab = Table.read(filename, hdu=("EXTRACT1D", ver))
    tab.meta.update(fits.getheader(filename, extname="PRIMARY"))
                  
# t = Table.read(filename, hdu=("EXTRACT1D", 1))
# t = Table.read(filename, format='fits', hdu=1)
             
ts_keys = ['TIMESYS', 'TIMEUNIT', 
                    'EXPSTART', 'EXPMID', 'EXPEND', 'EFFEXPTM', #exposures related
                    'BARTDELT', 'BSTRTIME', 'BENDTIME', 'BMIDTIME', 'HELIDELT', 'HSTRTIME', 'HENDTIME', 'HMIDTIME', #exposures related on reference positions
                    'NINTS', 'EFFINTTM', 'INTSTART', 'INTEND' #integrations related
                    ]

finfo = fits.info(filename, output=False) # list of tuples
tinfo = Table(rows=finfo, names=('No.', 'Name', 'Ver', 'Type', 'Cards', 'Dimensions', 'Format')) # Convert into table
hdr = fits.getheader(filename, extname="PRIMARY")


# Selection of cards related to time that we want to have in every table
hdr_time = dict(filter(lambda i:i[0] in ts_keys, hdr.items())) 
mask = tinfo['Name'] == 'EXTRACT1D'
extver = tinfo['Ver'][mask].data.tolist()
# extver = [ext[2] for ext in finfo for col in ext if  col  == 'EXTRACT1D'] # less readable
x1dints = []
for ver in extver:
   tab = Table.read(filename, hdu=("EXTRACT1D", ver))
   # propagate time information from the primary to all extract1d
   tab.meta.update(hdr_time)
   x1dints.append(t)

import matplotlib
import matplotlib.pyplot as plt
print(len(x1dints))

x1dfig, x1dax = plt.subplots(figsize=[12,4])

# Plot multiple spectra
for i in range(len(x1dints)):
    x1dax.plot(x1dints[i]['WAVELENGTH'], x1dints[i]['FLUX'])
    
x1dax.set_title('Extracted spectra (per integration)')
x1dax.set_xlabel('wavelength ($\mu$m)')
x1dax.set_ylabel('flux (Jy)')
plt.show()

# Time object
# Is there an INT_TIMES extension?
t = Table.read(filename, hdu="INT_TIMES")   

from astropy.time import Time, TimeDelta
# scale could come from hdr['TIMESYS'] is 'UTC'
t1 = Time(t['int_mid_MJD_UTC'].data, format='mjd', scale='utc')
t2 = Time(t['int_mid_BJD_TDB'].data, format='mjd', scale='tdb')

#checking some things
# t1.utc - t2.utc
dt = TimeDelta(hdr['BARTDELT'] * u.Unit(hdr['TIMEUNIT'])) 

format = 'BJD'.lower()
if any(format in col.lower() for col in  t.colnames):   
    sub_cols = [col for col in t.colnames if format in col.lower()]
    int_mid = [col for col in sub_cols if 'mid' in col]
    if len(int_mid) == 1:
        scale = int_mid[0].split('_')[-1].lower()
        print(format, scale)


if any( col == 'INT_TIMES' for ext in finfo for col in ext):
if any( col == 'SCI' for ext in finfo for col in ext):
    
if 'INT_TIMES' in tinfo['Name']: # It means one can explore per integration
    #
elif:

mask = tinfo['Name'] == 'EXTRACT1D'
tinfo['Ver'][mask]