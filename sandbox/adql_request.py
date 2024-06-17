import requests,sys,io
from astropy.io.votable import parse as vo_parse

# The server we're talking to, replace as appropriate
url = 'https://gaia.ari.uni-heidelberg.de/tap/async'
q = 'SELECT TOP 1000 * FROM gaiadr3.gaia_source_lite'

# Send the request
r = requests.post(url, data = {
    'REQUEST':'doQuery',
    'LANG':'ADQL',
    'FORMAT':'votable/td',
    'PHASE':'RUN',
    'QUERY':q
    })

# Get the job ID; proper XML parsing is for chumps
jobid = r.text.split('<jobId>')[1].split('</jobId>')[0]
url2 = '{:}/{:}'.format(url,jobid)
# Wait for a while until the job is done

 # Get the data and put it in a byte stream
res = requests.get(url2+'/results/result')
#f = io.BytesIO(bytes(res.text,'utf-8'))
f = io.BytesIO(res.content)

# Parse the VO table and return the results
data = vo_parse(f).get_first_table()

#---------------------------------------------

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
#resp.content is bytes

from astropy.table import Table
import gzip

format = 'fits'

try:
    result = Table.read(io.BytesIO(gzip.decompress(resp.content)), format=format)
except OSError:
    # data is not a valid gzip file by BadGzipFile.
    result = Table.read(io.BytesIO(resp.content), format=format)
    pass

fits_content = io.BytesIO(gzip.decompress(resp.content))
from tsview.io.fits_parser import ts_fits_reader
from astropy.io import registry, fits
time, data = ts_fits_reader(fits_content)
#--------------------------------


def smart_divide(func):
    def inner(a, b):
        print("I am going to divide", a, "and", b)
        if b == 0:
            print("Whoops! cannot divide")
            return

        return func(a, b)
    return inner

@smart_divide
def divide(a, b):
    print(a/b)

divide(2,5)

divide(2,0)

f = dp_request(url_str, id)
https://gea.esac.esa.int/data-server/data?RETRIEVAL_TYPE=EPOCH_PHOTOMETRY&ID=Gaia+DR3+4111834567779557376&FORMAT=votable&RELEASE=Gaia+DR3&DATA_STRUCTURE=INDIVIDUAL

https://jwst.esac.esa.int/server/data?ARTIFACTID=2d5553ff-d07a-41c9-9656-3185c432210a&RETRIEVAL_TYPE=PRODUCT