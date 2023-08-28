from flask import Flask, jsonify
import requests
import io
from astropy.io.votable.table import parse as vo_parse
from astropy.io.votable.table import is_votable
#from lxml import etree as ElementTree

app = Flask(__name__)

#The route function tells the application which URL should call the associated function.
@app.route('/', methods=['GET'])
def index():
    
    url = 'https://gea.esac.esa.int/data-server/data?RETRIEVAL_TYPE=EPOCH_PHOTOMETRY&ID=Gaia+DR3+4111834567779557376&FORMAT=votable&RELEASE=Gaia+DR3&DATA_STRUCTURE=INDIVIDUAL'
    
    #Get response from API
    try:
        resp = requests.get(url, timeout=1, verify=True)
        resp.raise_for_status()
    except requests.exceptions.HTTPError as errh:
        print("HTTP Error")
        print(errh.args[0])
    except requests.exceptions.ReadTimeout as errrt:
        print("Time out")
    except requests.exceptions.ConnectionError as conerr:
        print("Connection error")
    except requests.exceptions.RequestException as errex:
        print("Exception request")
        
    #What type of content do we have?
    print(type(resp.content)) # 'bytes'
    #If it is encoded in Bytes, then we need to wrap the data in a seekable file-like object for votable.parse to work
    if isinstance(resp.content, bytes):
        f = io.BytesIO(resp.content)
    else:
        f = resp.content
    print(type(f)) # '_io.BytesIO'
    
    # Parse the already VO table (only in the Gaia Case)
    
    if is_votable(f): # With Gaia we have 'astropy.io.votable.tree.VOTableFile'
        vot = vo_parse(f)
        # for info in vot.iter_info():
        #     print(info) 
        tbl = vot.get_first_table().to_table(use_names_over_ids=True)
    
        # Print the table column information
        tinfo = tbl.info(out=None)
        tinfo.pprint()
        
        # Present the data via basic json file display
        
        # Create an empty list for our results
        results = []
        
        for table in vot.iter_tables():
            tbl = table.to_table()
            d = {
                'source_id': tbl['source_id'],
                'transit_id': tbl['transit_id'],
                'band': tbl['band'],
                'time': tbl['time'].data.tolist(),
                'time_unit': tbl['time'].unit,
                'flux': tbl['flux'].data.tolist(),
                'flux_unit': tbl['flux'].unit,
                'flux_error': tbl['flux_error'].data.tolist(),
                'flux_error_unit': tbl['flux_error'].unit
            }
            results.append(d)

        
        # Use the jsonify function from Flask to convert our list of
        # Python dictionaries to the JSON format.
        return jsonify(results)

if __name__ == '__main__':
    app.run(host='0.0.0.0',port = 8000, threaded = True, debug = True)
    
#http://0.0.0.0:8000
