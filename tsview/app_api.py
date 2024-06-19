from flask import request, make_response
from flask import Flask
from flask_caching import Cache
import requests
import io
import os
from urllib.parse import urljoin
import gzip
import tempfile

import astropy.units as u
from astropy.utils.misc import JsonCustomEncoder

# config file stored in configfile
from tsview.access_config import DATA_ACCESS

from tsview.io.fits_parser import ts_fits_reader
from tsview.io.vo_parser import ts_votable_reader
from tsview.aggregate.data_processing import DataProcess, equivalent_units

TRUE_STR = ("yes", "true", "t", "1")
FALSE_STR = ("no", "false", "f", "0")

cacheConfig = {
    "DEBUG": True,          # some Flask specific configs
    "CACHE_TYPE": "SimpleCache",  # Flask-Caching related configs
    "CACHE_DEFAULT_TIMEOUT": 300
}

app = Flask(__name__)
app.json_encoder = JsonCustomEncoder
app.config.from_mapping(cacheConfig)
cache = Cache(app)

def str2bool(val):
    val = val.lower()
    if val in TRUE_STR:
        return True
    elif val in FALSE_STR:
        return False
    else:
        raise ValueError("invalid truth value %r" % (val,))

def mission_config(data_access, mission):
    '''Return the single dictionary corresponding to the mission'''
    # check if mission is in schema if not raise exception
    if not any(d['mission'] == mission for d in data_access): 
        raise NotImplementedError("Mission {} schema for time series is not implemented yet".format(mission))
    mission_access = [d for d in data_access if d['mission'] == mission][0]
    return mission_access
    
def build_dp_url(server_url, endpoint_path, query_params):
    '''Construct string with url with placeholders for parameters'''
    base_url = urljoin(server_url, endpoint_path)
    return base_url  + "?" + query_params

def adql_request(adql_info, obsID, prodType):
    '''Function to generate request to TAP Server using ADQL query to get data product id, e.g. artifact id'''
    
    r = requests.post(adql_info['url'], data = {
    'REQUEST': adql_info['REQUEST'],
    'LANG':adql_info['LANG'],
    'FORMAT':adql_info['FORMAT'],
    'PHASE':adql_info['PHASE'],
    'QUERY':adql_info['QUERY'].format(prodType, obsID)
    })
    print(r.json())
    id = r.json()['data'][0][0]
    return id
    
def dp_request(url_str, id):
    '''Funtion to request data product through server´s response. Outputs a list of Times and a list of Tables.'''

    url = url_str.format(id)
    #Get response from API
    try:
        resp = requests.get(url, timeout=20, verify=True)
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
    
    return resp
    
def get_data(resp):    
    '''Function to get response content in a byte array and then in time and data, depending of the content type'''
    
    try:
        res = gzip.decompress(resp.content)
    except OSError:
    # data is not a valid gzip file by BadGzipFile.
        res = resp.content
        pass
    #If it is encoded in Bytes, then we need to wrap the response in a seekable file-like object for any Table.read to work
    if isinstance(resp.content, bytes):
        f = io.BytesIO(res)
    else:
        f = res
    print(type(f)) # '_io.BytesIO'
    
    # Write the response to a temporary file
    # Delete needs to be False, otherwise the file isn't findable on Windows machines
    with tempfile.NamedTemporaryFile(delete=False) as fp:
        fp.write(resp.content)
        # Determine content-type in response (VOTable, FITS or csv)
        if resp.headers['content-type'] == 'application/fits':
            if os.path.exists(fp.name):
                time, data = ts_fits_reader(fp.name)
        elif resp.headers['content-type'] == 'application/x-votable+xml': 
            if os.path.exists(fp.name):
                time, data = ts_votable_reader(fp.name)
    return time, data


def _build_cors_preflight_response():
    response = make_response()
    response.headers.add("Access-Control-Allow-Origin", "*")
    response.headers.add('Access-Control-Allow-Headers', "*")
    response.headers.add('Access-Control-Allow-Methods', "*")
    return response

def _create_cors_response():
    response = make_response()
    response.headers.add("Access-Control-Allow-Origin", "*")
    return response

#The route function tells the application which API endpoint URL should call the associated function. 
@app.route('/ts/v1', methods=['GET', 'OPTIONS'])
def get_data_to_plot():
    if request.method == "OPTIONS": # CORS preflight
            return _build_cors_preflight_response()

    # Check if an ID was provided as part of the URL.
    # If ID is provided, assign it to a variable.
    # If no ID is provided, display an error in the browser.
    query_parameters = request.args

    mission = query_parameters.get('mission')
    system = query_parameters.get('system')
    time_view = query_parameters.get('timeView')

    try:
        sourceID = request.args['sourceID']
    except:
        sourceID = ''
    try:
        obsID = request.args['obsID']
    except:
        obsID = ''
    target_time_unit = query_parameters.get('target_time_unit')
    target_time_scale = query_parameters.get('target_time_scale')
    target_flux_unit = query_parameters.get('target_flux_unit')
    

    key = ''.join(filter(None, (mission, sourceID, obsID)))# mission + sourceID + obsID 
    print(key)
    cached_data = cache.get(key)
    if cached_data is None:
        print('Data not cached yet. Requesting data from server.')
        if mission:
            #select the mission access details
            mission_access = mission_config(DATA_ACCESS, mission)
            server_url = mission_access['server_url']
            endpoint_path = mission_access['endpoint']
            query_params = mission_access['query']

            #check if there is an adql query necessary
            if 'adql_info' in mission_access:
                adql_info = mission_access['adql_info']
                
                prodType = query_parameters['prodType']
                if 'obsID' in query_parameters:
                    id = adql_request(adql_info, obsID, prodType)
                else:
                    raise Exception('obsID is mandatory for mission {}'.format(mission))
            else:
                adql_info = None
                if 'sourceID' in query_parameters:
                    id = sourceID
                else:
                    raise Exception('sourceID is mandatory for mission {}'.format(mission))

            # Construct the url string to request https server data
            url_str = build_dp_url(server_url, endpoint_path, query_params)
            
            # Data request 
            resp = dp_request(url_str, id)
        
            # Data extraction
            time, data = get_data(resp)
            cache.set(key, [time, data])
        else:
            time = data = None
            raise NotImplementedError("Error: No valid mission field provided. Please specify an mission registered in config.")
    else:
        print('Data already cached. Using cached data.')
        [time, data] = cached_data
    if not isinstance(time, list): 
        time = [time]
    else:
        pass
    
    #  # Create a list for the objects
    # times_list = [obj.info._represent_as_dict() for obj in time]
    # times_final = json.dumps(times_list, cls=JsonCustomEncoder)
   
    # data_list = [obj.to_pandas().to_json() for obj in data]
    # data_final = json.dumps(data_list)

    # return {'time': times_final, 'times_len': len(times_list) , 'data': data_final, 'data_len': len(data_list)}
    # #return jsonify({'times': [obj.info._represent_as_dict() for obj in time]}) 
    if system:
        d = DataProcess(mission, time, data, system)
    else:
        d = DataProcess(mission, time, data) 

        
    if target_time_unit:
        d.convert_time(target_time_unit)# mjd
    if target_time_scale:
        d.convert_time_scale(target_time_scale)# tcb
    #WARNING: ad-hoc code to skip instrumental units for Gaia
    if target_flux_unit and target_flux_unit.strip() not in ("electron/s"):
        d.convert_flux(u.Unit(target_flux_unit))# u.mJy
    
    #WARNING: functional code to get alt_data_unit for instrumental units
    if not d.alt_data_unit:
        d.alt_data_unit =  equivalent_units(u.mJy.to_string()) + ['mJy']

    response = _create_cors_response()
    response.data = d.to_plotly(time=str2bool(time_view))
    return response

@app.route('/ts/v1/modifytime', methods=['GET', 'OPTIONS'])
def convert_time():
    if request.method == "OPTIONS": # CORS preflight
            return _build_cors_preflight_response()

    query_parameters = request.args

    mission = query_parameters.get('mission')
    system = query_parameters.get('system')
    
    try:
        sourceID = request.args['sourceID']
    except:
        sourceID = ''
    try:
        obsID = request.args['obsID']
    except:
        obsID = ''    
    target_time_unit = query_parameters.get('target_time_unit')
    target_time_scale = query_parameters.get('target_time_scale')
    
    key = ''.join(filter(None, (mission, sourceID, obsID)))# mission + sourceID + obsID 
    print(key)
    cached_data = cache.get(key)
    if cached_data is not None:
        [time, data] = cached_data
        if not isinstance(time, list): 
            time = [time]
        if system:
            d = DataProcess(mission, time, data, system)
        else:
            d = DataProcess(mission, time, data) 

        if target_time_unit:
            d.convert_time(target_time_unit)# mjd
        if target_time_scale:
            d.convert_time_scale(target_time_scale)# tcb

        response = _create_cors_response()
        response.data = d.to_plotly()
        return response
    else:
        print('You need to run the caching endpoint before')
        return
    
@app.route('/ts/v1/modifyflux', methods=['GET', 'OPTIONS'])
def convert_flux():
    if request.method == "OPTIONS": # CORS preflight
            return _build_cors_preflight_response()

    query_parameters = request.args

    mission = query_parameters.get('mission')
    system = query_parameters.get('system')
    time_view = query_parameters.get('timeView')
    
    try:
        sourceID = request.args['sourceID']
    except:
        sourceID = ''
    try:
        obsID = request.args['obsID']
    except:
        obsID = ''    
    target_flux_unit = query_parameters.get('target_flux_unit')
    
    
    key = ''.join(filter(None, (mission, sourceID, obsID)))# mission + sourceID + obsID 
    print(key)
    cached_data = cache.get(key)
    if cached_data is not None:
        [time, data] = cached_data
        if not isinstance(time, list): 
            time = [time]
        if system:
            d = DataProcess(mission, time, data, system)
        else:
            d = DataProcess(mission, time, data) 

        if target_flux_unit:
            d.convert_flux(u.Unit(target_flux_unit))# mjd 
        response = _create_cors_response()
        response.data = d.to_plotly(time=str2bool(time_view))
        return response
    else:
        print('You need to run the caching endpoint before')
        return
   
if __name__ == '__main__':
    app.run(host='0.0.0.0',port = 8000, threaded = True, debug = True)
    
    
#http://0.0.0.0:8000

'''# User parameters, via API
mission = 'gaia'
sourceID = 'Gaia+DR3+4111834567779557376'
native units
http://0.0.0.0:8000/ts/v1?mission=gaia&sourceID=Gaia+DR3+4111834567779557376
calibrated units
http://0.0.0.0:8000/ts/v1?mission=gaia&sourceID=Gaia+DR3+4111834567779557376&target_time_unit=mjd&target_flux_unit=mJy

http://0.0.0.0:8000/ts/v1?mission=gaia&sourceID=Gaia+DR3+4111834567779557376&target_time_unit=mjd&target_flux_unit=mJy&system=VEGAMAG
http://0.0.0.0:8000/ts/v1?mission=gaia&sourceID=Gaia+DR3+4111834567779557376&target_time_unit=mjd&system=VEGAMAG

http://0.0.0.0:8000/ts/v1/modifytime?mission=gaia&sourceID=Gaia+DR3+4111834567779557376&target_time_unit=jd&system=VEGAMAG
http://0.0.0.0:8000/ts/v1/modifyflux?mission=gaia&sourceID=Gaia+DR3+4111834567779557376&target_flux_unit=Jy&system=VEGAMAG

http://0.0.0.0:8000/ts/v1?mission=gaia&sourceID=Gaia%20DR3%203643046442207745280&target_time_unit=jd&target_flux_unit=electron%20/%20s&timeView=true
http://0.0.0.0:8000/ts/v1?mission=gaia&sourceID=Gaia+DR3+4111834567779557376&target_time_unit=mjd&target_time_scale=tdb&timeView=True

mission = 'jwst'
sourceID = None,
obsID = 'jw02783-o002_t001_miri_p750l-slitlessprism'
prodType = 'x1dints' 
http://0.0.0.0:8000/ts/v1?mission=jwst&obsID=jw02783-o002_t001_miri_p750l-slitlessprism&prodType=x1dints
http://0.0.0.0:8000/ts/v1?mission=jwst&obsID=jw02783-o002_t001_miri_p750l-slitlessprism&prodType=x1dints&target_time_unit=jd&target_flux_unit=Jy
http://0.0.0.0:8000/ts/v1/modifyflux?mission=jwst&obsID=jw02783-o002_t001_miri_p750l-slitlessprism&target_flux_unit=Jy&timeView=True
'''

app.config["SECRET_KEY"] = "any random string"
app.config["CACHE_TYPE"] = "SimpleCache" # better not use this type w. gunicorn
cache = Cache(app)

@app.route("/login", methods=["GET", "POST"])
def login():
    if request.method == "POST":
        cache.set("username", request.form["username"])
    # to get value use cache.get("username")