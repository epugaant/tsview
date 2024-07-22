DATA_ACCESS = [
    {'mission': 'gaia',   
     'server_url': 'https://gea.esac.esa.int',
     'endpoint': 'data-server/data',
     'query': 'RETRIEVAL_TYPE=EPOCH_PHOTOMETRY&ID={0}&FORMAT=votable&RELEASE=Gaia+DR3&DATA_STRUCTURE=INDIVIDUAL',
     # example url: 'https://gea.esac.esa.int/data-server/data?RETRIEVAL_TYPE=EPOCH_PHOTOMETRY&ID={0}&FORMAT=votable&RELEASE=Gaia+DR3&DATA_STRUCTURE=INDIVIDUAL'
    },
    {'mission': 'jwst',
     'server_url': 'https://jwst.esac.esa.int',
     'endpoint': 'server/data',
     'query': 'ARTIFACTID={0}&RETRIEVAL_TYPE=PRODUCT',
     # example url: 'https://jwst.esac.esa.int/server/data?ARTIFACTID={0}&RETRIEVAL_TYPE=PRODUCT'
     'adql_info': {
        'url': 'https://jwst.esac.esa.int/server/tap/sync',
        'REQUEST':'doQuery',
        'LANG':'ADQL',
        'FORMAT':'json',
        'PHASE':'RUN',
        'QUERY':'SELECT a.artifactid FROM  jwst.artifact AS a WHERE a.uri LIKE \'%{0}%.fits\' AND a.uri LIKE \'%{1}%\' ORDER BY a.filename DESC'
        },
    },
    {'mission': 'xmm-epic',   
     'server_url': 'https://nxsa.esac.esa.int',
     'endpoint': 'nxsa-sl/servlet/data-action-aio',
     'query': 'obsno={0}&sourceno={1:03X}&extension=FTZ&level=PPS&instname=PN&name=SRCTSR&expflag=X',
     'query_replacement': 'param[1:-4], int(param[-4:].strip("0"))'                                                                                                         
     }
    ]

# class Aio:
#     def __init__(self, mission):
#         self.mission = mission
#         self.server_url = None
#         self.endpoint = None
#         self.query = None
#         self.adql_info = None
#         self.query_replacement = None
#         for data in DATA_ACCESS:
#             if data['mission'] == mission:
#                 self.server_url = data['server_url']
#                 self.endpoint = data['endpoint']
#                 self.query = data['query']
#                 self.adql_info = data.get('adql_info', None)
#                 self.query_replacement = data.get('query_replacement', None)
#                 break
#         if self.server_url is None:
#             raise ValueError('Mission not found: {}'.format(mission))
            
#     def format_query(self,param):
#         if self.query_replacement:
#             return self.query.format(*eval(self.query_replacement))
#         else:
#             return self.query.format(param)
#     def adql_request(self, obsID, prodType):
#         '''Function to generate request to TAP Server using ADQL query to get data product id, e.g. artifact id'''
    
#         r = requests.post(self.adql_info['url'], data = {
#         'REQUEST': self.adql_info['REQUEST'],
#         'LANG': self.adql_info['LANG'],
#         'FORMAT': self.adql_info['FORMAT'],
#         'PHASE': self.adql_info['PHASE'],
#         'QUERY': self.adql_info['QUERY'].format(prodType, obsID)
#         })
#         #print(r.json())
#         id = r.json()['data'][0][0]
#         return id
    
#     def dp_request(url):
#         '''Funtion to request data product through server´s response. Outputs a list of Times and a list of Tables.'''

#         #Get response from API
#         try:
#             resp = requests.get(url, timeout=20, verify=True)
#             resp.raise_for_status()
#         except requests.exceptions.HTTPError as errh:
#             print("HTTP Error")
#             print(errh.args[0])
#         except requests.exceptions.ReadTimeout as errrt:
#             print("Time out")
#         except requests.exceptions.ConnectionError as conerr:
#             print("Connection error")
#         except requests.exceptions.RequestException as errex:
#             print("Exception request")
        
#         return resp
        
# Aio('xmm-epic').format_query('106945101010025')
# Aio('gaia').format_query('Gaia+DR3+4111834567779557376')
# Aio('jwst').format_query('jw02783-o002_t001_miri_p750l-slitlessprism')
        
        