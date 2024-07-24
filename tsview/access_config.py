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
     'query': 'obsno={0}&sourceno={1:03X}&extension=FTZ&level=PPS&instname=PN,M1,M2&name=SRCTSR&expflag=X',
     'query_replacement': 'param[1:-4], int(param[-4:].strip("0"))'                                                                                                         
     }
    ]

        
        