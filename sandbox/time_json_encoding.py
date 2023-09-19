# https://github.com/astropy/astropy/issues/9774

from astropy.time import Time
import json
start = Time.now()
json.dumps(start)

from astropy.utils.misc import JsonCustomEncoder
json.dumps(start, cls=JsonCustomEncoder)

class JsonCustomEncoder(json.JSONEncoder):
    """Support for data types that JSON default encoder
    does not do.

    This includes:

        * Numpy array or number
        * Complex number
        * Set
        * Bytes
        * astropy.UnitBase
        * astropy.Quantity

    Examples
    --------
    >>> import json
    >>> import numpy as np
    >>> from astropy.utils.misc import JsonCustomEncoder
    >>> json.dumps(np.arange(3), cls=JsonCustomEncoder)
    '[0, 1, 2]'

    """

    def default(self, obj):
        from astropy import units as u
        import numpy as np
        if isinstance(obj, u.Quantity):
            return dict(value=obj.value, unit=obj.unit.to_string())
        if isinstance(obj, (np.number, np.ndarray)):
            return obj.tolist()
        elif isinstance(obj, complex):
            return [obj.real, obj.imag]
        elif isinstance(obj, set):
            return list(obj)
        elif isinstance(obj, bytes):  # pragma: py3
            return obj.decode()
        elif isinstance(obj, (u.UnitBase, u.FunctionUnitBase)):
            if obj == u.dimensionless_unscaled:
                obj = 'dimensionless_unit'
            else:
                return obj.to_string()

        return json.JSONEncoder.default(self, obj)
    
#example for SkyCoord
import json
from astropy.coordinates import SkyCoord
from astropy.utils.misc import JsonCustomEncoder
sc = SkyCoord(1, 2, unit="deg")
json.dumps(sc.info._represent_as_dict(), cls=JsonCustomEncoder)
'{"ra": {"value": 1.0, "unit": "deg"}, "dec": {"value": 2.0, "unit": "deg"}, "representation_type": "spherical", "frame": "icrs"}'
[obj.__dict__ for obj in list_name]

import json
import numpy as np
from astropy.time import Time
from astropy.utils.misc import JsonCustomEncoder
t = Time.now()
json.dumps(t.info._represent_as_dict(), cls=JsonCustomEncoder)
'{"jd1": 2460206.0, "jd2": -0.16213025939814818, "format": "datetime", "scale": "utc", "precision": 3, "in_subfmt": "*", "out_subfmt": "*"}'
t = Time(['2001:020', '2001:040', '2001:060', '2001:080'], out_subfmt='date')
t[2] = np.ma.masked
json.dumps(t.info._represent_as_dict(), cls=JsonCustomEncoder)
'{"jd1": [2451930.0, 2451950.0, 2451970.0, 2451990.0], "jd2": [-0.5, -0.5, -0.5, -0.5], "format": "yday", "scale": "utc", "precision": 3, "in_subfmt": "*", "out_subfmt": "date"}'
t = [Time.now(), Time.now(), Time.now()]
json.dumps([obj.info._represent_as_dict() for obj in t], cls=JsonCustomEncoder)

'''From https://stackoverflow.com/questions/21411497/flask-jsonify-a-list-of-objects'''
app = Flask(__name__)
app.json_encoder = MyJSONEncoder
# and just pass in your list directly to jsonify():
return jsonify(my_list_of_obj)



import json
from json import JSONEncoder
import numpy as np

#all numpy types encoder
class NumpyArrayEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        else:
            return super(NumpyArrayEncoder, self).default(obj)

# Serialization
numPyData = {"id": 25, "floatSample": np.float32(1.2), "intSample":np.int32(42), "arangeSample": np.arange(12)}
encodedNumpyData = json.dumps(numPyData, cls=NumpyArrayEncoder)  # use dump() to write array into file
print("Printing JSON serialized NumPy array")