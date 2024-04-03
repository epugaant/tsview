'''tsview is a web application to visualize time series data products served by ESA Archives APIs'''
import os

__version__ = '0.1'


# set Python env variable to keep track of example data dir (e.g orbitize)
data_dir = os.path.dirname(__file__)
DATADIR = os.path.join(data_dir, 'example_data/')

#methods that the interpreter will use in the wild card import statement
# __all__ = ["attitude",
#            "utils"]

from .app_api import app