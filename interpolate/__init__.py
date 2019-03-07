from __future__ import division, print_function, absolute_import


"""
individual module of interpolation modified from package wafo

Dedicated to WAFO

"""

#from . import misc

from . import polynomial

from . import interpolate


try:
    from wafo.version import version as __version__
except ImportError:
    __version__ = 'nobuilt'

from numpy.testing import Tester
test = Tester().test
