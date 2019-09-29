#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
dir_path = os.path.dirname(os.path.realpath(__file__))
root_dir = dir_path.replace("/util", '', 1)
sys.path.insert(0, root_dir)



from src.reader.Reader import *
from src.reader.LAMMPSReader import *
from src.reader.XYZReader import *

from src.utility.system import *
from src.utility.application import *

from src.list.multibody_list import *
from src.list.atom_list import *

from src.method.analysis import *
from src.method.wave_vectors import *

from src.method.dynamics.mean_squared_displacement import *
from src.method.dynamics.intermediate_scattering_function import *
from src.method.dynamics.vector_autocorrelation_function import *
from src.method.dynamics.zeroth_modulus import *


from src.method.statics.radial_distribution_function import *
from src.method.statics.end_end_distance import *
from src.method.statics.gyration_tensor import *
from src.method.statics.mean_squeared_internal_distance import *
from src.method.statics.structure_factor import *


#from pyZ1.pyZ1 import *



print("Soft Matter Molecular Simulation Analysis Tool(SMMSAT) version "+System.Version)