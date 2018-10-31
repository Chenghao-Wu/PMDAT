#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
dir_path = os.path.dirname(os.path.realpath(__file__))
root_dir = dir_path.replace("/util", '', 1)
sys.path.insert(0, root_dir)

from reader.Reader import *
from reader.LAMMPSReader import *
from reader.XYZReader import *

from analysis_list.build_list import *
from analysis_list.multi_body import *
from system import *
from application import *

from analysis.analysis import *
from analysis.mean_squared_displacement import *
from analysis.radial_distribution_function import *
from analysis.wave_vectors import *
from analysis.intermediate_scattering_function import *
from analysis.end_end_distance import *
from analysis.gyration_tensor import *
from analysis.vector_autocorrelation_function import *
from analysis.bond_autocorrelation_function import *
from analysis.order_parameter import *

print("Soft Matter Molecular Simulation Analysis Tool(SMMSAT) version "+System.Version)