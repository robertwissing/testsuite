#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 25 10:42:27 2024

@author: robertwi
"""

import numpy as np

# Precision setting (switch between float32 and float64)
PRECISION = np.float64  # Change to np.float64 when needed

# Boundary type constants
PERIODIC = 0
REFLECTING = 1
LATTICE = 2
ALPHA = 0.1
BETA = 2.0
ETA = 0.01
PARTICLE_ARRAY_SIZE=21

# Utility functions for creating arrays
def float_array(values):
    return np.array(values, dtype=PRECISION)

def scalar(value):
    return PRECISION(value)