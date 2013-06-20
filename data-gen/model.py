#! /usr/bin/python2.7

from __future__ import print_function, division

import pyximport; pyximport.install()
import signal_model as sm
import numpy as np

# I generate N pulsars in the solid angle defined by \theta and \phi randomly
# scattered accross some length
N = 8
PulsarCoordinateRanges = np.array([ [10,40],[0.2,1.],[-0.5,0.5]])
t = 0.
# Stop the generation here
t_final = 20*365.25*24*3600
# Interval when we do not take any measurements
t_interval_min = 0.9*7*24*3600  # 0.9 week gap is minimum
t_interval_max = 1.3*7*24*3600  # 1.3 week gap is maximum
# Max number of timings
N_t_max = 1
N_t_min = 1
#

# Time schedule of the pulsar timing:
schedule = np.random.rand(N) * 0

# Define our pulsar grid and the sources
pulsars = sm.PulsarGrid (N, PulsarCoordinateRanges, pulsarNoise)
sources = np.array([ sm.Source( M=5, D=10e9*5, iota=2, Phi0=1, psi=2,
                                theta=0.5, phi=0.2, omega0=1e-8)
                    ])

# Generate the actual data.
data, dates = sm.dataGeneration(schedule, sources, pulsars, t_final, 
                          t_interval_min, t_interval_max)

# Save the data in a gz format. Numpy load txt understands gzipped files
# transparently
np.savetxt("../data-crunch/pulsardata.txt.gz", data)
np.savetxt("../data-crunch/pulsarschedule.txt.gz", dates)
