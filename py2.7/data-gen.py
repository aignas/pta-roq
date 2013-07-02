#! /usr/bin/env python2.7

from __future__ import print_function, division

import pyximport; pyximport.install()
import signal_model as sm
import numpy as np

# I generate N pulsars in the solid angle defined by \theta and \phi randomly
# scattered accross some length
N = 36
PulsarCoordinateRanges = np.array([ [10,40],[0, np.pi],[-np.pi, np.pi]])
pulsarNoise = np.array([0.003, 0, 0])

# Some time variables so that it is easier to read the code
week = 7*24*3600
yr = 52 * week

# Stop the generation here
t_final = 5*yr
# Interval when we do not take any measurements
dt_min = 2 * week
dt_max = 2 * week
# Max number of timings
N_t_max = 1
N_t_min = 1
#

# Initial times for each pulsar measurement:
time_init = np.random.rand(N) * 0

# Define our pulsar grid and the sources
pulsars = sm.PulsarGrid (N, PulsarCoordinateRanges, pulsarNoise)
# Define the parameters for each source. The sources matrix is of the form:
# sources = np.array([
#   M, D, iota, Phi0, psi, theta, phi, omega0,         1st source
#   M, D, iota, Phi0, psi, theta, phi, omega0,         2nd source
#   ])

sources = np.array([
   5, 10e9*5, 2, 1, 2, 0.5, 0.2, 2*np.pi*1e-8  # 1st source
  ])

# Generate the actual data.
print("Generate a measurement schedule")
schedule = sm.genSchedule(time_init, t_final, dt_min, dt_max)
print("Generate the data according to the schedule")
data = sm.dataGeneration(schedule, sources, pulsars)

# Save the data in a gz format. Numpy load txt understands gzipped files
# transparently
print("Save the generated data and the schedule")
np.savetxt("pulsardata.txt.gz", data)
np.savetxt("pulsarschedule.txt.gz", schedule)
