#! /usr/bin/python2.7

from __future__ import print_function, division

import pyximport; pyximport.install()
import signal_model as sm
import numpy as np

# I generate N pulsars in the solid angle defined by \theta and \phi randomly
# scattered accross some length
N = 25
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

# Time schedule of the pulsar timing:
schedule = np.random.rand(N) * t_interval_min

# Define our pulsar grid
pulsars = sm.PulsarGrid (N, PulsarCoordinateRanges)
a = []
for i in range(N):
  a += [ np.array([]) ]

# Define the constants for use in the data generation
t = 5
M = 5
D = 10e9*M
iota = 2
Phi0 = 1
psi = 2
theta = 0.5
phi = 0.2
omega0 = 1e-8
CONTINUE = True

# Generate the actual data.
while CONTINUE:
  CONTINUE = False

  for j in range(N):
    t = schedule[j]

    if t < t_final:
      L = pulsars.getLength(j)
      u_p = pulsars.getUnitVector(j)
      a[j] = np.append(a[j], sm.residual (t, M, D, iota, Phi0, psi, theta, phi, omega0, L, u_p))
      schedule[j] += np.random.rand()*abs(t_interval_max - t_interval_min) + \
                     t_interval_min
    t = schedule[j]
    # Do not stop generating data if there is at least one pulsar, which had not "gone"
    # in to the future
    if t < t_final:
      CONTINUE = True

# Contract the data into a one vector
data = np.array([])
for i in a:
  data = np.append(data, i)

# Save the data in a gz format. Numpy load txt understands gzipped files
# transparently
np.savetxt("pulsardata.txt.gz", data)
