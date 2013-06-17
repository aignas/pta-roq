#! /usr/bin/python2.7

import pyximport; pyximport.install()
import signal_model as sm
import numpy as np

N = 8
PulsarCoordinateRanges = np.array([ [10,40],[0.2,1.],[-0.5,0.5]])
t_max = 10.
t_N = 255

# Define our pulsar grid
pulsars = sm.PulsarGrid (N, PulsarCoordinateRanges)
a = np.zeros((t_N,N))

#def residual (t, zeta, iota, Phi0, psi, theta, phi, omega0, L, u_p):
t = 5
M = 1
D = 20
iota = 2
Phi0 = 1
psi = 2
theta = 0.5
phi = 0.2
omega0 = 0.00004

for i in range(t_N):
  t = t_max/t_N * i
  for j in range(N):
    L = pulsars.getLength(j)
    u_p = pulsars.getUnitVector(j)
    a[i,j] = sm.residual (t, M, D, iota, Phi0, psi, theta, phi, omega0, L, u_p)

print a
