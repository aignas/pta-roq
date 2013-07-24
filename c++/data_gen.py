#! /usr/bin/env python2.7

# Import some things from Python 3
from __future__ import print_function, division
from time import gmtime, strftime
from subprocess import call

# Import some os stuff
stamp = strftime("%F", gmtime())
pulsarNumber="36"
t_final="5"
dt_min="2"
dt_max="2"

data_gen = "bin/data_generator"
data_gen_rc = "rc/data-IPTA-like"

# Make the binaries, just in case
call(["make"])

for i in xrange(1,7):
  pulsarNumber = str(6*i)
  stamp_tmp = stamp + "-" + str(i)
  call([data_gen, data_gen_rc, stamp_tmp, pulsarNumber, t_final, dt_min, dt_max])
