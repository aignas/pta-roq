#! /usr/bin/env python2.7

# Import some things from Python 3
from __future__ import print_function, division
from time import gmtime, strftime
from subprocess import call
import sys

# Import some os stuff
#stamp = strftime("%F", gmtime())
stamp_dat = "2013-07-29"
stamp_roq = "2013-07-29"
pulsarNumber="36"
t_final="5"
dt_min="2"
dt_max="2"

# ROQ Error
error = "1e-33"

roq_gen = "bin/roq_generator"
roq_gen_rc = "rc/roq-IPTA-like"
roq_gen_param_rc = roq_gen_rc + "-params"
data_gen = "bin/data_generator"
data_gen_rc = "rc/data-IPTA-like"
data_gen_source_rc = data_gen_rc + "-sources-1"

# Make the binaries, just in case
call(["make"])

def gen_data():
  for i in xrange(1,7):
    pulsarNumber = str(6*i)
    stamp_tmp = stamp_dat + "-" + str(i)
    call([data_gen, data_gen_rc, stamp_tmp, data_gen_source_rc, pulsarNumber, t_final, dt_min, dt_max])

def gen_roq():
  for j in xrange(1,2):
    roq_gen_param_rc_tmp = roq_gen_param_rc + "-" + str(j)
    for i in xrange(1,7):
      pulsarNumber = str(6*i)
      stamp_tmp_in = stamp_dat + "-" + str(i)
      stamp_tmp = stamp_roq + "-" + str(i) + "-param-" + str(j)
      call([roq_gen, roq_gen_rc, roq_gen_param_rc_tmp, stamp_tmp, stamp_tmp_in, error ])

gen_data()
gen_roq()
