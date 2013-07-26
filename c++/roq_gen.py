#! /usr/bin/env python2.7

# Import some things from Python 3
from __future__ import print_function, division
from time import gmtime, strftime
from subprocess import call

# Import some os stuff
stamp = "2013-07-26"
stamp_in = "2013-07-26"

roq_gen = "bin/roq_generator"
roq_gen_rc = "rc/roq-IPTA-like"
roq_gen_param_rc = "rc/roq-IPTA-like-params"
error = "1e-33"

# Make the binaries, just in case
call(["make"])

for j in xrange(1,2):
  roq_gen_param_rc_tmp = roq_gen_param_rc + "-" + str(j)
  for i in xrange(1,7):
    pulsarNumber = str(6*i)
    stamp_tmp = stamp + "-" + str(i) + "-param-" + str(j)
    stamp_tmp_in = stamp_in + "-" + str(i)
    call([roq_gen, roq_gen_rc, roq_gen_param_rc_tmp, stamp_tmp, stamp_tmp_in, error ])
