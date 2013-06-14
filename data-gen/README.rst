=======================
A simple data generator
=======================

This is a simple data generator using the Ellis et al. paper [Ellis2010a]_. You need
Cython in order to run the code. Also, reffer to the paper by Corbin and Cornish
[Corbin2008a]_ for some detail on the formulae derivations. There is also the Jaranowski
et al. [Jaranowski2008a]_ paper for explanation of intrinsic and extrinsic parameter
classification.

The original paper is slightly unclear about some of the parameters
The parameter names:

* :math:`\zeta := M^{5/3}/D`
* :math:`\cos(\iota)` - orbital inclination (? is it the cosine, or the angle itself?)
* :math:`\psi` - polarisation angle of the GW
* :math:`t` - time
* :math:`\theta, \phi` - spherical polar angles for a vector pointing from Solar System
  Baricentre to the GW source
* :math:`\theta_p, \phi_p` - angular pulsar coordinates
  This is an added parameter by me, which I will probably substitute to something
  better
* :math:`\omega_0` - Initial frequency
* :math:`L` - Distance to the pulsar
* :math:`D` - Distance to the Blackhole binary
* :math:`M` - the chirp mass of the binary
  .. math:: M^{5/3} = \frac{m_1 m_2}{(m_1 + m_2)^{1/3}}
* :math:`\omega(t)` - orbital frequency
* :math:`\Phi(t), \Phi_0` - orbital phase

 .. [Ellis2010a] http://arxiv.org/abs/1204.4218
 .. [Corbin2008a] http://arxiv.org/abs/1008.1782
 .. [Jaranowski2008a] http://arxiv.org/abs/gr-qc/9804014
