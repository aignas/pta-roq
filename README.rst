=======
Pulsar Timing Array - Reduced Order Methods
=======

This is my summer project on Pulsar Timming Array data fitting using Reduced Order
Quadrature methods

The structure of the project:

1. This is a simple data generator using the Ellis et al. paper [Ellis2010a]_. You need
   Cython in order to run the code. Also, reffer to the paper by Corbin and Cornish
   [Corbin2008a]_ for some detail on the formulae derivations. There is also the
   Jaranowski et al. [Jaranowski2008a]_ paper for explanation of intrinsic and extrinsic
   parameter classification.

   TODO:

 a) add red noise

2. Generation of the Reduced Basis for the model.

3. Analysis of the usual problem, so that proper comparisons can be done.

   Notes on the original method:

   There is also a speed-up version of the [Ellis2010a]_ formalism presented by
   [vanHaasteren2012a]_. Also, there is a paper by [Petiteau2012a]_, which presents a way
   how to validate data processing methods by means of injection of new sources and then
   trying to recover them.

====
Meanings of some symbols oftenly used in the cited papers and in this project.
====

The original paper is slightly unclear about some of the parameters
The parameter names:

* :math:`\zeta := M^{5/3}/D`
* :math:`\iota` - orbital inclination
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
 .. [vanHaasteren2012a] http://arxiv.org/abs/1210.0584v2
 .. [Petiteau2012a] http://arxiv.org/abs/1210.2396
