=============================================
 Pulsar Timing Array - Reduced Order Methods
=============================================

This is my summer project on Pulsar Timming Array data fitting using Reduced Order
Quadrature methods

The structure of the project:
-----------------------------

1. This is a simple data generator using the Ellis et al. paper [Ellis2010a]_. You need
   Cython in order to run the code. Also, reffer to the paper by Corbin and Cornish
   [Corbin2008a]_ for some detail on the formulae derivations. There is also the
   Jaranowski et al. [Jaranowski2008a]_ paper for explanation of intrinsic and extrinsic
   parameter classification.

   TODO:

 a) add red noise

 b) Add a stochastic gravitational wave background

2. Analysis of the usual problem, so that proper comparisons can be done.

   Notes on the original method:

   There is also a speed-up version of the [Ellis2010a]_ formalism presented by
   [vanHaasteren2012a]_. Also, there is a paper by [Petiteau2012a]_, which presents a way
   how to validate data processing methods by means of injection of new sources and then
   trying to recover them.

3. Generation of the Reduced Basis for the :math:`n^{T} \hat{s}`, where :math:`\hat{s} =
   C^{-1} s`.

4. Implement the covariance matrix to calculate the actual values of :math:`\hat{s}`

5. Test the above code for specific noise realisations and search for the gravitational
   waves in the mock data.

6. Think about the ROQ for the :math:`mathcal{F}_e` and :math:`mathcal{F}_p` statistics
   in three different regime:

 a) When the Gravitational Wave Background (GWB) is known

 b) When the GWB slope is known (find the actual references on Arxiv.org for this)

 c) For a general GWB

7. Check the model for multiple sources

Meanings of symbols
-------------------

Meanings of some symbols oftenly used in the cited papers and in this project.
The original paper (by [Ellis2010a]_) was slightly unclear about some of the parameters,
hence the names are enumerated bellow:

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
* :math:`s` - signal from the pulsar timing data
* :math:`h \left(\vec{\lambda}\right)` - the signal template depending on the parameters
  :math:`\vec{\lambda}`
* :math:`C` - The covariance matrix, which is used to calculate the multivariate noise
  gaussian.

 .. [Ellis2010a] http://arxiv.org/abs/1204.4218
 .. [Corbin2008a] http://arxiv.org/abs/1008.1782
 .. [Jaranowski2008a] http://arxiv.org/abs/gr-qc/9804014
 .. [vanHaasteren2012a] http://arxiv.org/abs/1210.0584v2
 .. [Petiteau2012a] http://arxiv.org/abs/1210.2396
