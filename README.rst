=============================================
 Pulsar Timing Array - Reduced Order Methods
=============================================

This is my summer project on Pulsar Timing Array data fitting using Reduced Order
Quadrature methods

The structure of the project:
-----------------------------

1. A Data generator (94%)

2. A Reduced Basis (RB) generator (50%)

3. A parameter fitting program:

 a) Using the conventional methods, described by [vanHaasteren2012a]_ and [Ellis2010a]_. (30%)

 b) Using RB methods (0 %)

Meanings of symbols used in the project
---------------------------------------

Meanings of some symbols often used in the cited papers and in this project. The
original paper (by [Ellis2010a]_) was slightly unclear about some of the parameters,
hence the names are enumerated bellow:

* :math:`\zeta := \mathcal{M}^{5/3}/D`
* :math:`\iota` - orbital inclination
* :math:`\psi` - polarisation angle of the GW
* :math:`t` - time
* :math:`\theta, \phi` - spherical polar angles for a vector pointing from Solar System
  Barycentre to the GW source
* :math:`\theta_p, \phi_p` - angular pulsar coordinates
  This is an added parameter by me, which I will probably substitute to something
  better
* :math:`\omega_0` - Initial frequency
* :math:`L` - Distance to the pulsar
* :math:`D` - Distance to the Blackhole binary
* :math:`\mathcal{M}` - the chirp mass of the binary
  .. math:: \mathcal{M}^{5/3} = \frac{m_1 m_2}{(m_1 + m_2)^{1/3}}
* :math:`\omega(t)` - orbital frequency
* :math:`\Phi(t), \Phi_0` - orbital phase
* :math:`s` - signal from the pulsar timing data
* :math:`h \left(\vec{\lambda}\right)` - the signal template depending on the parameters
  :math:`\vec{\lambda}`
* :math:`C` - The covariance matrix, which is used to calculate the multivariate noise
  Gaussian.

 .. [Ellis2010a] http://arxiv.org/abs/1204.4218
 .. [Corbin2008a] http://arxiv.org/abs/1008.1782
 .. [Jaranowski2008a] http://arxiv.org/abs/gr-qc/9804014
 .. [vanHaasteren2012a] http://arxiv.org/abs/1210.0584v2
 .. [Petiteau2012a] http://arxiv.org/abs/1210.2396

.. vim: tw=88:spell:spelllang=en_gb
