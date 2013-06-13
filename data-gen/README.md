A simple data generator
=====

This is a simple data generator using the Ellis et al. paper [1]. You need Cython in
order to run the code. Also, reffer to the paper by Corbin and Cornish [2] for some
detail on the formulae derivations. There is also the Jaranowski et al. [3] paper for
explanation of intrinsic and extrinsic parameter classification.

The original paper is slightly unclear about some of the parameters
The parameter names:
  * $\zeta$ - $M^{5/3}/D$
  * $\cos(\iota)$ - orbital inclination (? is it the cosine)
  * $\psi$ - polarisation angle of the GW
  * $t$ - time
  * $\theta, \phi$ - spherical polar angles for a vector pointing from Solar System
     Baricentre to the GW source
  * $\theta_p, \phi_p$ - angular pulsar coordinates
    This is an added parameter by me, which I will probably substitute to something
    better
  * $\omega_0$ - Initial frequency
  * $L$ - Distance to the pulsar
  * $D$ - Distance to the Blackhole binary
  * $M$ - the chirp mass of the binary
    $M^{5/3} = \frac{m_1 m_2}{(m_1 + m_2)^{1/3}}$
  * $\omega(t)$ - orbital frequency
  * $\Phi(t), \Phi_0$ - orbital phase

 [1] - http://arxiv.org/abs/1204.4218
 [2] - http://arxiv.org/abs/1008.1782
 [3] - http://arxiv.org/abs/gr-qc/9804014
