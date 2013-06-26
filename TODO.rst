=============================================
 Pulsar Timing Array - Reduced Order Methods
=============================================

My TODO list:
-------------

1. Data generator:

 a) add red noise

 b) Add a stochastic gravitational wave background

 c) Use CUDA or C++ to speed things up?

2. Conventional data analysis

 a) Implement a method to calculate the covariance matrix

 b) Write an algorithm for the usual data analysis algorithm

3. Generation of my reduced basis:

 a) Implement the inner product for finding the reduced basis

 b) Implement a way to construct the RB by traversing the whole parameter space and
 eliminating the points, which constitute the RB.

 c) Implement a way to construct the RB by randomly selecting M vectors all the time and
 finding the reduced basis in this set.

 The second method should be better when the set of RB is sufficiently large.

4. Use the Empyrical Interpolation Method with the Reduced basis.

 a) read the paper by [Canizares2013a]_ again on EIM usage in ROM modelling and apply
 the method to my problem.

5. Test the above code for specific noise realisations and search for the gravitational
   waves in the mock data.

 a) Use a simple MCMC algorithm for finding the sources.

 b) Check how the method presented by [Petiteau2012a]_ scales up.

6. Think about the ROQ for the :math:`\mathcal{F}_e` and :math:`\mathcal{F}_p` statistics
   in three different regime:

 a) When the Gravitational Wave Background (GWB) is known

 b) When the GWB slope is known (find the actual references on Arxiv.org for this)

 c) For a general GWB

7. Check the model for multiple sources

 .. [Petiteau2012a] http://arxiv.org/abs/1210.2396
 .. [Canizares2013a] http://arxiv.org/abs/1304.0462

.. vim: tw=88:spell:spelllang=en_gb
