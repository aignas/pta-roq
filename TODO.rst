=============================================
 Pulsar Timing Array - Reduced Order Methods
=============================================

My TODO list:
-------------

1. Data generator:

 a) add a way to generate noisy data:

  1. Red noise

  2. GWB

  3. Power Law noise

2. Conventional data analysis

 a) Finish the covariance matrix calculation function (red and power law noise terms are
    missing at the moment). I need to understand what exactly they mean by them and how
    to store the parameters needed for it the best. Probably it should be the pulsar
    class.

 b) Write an algorithm for the usual data analysis procedure. This should be MCMC with
    maximizing the likelihood via brute force method.

3. Generation of my reduced basis:

 a) Implement two ways of selecting points in the parameter space:
  
  1. Traverse the whole parameter space and eliminate the points, which are in the RB.

  2. Implement a way to construct the RB by randomly selecting M vectors all the time and
     finding the reduced basis in this set.
    
    The second method should be better when the set of RB is sufficiently large.

 b) How to calculate the projections efficiently, so that we do not have lots of memory
    occupied?

4. Use the Empirical Interpolation Method with the Reduced basis.

 a) read the paper by [Canizares2013a]_ again on EIM usage in ROM modelling and apply
    the method to my problem.

5. Test the above code for specific noise realisations and search for the gravitational
   waves in the mock data.

 a) Test several methods:
 
  1. Use a simple MCMC algorithm for finding the sources.

  2. Check how the method presented by [Petiteau2012a]_ scales up.

 b) Generate IPTA challenge like data to do this.

6. Think about the ROQ for the :math:`\mathcal{F}_e` and :math:`\mathcal{F}_p` statistics
   in three different regime:

 a) When the Gravitational Wave Background (GWB) is known

 b) When the GWB slope is known (find the actual references on Arxiv.org for this)

 c) For a general GWB

7. Check the model for multiple sources

 .. [Petiteau2012a] http://arxiv.org/abs/1210.2396
 .. [Canizares2013a] http://arxiv.org/abs/1304.0462

.. vim: tw=88:spell:spelllang=en_gb
