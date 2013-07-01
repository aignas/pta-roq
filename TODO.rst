=============================================
 Pulsar Timing Array - Reduced Order Methods
=============================================

List of TODOs
=============

Data generator
--------------

1. Add a way to generate noisy data:

 1. Red noise

 2. GWB

 3. Power Law noise

Conventional data analysis
--------------------------

1. Finish the covariance matrix calculation function (red and power law noise terms are
   missing at the moment). I need to understand what exactly they mean by various
   parameters in the spectrum formulae and how to store them best in the pulsar data
   structure. Also I need an algorithm, which would convert these frequency domain
   formulae into the time domain.

2. Write an algorithm for the usual data analysis procedure. This should be MCMC with
   maximizing the likelihood via brute force method.

Generation of my reduced basis
------------------------------

1. Implementing ways of selecting points in the parameter space:
 
 1. Traverse the whole parameter space and eliminate the points, which are in the RB.

 2. Traverse the whole parameter space during the first run and then use the properties
    of the norm to remove most of the points which do not lie in a circle defined by the
    current error.

    The current question which is not clear for me before implementing this procedure is
    how to know where to centre the circle (or a some multidimensional spheroid) and
    also why can we do that?

 3. Implement a way to construct the RB by randomly selecting M vectors all the time and
    finding the reduced basis in this set.
   
    The second method should be better when the set of RB is sufficiently large, but it
    is still not as precise. And, because this is the off-line stage, we do not care
    about the execution time that much, but we do care about the precision.

2. How to calculate the projections efficiently, so that we do not have lots of memory
   occupied?

   Current problem is that we need to evaluate lots of products. Theoretically, we could
   store everything, i.e. a grid of parameter space, the training set and projections,
   but this doesn't smell right.

   If our parameter space is large N-dim space, then this method is bound to have memory
   problems.

3. I talked with P. Canizares and the thing, which she suggested, was to use the RB and
   EIM methods on "each" pulsar, because otherwise, the method might not be valid as
   the waveforms depend on the pulsar positions and a single set of RB might not work.

   However, this is still better, than doing the product the hard way. Although, I need
   to get any numbers before stating it with confidence.

   But first I will try the old method where I have the inner product well defined.

Use the Empirical Interpolation Method with the Reduced basis
-------------------------------------------------------------

1. EIM is the reduction of points in the time domain. 

   The method here might be different in a sense, that instead of time series, I would
   have a time/pulsar series, which makes this sort of a different type of problem.
   Also there might be some reasons why I need to change the method of interpolation in
   some way.

2. If Priscilla is correct, then I have several time series for different pulsars, I
   will need to use the EIM and RB for each time series for each pulsar.

3. After doing it I will be able (hopefully) to weed out the Covariance matrix.

Test the above code for specific noise realisations
---------------------------------------------------

1. Test several methods:

 1. Use a simple MCMC algorithm for finding the sources.

 2. Check how the method presented by [Petiteau2012a]_ scales up.

Think about the ROQ for the :math:`\mathcal{F}_e` and :math:`\mathcal{F}_p` statistics   in three different regimes
-------------------------------------------------------------------------------------------------------------------

1. When the Gravitational Wave Background (GWB) is known

2. When the GWB slope is known (find the actual references on Arxiv.org for this)

3. For a general GWB

Check the model for multiple sources
------------------------------------

References
----------

 .. [Petiteau2012a] http://arxiv.org/abs/1210.2396
 .. [Canizares2013a] http://arxiv.org/abs/1304.0462

.. vim: tw=88:spell:spelllang=en_gb
