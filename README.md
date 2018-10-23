# Shallow Water Solitons #

## Problem Description ##

This is a nice example of the n = 1 Rossby wave soliton found in the equatorial beta-plane model given by Boyd (1980). It demonstrates the non-linear, equatorial shallow water equations in Dedalus on a non-trivial test problem. More details can be found in the notes file. 

## How to run ##

1. Run the script `rossby.py`. On my 4-core 2018 chromebook laptop, it takes around 10 minutes to run. E.g.
```
$ mpirun -np 4 python3 rossby.py
```
2. Run the script `plot_2d_series.py`. This will create a new directory `frames`, which will contain a bunch of `png` images showing the height field phi and the x and y velocities u and v, respectively.

## Suggested exercises ##

The Rossby soliton should propagate westward (the negative x direction). However, if you look closely, you'll see a small height perturbation propagate eastward. This is probably a Kelvin wave excited because discretized initial conditions aren't exactly those of the analytic solution or because the analytic solution itself is only a leading order approximation to the true non-linear dynamics. One could test which effect (limited resolution vs limited fidelity of initial conditions) is more important by performing a spatial convergence test. That is, raise and lower nx and ny and see how the Kelvin* wave

The timestep is fixed to `dt = 0.025`, which is a little less than 10% of Boyd's computation of the wave speed at leading order, c = -1/3. 

* I'm not 1oo% sure this is really a Kelvin wave. But it probably is. They go eastward. 

## References ##
Boyd, J. P. "Equatorial Solitary Waves. Part I: Rossby Solitons", Journal of Physical Oceanography (1980), v 10, p 1699

