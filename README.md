**MCContrails** is a stochastic 0D population balance model for contrails. It represents the particle population with a set of *superparticles* and solves droplet growth, freezing, ice crystal growth, and coagulation.


# Compiling

The code is compiled with `make` in the main directory. This creates the file `MCContrails`. If recompiling, `make -B` is safest.


# Executing

Desired initial particle properties, environmental conditions, and simulation settings should be input by amending the files `simulation.in`, `population.in`, `species.in`, and `environment.in`.

Due to the volatility of growth rates in the early plume, it is likely that the dry volume fraction in an interval will try to go above 1 or below 0 at some point. For this reason, a tolerance on this property can be chosen where $1 < f_{dry} < 1+tolerance$ results in $f_{dry} = 1$ and similarly for below 0. The tolerance should be as small as possible.

The program is executed with `./MCContrails`.

Output is in the output directory. The columns of the particle size distribution (`psd.out`) are:
1. Time (s)
2. Particle radius at centre of interval (m)
3. Droplet number density $n$ (m-3)
4. Droplet $dn/dlogr$ (m-3)
5. Crystal number density $n$ (m-3)
6. Crystal $dn/dlogr$ (m-3)

The columns of the environmental variables (`environment.out`) are:
1. Time (s)
2. Temperature (K)
3. Vapour pressure (Pa)
4. Saturation vapour pressure with respect to liquid water (Pa)
5. Saturation vapour pressure with respect to ice (Pa)


# Notes on model dynamics

Plume dynamics currently follow [Kärcher et al. (2015)](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2015JD023491). It can be amended in the `Environment` class as desired.

The droplet growth equation follows Pruppacher and Klett / Seinfeld and Pandis. Freezing and ice crystal growth is adapted from [Kärcher et al. (2015)](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2015JD023491).

Superparticles are initialised with equal fractions of the total number density. A number of superparticles proportional to an input species' number density are initialised with the properties of that species. Radii are randomly assigned to superparticles using a lognormal radius distribution as weights.


# Authors

Code developed and written by Jack Bartlett.