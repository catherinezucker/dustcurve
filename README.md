# dustcurve:
This package uses the stellar posteriors on distance and reddening for individual stars (derived from near-infrared and optical photometry) to refine the distances to carbon monoxide emission features given by a Galactic rotation curve. Within a Bayesian framework, we use parallel tempering MCMC to sample from a 12-dimensional posterior space, consisting of a distance to each velocity slice of a spectral cube. The required data includes:

1) The joint posterior of distance and reddening for each star. This is a 700x120 array, with 700 reddening bins (from 0 to 7 magnitudes in reddening) and 120 distance modulus bins (from 4-19 in distance modulus). These posteriors are pulled from the work of Green et al. 2014 and Green et al. 2015.  Each bin is filled with a value representing the probability in that bin.  
2) The CO emission intensity for every pixel for every velocity slice in our field of interest (CO cubes are from the Dame et al. 2001 CO survey).  

## Files:

- model.ipynb : A jupyter notebook containing details of our statistical modeling, and the application of our analysis to simulated data
- tutorial. ipynb: A jupyter notebook containing sample analysis on real data, and walking through the required steps to run the MCMC code yourself
