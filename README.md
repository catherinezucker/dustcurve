# dustcurve:
This package uses the stellar posteriors on distance and reddening for individual stars (derived from near-infrared and optical photometry) to refine the distances to carbon monoxide emission features given by a Galactic rotation curve. Within a Bayesian framework, we use parallel tempering MCMC to sample from a 12-dimensional posterior space, consisting of a distance to each velocity slice of a spectral cube. The required data includes:

1) The joint posterior of distance and reddening for each star. This is a 700x120 array, with 700 reddening bins (from 0 to 7 magnitudes in reddening) and 120 distance modulus bins (from 4-19 in distance modulus). These posteriors are pulled from the work of Green et al. 2014 and Green et al. 2015.  Each bin is filled with a value representing the probability in that bin.  
2) The CO emission intensity for every pixel for every velocity slice in our field of interest (CO cubes are from the Dame et al. 2001 CO survey).  

## Files:

- model.ipynb : A jupyter notebook containing details of our statistical modeling, and the application of our analysis to simulated data
- tutorial. ipynb: A jupyter notebook containing sample analysis on real data, and walking through the required steps to run the MCMC code yourself
- simulated_data.h5: the required data file needed to run model.ipynb. Model.ipynb creates this file as part of its script, but if you download this directly, you only need to run the part of the notebook containing the MCMC code
- pixclass.py : a class created to hold our data and metadata, called a PixStars object. In our case, it holds all the stellar posterior arrays, stellar coordinates, and CO intensity information for a single pixel 
- io.py: A module that creates an instance of our class, a PixStar, and returns all the required arguments to run our parallel tempering sampler (PTSampler from the emcee class). It's primary function, fetch_args, accepts a filename, an array containing the lower and upper bounds of our uniform distance prior, and a desired gas-to-dust coefficient. 
- hputils.py & repackage_data.py: packages used to convert the raw data into a more MCMC-friendly format. These modules were used to create our 
