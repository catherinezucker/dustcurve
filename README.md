# dustcurve:

Update on progress as of April 18th, 2016: 
- Wrote code to repackage all the data from Odyssey into a new, more convenient file for running the MCMC calculations. Ran the code to produce "testdata.h5," which contains the stellar posteriors and the co intensities for all the stars in a small area on the sky (56 stars total). This data is 19 MB (sorry!) but it was the smallest convenient data package I could produce. The test data is now stored in the root directory. 
- Created a class "PixStars" to hold the data from my hdf5 file, including both the CO intensity array and the stellar posterior array for all the stars in a single pixel on the sky
- Created a test for the PixStars class and ensured it ran successfully from my root directory with nosetests
- Having a separate io module is really unnecessary at the moment due to how I'm setting up my data in the PixStars class. All you need to due to load your data is to send PixStars the strings containing the file name and pixel of interest within that file name, and you have the data loaded in one line.  
- Coded up the model.py function which contains functions for the log_prior, log_likelihood, and log_posterior functions. Also contains all the functions necessary to take the line integrals over the stellar posterior arrays
- Set up runmcmc.py code in my root directory, which will eventually be merged into my model.ipynb code once I work through all my bugs 
- There is currently a bad bug in the code I use to take my line integrals over the probability surfaces

This package uses the stellar posteriors on distance and reddening for individual stars (derived from near-infrared and optical photometry) to refine the distances to carbon monoxide emission features given by a Galactic rotation curve. Within a Bayesian framework, we use affine invariant MCMC to sample from the posterior of distance and gas-to-dust coefficient in each velocity slice. The inputs into the code include:  

1) The joint posterior of distance and reddening for each star. This is a 700x120 array, with 700 reddening bins (from 0 to 7 magnitudes in reddening) and 120 distance modulus bins (from 4-19 in distance modulus). These posteriors are pulled from the work of Green et al. 2014 and Green et al. 2015.  Each bin is filled with a value representing the probability in that bin.  
2) The CO emission intensity for every pixel for every velocity slice in our field of interest (CO cubes are from the Dame et al. 2001 CO survey).  

