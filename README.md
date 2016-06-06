# dustcurve:
This package uses the stellar posteriors on distance and reddening for individual stars (derived from near-infrared and optical photometry) to refine the distances to carbon monoxide emission features given by a Galactic rotation curve. Within a Bayesian framework, we use parallel tempering MCMC to sample from a 12-dimensional posterior space, consisting of a distance to each velocity slice of a spectral cube. The required data includes:

1) The joint posterior of distance and reddening for each star. This is a 700x120 array, with 700 reddening bins (from 0 to 7 magnitudes in reddening) and 120 distance modulus bins (from 4-19 in distance modulus). These posteriors are pulled from the work of Green et al. 2014 and Green et al. 2015.  Each bin is filled with a value representing the probability in that bin.  
2) The CO emission intensity for every pixel for every velocity slice in our field of interest (CO cubes are from the Dame et al. 2001 CO survey).  

## Installation
You can install this package by cloning the repository via the terminal ("git clone https://github.com/catherinezucker/dustcurve.git") and then running "python setup.py install" within the root directory. However, you need to first install the python package h5py before you install dustcurve. H5py is dependent on a c library which can only be installed through "conda install h5py" and not through setup.py. If you do not first conda install h5py, you will get import h5py errors when attempting to run files like model.ipynb

Moreover, there are two python packages required to run repackage_data.py that are not included in setup.py: healpy and spectral-cube. Both of these packages are notoriously difficult to install. Repackage_data.py and hputils.py are the only scripts that are dependent on these packages. As the sole purpose of repackage_data.py is to convert the raw data into a more MCMC friendly format, there is no need for a general user to run repackage_data.py, as you would need access to the private Odyssey server with the raw data files. In the future, we plan to make all the repackaged data files available to the public, but for now setup.py only installs the packages required to run all the MCMC code and related files. 

## Files:

- model.ipynb : A jupyter notebook containing details of our statistical modeling, and the application of our analysis to simulated data
- tutorial. ipynb: A jupyter notebook containing sample analysis on real data, and walking through the required steps to run the MCMC code yourself
- simulated_data.h5: the required data file needed to run model.ipynb. Model.ipynb creates this file as part of its script, but if you download this directly, you only need to run the part of the notebook containing the MCMC code
- model.py and modelunc.py: the module contains all the necessary functions to run our MCMC code, including a log_posterior function, a log_likelihood function, a log_prior function, and all the functions that these depend on. See note below on the difference between the two modules
- pixclass.py : a class created to hold our data and metadata, called a PixStars object. In our case, it holds all the stellar posterior arrays, stellar coordinates, and CO intensity information for a single pixel 
- io.py: A module that creates an instance of our class, a PixStar, and returns all the required arguments to run our parallel tempering sampler (PTSampler from the emcee class). It's primary function, fetch_args, accepts a filename, an array containing the lower and upper bounds of our uniform distance prior, and a desired gas-to-dust coefficient. 
- hputils.py & repackage_data.py: modules used to convert the raw data into a more MCMC-friendly format. These modules were used to create the data used in tutorial.ipynb, 89996.h5, which is all the data for the center of the Cepheus cloud. These files are provided for reference, and are not meant to be run by a general user. Repackage_data.py in particular contains user specific directory paths that won't work unless on a specific machine. See our note on setup.py below for more.  
- Cepheus_hires.fits: spectral cube data used in our tutorial.ipynb file. Not necessary to run the MCMC code unless you want to use the script repackage_data.py
- Plot_posterior.py: plots the best fit reddening profiles overtop the stacked stellar posterior surface for the given pixel




